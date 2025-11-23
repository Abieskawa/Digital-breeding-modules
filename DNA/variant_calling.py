#!/usr/bin/env python3
# DNA/variant_calling.py

import sys
import shutil
import subprocess as sbp
from pathlib import Path
from typing import List, Tuple, Optional

from Utils.utils import clean_cmd, call_log


class VariantCalling:
    """
    Steps:
      1) DeepVariant (Docker) per-sample → per-sample VCF + gVCF
      2) GLnexus joint-call → cohort.merged.vcf.gz (+.tbi)
      3) plink QC → cohort.filtered.vcf.gz (+.tbi)
      4) plink LD-prune → cohort.filtered.ldpruned.vcf.gz (+.tbi)
      5) pigz archive BAMs

    Required args:
        variant_outdir, align_outdir, ref_fasta, threads

    Optional defaults:
        deepvariant_image='google/deepvariant:1.9.0'
        model_type='WGS'  # or 'WES'
        capture_bed=None
        hwe_p='1e-5'
        geno_missing='0.10'
        ld_method='indep'        # or 'indep-pairwise'
        ld_window=50
        ld_step=5
        ld_param=2.0             # VIF (indep) or r^2 (pairwise)
    """

    def __init__(self, args):
        self.variant_outdir = Path(args.variant_outdir).resolve()
        self.align_outdir   = Path(args.align_outdir).resolve()
        self.ref_fasta      = Path(args.ref_fasta).resolve()
        self.threads        = int(args.threads)
        self.skip_samples   = set(getattr(args, 'skip_variant_samples', set()) or getattr(args, 'skip_align_samples', set()) or [])

        self.variant_outdir.mkdir(parents=True, exist_ok=True)

        self.deepvariant_image = getattr(args, 'deepvariant_image', 'google/deepvariant:1.9.0')
        self.model_type        = getattr(args, 'model_type', 'WGS')
        cap                    = getattr(args, 'capture_bed', None)
        self.capture_bed       = Path(cap).resolve() if cap else None

        self.hwe_p        = str(getattr(args, 'hwe_p', '1e-5'))
        self.geno_missing = str(getattr(args, 'geno_missing', '0.10'))

        self.ld_method = str(getattr(args, 'ld_method', 'indep'))
        self.ld_window = int(getattr(args, 'ld_window', 50))
        self.ld_step   = int(getattr(args, 'ld_step', 5))
        self.ld_param  = float(getattr(args, 'ld_param', 2.0))

        if not self.ref_fasta.exists():
            raise FileNotFoundError('ref_fasta not found: {0}'.format(self.ref_fasta))
        if self.capture_bed and not self.capture_bed.exists():
            self.capture_bed = None  # ignore silently

    # ---------------- helpers ----------------

    @staticmethod
    def _discover_bams(wd: Path) -> List[Tuple[str, Path]]:
        return [(p.stem, p.resolve()) for p in sorted(wd.glob('*.bam'))]

    @staticmethod
    def _vcf_paths(out_base: Path) -> Tuple[Path, Path]:
        raw_vcf = out_base.parent / f"{out_base.name}.vcf"
        gz_vcf = raw_vcf.with_suffix('.vcf.gz')
        return raw_vcf, gz_vcf

    @staticmethod
    def _run_quiet(cmd: str, workdir: Optional[Path] = None, step: str = 'cmd',
                   logdir: Optional[Path] = None):
        cmd = clean_cmd(cmd)
        log_handle = None
        stdout = sbp.DEVNULL
        stderr = sbp.DEVNULL

        if logdir:
            log_path_parent = Path(logdir)
            log_path_parent.mkdir(parents=True, exist_ok=True)
            log_path = log_path_parent / f"{step}.log"
            # open in binary to safely capture arbitrary subprocess output
            log_handle = open(log_path, 'wb')
            stdout = log_handle
            stderr = sbp.STDOUT

        try:
            sbp.run(cmd,
                    shell=True,
                    check=True,
                    stdout=stdout,
                    stderr=stderr,
                    cwd=str(workdir) if workdir else None)
        except sbp.CalledProcessError:
            if logdir:
                call_log(logdir, step, cmd)
            sys.exit(1)
        finally:
            if log_handle:
                log_handle.close()

    # ---------------- DeepVariant (per sample) ----------------

    def run_deepvariant_one(self, sample: str, bam_path: Path) -> Tuple[Path, Path]:
        ref_dir = self.ref_fasta.parent
        bam_dir = bam_path.parent
        v_out   = self.variant_outdir / '{0}.vcf.gz'.format(sample)
        g_out   = self.variant_outdir / '{0}.g.vcf.gz'.format(sample)

        cmd = '/opt/deepvariant/bin/run_deepvariant \
               --model_type={4} \
               --ref={5} \
               --reads={1}/{6} \
               --output_vcf={2}/{7} \
               --output_gvcf={2}/{8} \
               --num_shards={9}'.format(
                   ref_dir,                # {0}
                   bam_dir,                # {1}
                   self.variant_outdir,    # {2}
                   self.deepvariant_image, # {3}
                   self.model_type,        # {4}
                   self.ref_fasta,    # {5}
                   bam_path.name,          # {6}
                   v_out.name,             # {7}
                   g_out.name,             # {8}
                   self.threads            # {9}
               )
        self._run_quiet(cmd, step='deepvariant', logdir=self.variant_outdir)

        # index outputs with host tabix
        for out in (v_out, g_out):
            if out.exists() and not out.with_suffix(out.suffix + '.tbi').exists():
                self._run_quiet('tabix -p vcf {0} -@ {1}'.format(out, self.threads),
                                step='tabix_index', logdir=self.variant_outdir)
        return v_out, g_out

    # ---------------- GLnexus (joint-call) ----------------

    def joint_call_glnexus(self) -> Path:
        gvcfs = sorted(self.variant_outdir.glob('*.g.vcf.gz'))
        if not gvcfs:
            raise RuntimeError('No gVCFs found in {0}'.format(self.variant_outdir))

        listfile = self.variant_outdir / 'gvcf.list'
        listfile.write_text('\n'.join(str(p) for p in gvcfs) + '\n')

        merged_vcf = self.variant_outdir / 'cohort.merged.vcf.gz'
        cfg = 'DeepVariantWES' if self.model_type.upper() == 'WES' else 'DeepVariantWGS'
        bed_arg = '--bed {0}'.format(self.capture_bed) if self.capture_bed else ''
        glnexus_db = self.variant_outdir / 'GLnexus.DB'

        # GLnexus refuses to start if the DB dir from a previous attempt remains.
        if glnexus_db.exists():
            if glnexus_db.is_dir():
                shutil.rmtree(glnexus_db)
            else:
                glnexus_db.unlink()

        cmd = 'glnexus_cli --config {0} \
                           --threads {1} {2} \
                           --dir {5} \
                           --list {3} | \
                           bcftools view -Oz -o {4}'.format(
                               cfg, self.threads, bed_arg, listfile, merged_vcf, glnexus_db
                           )
        self._run_quiet(cmd, workdir=self.variant_outdir, step='glnexus', logdir=self.variant_outdir)
        self._run_quiet('tabix -p vcf {0}'.format(merged_vcf),
                        step='tabix_index', logdir=self.variant_outdir)
        return merged_vcf

    # ---------------- plink QC ----------------

    def plink_qc(self, input_vcf: Path, out_prefix: str = 'cohort.filtered') -> Path:
        out_base = self.variant_outdir / out_prefix
        raw_vcf, gz_vcf = self._vcf_paths(out_base)

        cmd = (
            f"plink --vcf {input_vcf} "
            "--vcf-half-call missing "
            "--allow-extra-chr "
            "--snps-only just-acgt "
            "--biallelic-only strict "
            f"--geno {self.geno_missing} "
            f"--hwe {self.hwe_p} midp "
            f"--threads {self.threads} "
            "--recode vcf "
            f"--out {out_base}"
        )
        self._run_quiet(cmd, step='plink_qc', logdir=self.variant_outdir)
        if not raw_vcf.exists():
            raise RuntimeError(f"plink did not produce expected VCF: {raw_vcf}")
        self._run_quiet(f'bgzip -f {raw_vcf}', step='bgzip_qc_vcf', logdir=self.variant_outdir)
        self._run_quiet(f'tabix -p vcf {gz_vcf}',
                        step='tabix_index', logdir=self.variant_outdir)
        return gz_vcf

    # ---------------- plink LD-pruning ----------------

    def ld_prune(self, input_vcf: Path, out_prefix: str = 'cohort.filtered.ldpruned') -> Path:
        prune_prefix = self.variant_outdir / 'ldprune'

        if self.ld_method == 'indep-pairwise':
            indep_flag = f"--indep-pairwise {self.ld_window} {self.ld_step} {self.ld_param}"
        else:
            indep_flag = f"--indep {self.ld_window} {self.ld_step} {self.ld_param}"

        cmd_indep = (
            f"plink --vcf {input_vcf} "
            "--allow-extra-chr "
            f"{indep_flag} "
            f"--threads {self.threads} "
            f"--out {prune_prefix}"
        )
        self._run_quiet(cmd_indep, step='plink_indep', logdir=self.variant_outdir)

        prune_in = prune_prefix.with_suffix('.prune.in')
        if not prune_in.exists():
            raise RuntimeError('plink did not produce prune.in')

        out_base   = self.variant_outdir / out_prefix
        raw_vcf, pruned_vcf = self._vcf_paths(out_base)

        cmd_export = (
            f"plink --vcf {input_vcf} "
            "--allow-extra-chr "
            f"--extract {prune_in} "
            "--recode vcf "
            f"--out {out_base}"
        )
        self._run_quiet(cmd_export, step='plink_export', logdir=self.variant_outdir)
        if not raw_vcf.exists():
            raise RuntimeError(f"plink did not produce expected VCF: {raw_vcf}")
        self._run_quiet(f'bgzip -f {raw_vcf}', step='bgzip_ld_vcf', logdir=self.variant_outdir)
        self._run_quiet(f'tabix -p vcf {pruned_vcf}',
                        step='tabix_index', logdir=self.variant_outdir)
        return pruned_vcf

    # ---------------- pigz archive BAMs ----------------

    def pigz_archive_bams(self, archive_name: str = 'bam_archive.tar.gz') -> Path:
        bams = [p.name for p in sorted(self.align_outdir.glob('*.bam'))]
        tar_path = self.align_outdir / archive_name
        if not bams:
            return tar_path
        cmd = 'cd {0} && tar -I pigz -cf {1} {2}'.format(
            self.align_outdir, tar_path.name, ' '.join(bams)
        )
        self._run_quiet(cmd, workdir=self.align_outdir, step='pigz_bam', logdir=self.align_outdir)
        return tar_path

    # ---------------- full pipeline ----------------

    def run_all(self) -> None:
        bams = self._discover_bams(self.align_outdir)
        if self.skip_samples:
            bams = [(s, p) for s, p in bams if s not in self.skip_samples]
        if not bams:
            raise RuntimeError('No BAM files found in {0}'.format(self.align_outdir))

        for sample, bam in bams:
            self.run_deepvariant_one(sample, bam)

        merged   = self.joint_call_glnexus()
        filtered = self.plink_qc(merged, out_prefix='cohort.filtered')
        _pruned  = self.ld_prune(filtered, out_prefix='cohort.filtered.ldpruned')
        _archive = self.pigz_archive_bams()
