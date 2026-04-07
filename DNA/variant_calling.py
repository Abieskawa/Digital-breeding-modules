#!/usr/bin/env python3
# DNA/variant_calling.py

import shutil
import shlex
import subprocess as sbp
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple

from Utils.utils import run_quiet, time_stamp, coerce_bool


class VariantCalling:
    """
    Steps:
      1) DeepVariant (Docker) per-sample → per-sample VCF + gVCF
      2) GLnexus joint-call → cohort.merged.vcf.gz (+.tbi)
      3) plink2 QC (+MAF) → cohort.filtered.vcf.gz (+.tbi)
      4) plink2 LD-prune → cohort.filtered.ldpruned.vcf.gz (+.tbi)
      5) pigz archive BAMs

    Required args:
        variant_outdir, align_outdir, ref_fasta, threads

    Optional defaults:
        model_type='WGS'
        capture_bed=None
        hwe_p='1e-5'
        geno_missing='0.10'
        maf=0.05
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

        self.deepvariant_use_gpu = coerce_bool(getattr(args, 'deepvariant_use_gpu', None), False)
        raw_mode = str(getattr(args, 'deepvariant_mode', '') or '').strip().lower()
        self.deepvariant_mode = raw_mode or ('gpu' if self.deepvariant_use_gpu else 'cpu')
        if self.deepvariant_mode not in {'cpu', 'gpu', 'auto'}:
            raise ValueError(f"Unsupported deepvariant_mode: {self.deepvariant_mode}")
        self.deepvariant_docker_mode = coerce_bool(getattr(args, 'deepvariant_docker_mode', None), False)
        self.deepvariant_gpu_image = str(
            getattr(args, 'deepvariant_gpu_image', 'google/deepvariant:1.10.0-gpu') or 'google/deepvariant:1.10.0-gpu'
        ).strip()
        self.deepvariant_cpu_image = str(
            getattr(args, 'deepvariant_cpu_image', 'google/deepvariant:1.10.0-beta') or 'google/deepvariant:1.10.0-beta'
        ).strip()
        self.gpu_devices = str(getattr(args, 'gpu_devices', '0') or '0').strip()
        self.step3_max_concurrent_samples = max(1, int(getattr(args, 'step3_max_concurrent_samples', 1) or 1))
        self.step3_gpu_jobs = max(1, int(getattr(args, 'step3_gpu_jobs', 1) or 1))
        self.step3_cpu_jobs = max(1, int(getattr(args, 'step3_cpu_jobs', 1) or 1))
        self.deepvariant_num_shards = max(1, int(getattr(args, 'deepvariant_num_shards', self.threads) or self.threads))
        self._deepvariant_use_gpu_supported = None
        self._deepvariant_backend = None
        self.model_type        = 'WGS'
        self.deepvariant_extra_args = str(getattr(args, 'deepvariant_extra_args', '') or '').strip()
        cap                    = getattr(args, 'capture_bed', None)
        self.capture_bed       = Path(cap).resolve() if cap else None

        self.hwe_p        = str(getattr(args, 'hwe_p', '1e-5'))
        self.geno_missing = str(getattr(args, 'geno_missing', '0.10'))
        maf_raw = getattr(args, 'maf', None)
        try:
            self.maf = float(maf_raw) if maf_raw not in (None, "", False) else 0.05
        except (TypeError, ValueError):
            self.maf = 0.05
        if not (0.0 < self.maf < 0.5):
            time_stamp(f"[variant] Invalid maf={self.maf}; using default 0.05")
            self.maf = 0.05

        ld_prune_raw = getattr(args, 'enable_ld_pruning', None)
        if ld_prune_raw is None:
            ld_prune_raw = getattr(args, 'ld_prune', None)
        if ld_prune_raw is None:
            ld_prune_raw = getattr(args, 'ld_pruning', None)
        if ld_prune_raw is not None and not coerce_bool(ld_prune_raw, True):
            time_stamp("[variant] LD pruning disabled in config; forcing enabled.")
        self.enable_ld_pruning = True

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

    def _deepvariant_supports_use_gpu(self, deepvariant_cmd: str) -> bool:
        if self._deepvariant_use_gpu_supported is not None:
            return self._deepvariant_use_gpu_supported
        help_cmd = f"{deepvariant_cmd} --help"
        try:
            out = sbp.check_output(help_cmd, shell=True, text=True, stderr=sbp.STDOUT)
        except Exception as exc:
            time_stamp(f"[deepvariant] unable to verify --use_gpu support ({exc}); running without GPU flag.")
            self._deepvariant_use_gpu_supported = False
            return False
        self._deepvariant_use_gpu_supported = "--use_gpu" in out
        if not self._deepvariant_use_gpu_supported:
            time_stamp("[deepvariant] --use_gpu not supported by this DeepVariant build; running without GPU flag.")
        return self._deepvariant_use_gpu_supported

    def _local_deepvariant_bin(self) -> str:
        jemalloc_path = Path("/usr/lib/x86_64-linux-gnu/libjemalloc.so.2")
        deepvariant_bin = "/opt/deepvariant/bin/run_deepvariant"
        if jemalloc_path.exists():
            deepvariant_bin = f"LD_PRELOAD={jemalloc_path} {deepvariant_bin}"
        return deepvariant_bin

    @staticmethod
    def _probe_command(cmd: str) -> bool:
        try:
            sbp.run(
                cmd,
                shell=True,
                check=True,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
            )
            return True
        except Exception:
            return False

    def _probe_docker_backend(self, use_gpu: bool) -> bool:
        if not self.deepvariant_docker_mode or shutil.which("docker") is None:
            return False
        gpus_arg = f"--gpus device={shlex.quote(self.gpu_devices)} " if use_gpu else ""
        image = self.deepvariant_gpu_image if use_gpu else self.deepvariant_cpu_image
        cmd = (
            f"docker run --rm {gpus_arg}{shlex.quote(image)} "
            "/opt/deepvariant/bin/run_deepvariant --help"
        )
        return self._probe_command(cmd)

    def _probe_local_backend(self, use_gpu: bool) -> bool:
        deepvariant_bin = self._local_deepvariant_bin()
        if use_gpu:
            return self._deepvariant_supports_use_gpu(deepvariant_bin)
        return self._probe_command(f"{deepvariant_bin} --help")

    def _select_deepvariant_backend(self) -> str:
        if self._deepvariant_backend is not None:
            return self._deepvariant_backend

        if self.deepvariant_mode == 'cpu':
            candidates = ['docker_cpu', 'local_cpu'] if self.deepvariant_docker_mode else ['local_cpu']
        elif self.deepvariant_mode == 'gpu':
            candidates = ['docker_gpu', 'local_gpu', 'docker_cpu', 'local_cpu'] if self.deepvariant_docker_mode else ['local_gpu', 'local_cpu']
        else:
            candidates = ['docker_gpu', 'local_gpu', 'docker_cpu', 'local_cpu'] if self.deepvariant_docker_mode else ['local_gpu', 'local_cpu']

        for backend in candidates:
            if backend == 'docker_gpu' and self._probe_docker_backend(True):
                self._deepvariant_backend = backend
                break
            if backend == 'docker_cpu' and self._probe_docker_backend(False):
                self._deepvariant_backend = backend
                break
            if backend == 'local_gpu' and self._probe_local_backend(True):
                self._deepvariant_backend = backend
                break
            if backend == 'local_cpu' and self._probe_local_backend(False):
                self._deepvariant_backend = backend
                break

        if self._deepvariant_backend is None:
            raise RuntimeError("No usable DeepVariant backend was found (checked Docker/local CPU/GPU paths).")

        time_stamp(
            f"[deepvariant] Selected backend={self._deepvariant_backend} "
            f"(mode={self.deepvariant_mode}, docker_mode={self.deepvariant_docker_mode})"
        )
        return self._deepvariant_backend

    def _docker_mounts(self, bam_path: Path) -> str:
        mounts = [
            (self.ref_fasta.parent.resolve(), "ro"),
            (bam_path.parent.resolve(), "ro"),
            (self.variant_outdir.resolve(), "rw"),
        ]
        if self.capture_bed:
            mounts.append((self.capture_bed.parent.resolve(), "ro"))

        seen = set()
        parts = []
        for path, mode in mounts:
            key = (str(path), mode)
            if key in seen:
                continue
            seen.add(key)
            parts.append(f"-v {shlex.quote(f'{path}:{path}:{mode}')}")
        return " ".join(parts)

    def _build_deepvariant_cmd(self, sample: str, bam_path: Path, backend: str) -> Tuple[str, Path, Path]:
        v_out = self.variant_outdir / '{0}.vcf.gz'.format(sample)
        g_out = self.variant_outdir / '{0}.g.vcf.gz'.format(sample)
        regions_arg = f"--regions={shlex.quote(str(self.capture_bed))}" if self.capture_bed else ""
        extra_args = f"{self.deepvariant_extra_args} " if self.deepvariant_extra_args else ""

        common_args = (
            f"--model_type={shlex.quote(str(self.model_type))} "
            f"--ref={shlex.quote(str(self.ref_fasta))} "
            f"--reads={shlex.quote(str(bam_path))} "
            f"{regions_arg} "
            f"{extra_args}"
            f"--output_vcf={shlex.quote(str(v_out))} "
            f"--output_gvcf={shlex.quote(str(g_out))} "
            f"--num_shards={self.deepvariant_num_shards}"
        )

        if backend == 'docker_gpu':
            cmd = (
                f"docker run --rm --gpus device={shlex.quote(self.gpu_devices)} "
                f"{self._docker_mounts(bam_path)} "
                f"{shlex.quote(self.deepvariant_gpu_image)} "
                f"/opt/deepvariant/bin/run_deepvariant {common_args}"
            )
            return cmd, v_out, g_out

        if backend == 'docker_cpu':
            cmd = (
                f"docker run --rm "
                f"{self._docker_mounts(bam_path)} "
                f"{shlex.quote(self.deepvariant_cpu_image)} "
                f"/opt/deepvariant/bin/run_deepvariant {common_args}"
            )
            return cmd, v_out, g_out

        env_prefix = "TF_CPP_MIN_LOG_LEVEL=2 PYTHONWARNINGS=ignore"
        gpu_arg = ""
        deepvariant_bin = self._local_deepvariant_bin()
        if backend == 'local_gpu' and self._deepvariant_supports_use_gpu(deepvariant_bin):
            gpu_arg = "--use_gpu "
        cmd = f"{env_prefix} {deepvariant_bin} {gpu_arg}{common_args}"
        return cmd, v_out, g_out

    # ---------------- DeepVariant (per sample) ----------------

    def run_deepvariant_one(self, sample: str, bam_path: Path, backend: str) -> Tuple[Path, Path]:
        cmd, v_out, g_out = self._build_deepvariant_cmd(sample, bam_path, backend)
        run_quiet(cmd, step=f'deepvariant_{sample}', logdir=self.variant_outdir)

        # index outputs with host tabix
        for out in (v_out, g_out):
            if out.exists() and not out.with_suffix(out.suffix + '.tbi').exists():
                label = 'gvcf' if '.g.vcf.gz' in out.name else 'vcf'
                run_quiet(
                    'tabix -p vcf {0} -@ {1}'.format(out, self.threads),
                    step=f'tabix_index_{sample}_{label}',
                    logdir=self.variant_outdir,
                )
        return v_out, g_out

    # ---------------- GLnexus (joint-call) ----------------

    def joint_call_glnexus(self) -> Path:
        gvcfs = sorted(self.variant_outdir.glob('*.g.vcf.gz'))
        if not gvcfs:
            raise RuntimeError('No gVCFs found in {0}'.format(self.variant_outdir))

        listfile = self.variant_outdir / 'gvcf.list'
        listfile.write_text('\n'.join(str(p) for p in gvcfs) + '\n')

        merged_vcf = self.variant_outdir / 'cohort.merged.vcf.gz'
        cfg = 'DeepVariantWGS'
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
        run_quiet(cmd, workdir=self.variant_outdir, step='glnexus', logdir=self.variant_outdir)
        run_quiet('tabix -p vcf {0}'.format(merged_vcf),
                  step='tabix_index', logdir=self.variant_outdir)
        return merged_vcf

    # ---------------- plink2 QC ----------------

    def plink_qc(self, input_vcf: Path, out_prefix: str = 'cohort.filtered') -> Path:
        out_base = self.variant_outdir / out_prefix
        raw_vcf, gz_vcf = self._vcf_paths(out_base)

        cmd = (
            f"plink2 --vcf {input_vcf} "
            "--vcf-half-call missing "
            "--allow-extra-chr "
            "--snps-only just-acgt "
            "--biallelic-only strict "
            f"--geno {self.geno_missing} "
            f"--maf {self.maf} "
            f"--hwe {self.hwe_p} midp "
            f"--threads {self.threads} "
            "--export vcf "
            f"--out {out_base}"
        )
        run_quiet(cmd, step='plink_qc', logdir=self.variant_outdir)
        if not raw_vcf.exists():
            raise RuntimeError(f"plink2 did not produce expected VCF: {raw_vcf}")
        run_quiet(f'bgzip -f {raw_vcf}', step='bgzip_qc_vcf', logdir=self.variant_outdir)
        run_quiet(f'tabix -p vcf {gz_vcf}',
                  step='tabix_index', logdir=self.variant_outdir)
        return gz_vcf

    # ---------------- plink2 LD-pruning ----------------

    def ld_prune(self, input_vcf: Path, out_prefix: str = 'cohort.filtered.ldpruned') -> Path:
        prune_prefix = self.variant_outdir / 'ldprune'

        if self.ld_method == 'indep-pairwise':
            indep_flag = f"--indep-pairwise {self.ld_window} {self.ld_step} {self.ld_param}"
        else:
            indep_flag = f"--indep {self.ld_window} {self.ld_step} {self.ld_param}"

        cmd_indep = (
            f"plink2 --vcf {input_vcf} "
            "--allow-extra-chr "
            f"{indep_flag} "
            f"--threads {self.threads} "
            f"--out {prune_prefix}"
        )
        run_quiet(cmd_indep, step='plink_indep', logdir=self.variant_outdir)

        prune_in = prune_prefix.with_suffix('.prune.in')
        if not prune_in.exists():
            raise RuntimeError('plink2 did not produce prune.in')

        out_base   = self.variant_outdir / out_prefix
        raw_vcf, pruned_vcf = self._vcf_paths(out_base)

        cmd_export = (
            f"plink2 --vcf {input_vcf} "
            "--allow-extra-chr "
            f"--extract {prune_in} "
            "--export vcf "
            f"--out {out_base}"
        )
        run_quiet(cmd_export, step='plink_export', logdir=self.variant_outdir)
        if not raw_vcf.exists():
            raise RuntimeError(f"plink2 did not produce expected VCF: {raw_vcf}")
        run_quiet(f'bgzip -f {raw_vcf}', step='bgzip_ld_vcf', logdir=self.variant_outdir)
        run_quiet(f'tabix -p vcf {pruned_vcf}',
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
        run_quiet(cmd, workdir=self.align_outdir, step='pigz_bam', logdir=self.align_outdir)
        return tar_path

    # ---------------- full pipeline ----------------

    def run_all(self) -> None:
        bams = self._discover_bams(self.align_outdir)
        if self.skip_samples:
            bams = [(s, p) for s, p in bams if s not in self.skip_samples]
        if not bams:
            raise RuntimeError('No BAM files found in {0}'.format(self.align_outdir))

        backend = self._select_deepvariant_backend()
        if 'gpu' in backend:
            deepvariant_workers = min(len(bams), self.step3_max_concurrent_samples, self.step3_gpu_jobs)
        else:
            deepvariant_workers = min(len(bams), self.step3_max_concurrent_samples, self.step3_cpu_jobs)
        deepvariant_workers = max(1, deepvariant_workers)
        time_stamp(
            f"[deepvariant] Running {len(bams)} sample(s) with backend={backend}, "
            f"max_concurrent_samples={deepvariant_workers}, num_shards={self.deepvariant_num_shards}"
        )

        if deepvariant_workers == 1:
            for sample, bam in bams:
                self.run_deepvariant_one(sample, bam, backend)
        else:
            with ThreadPoolExecutor(max_workers=deepvariant_workers) as ex:
                futures = [ex.submit(self.run_deepvariant_one, sample, bam, backend) for sample, bam in bams]
                for fut in as_completed(futures):
                    fut.result()

        merged   = self.joint_call_glnexus()
        filtered = self.plink_qc(merged, out_prefix='cohort.filtered')
        _pruned  = self.ld_prune(filtered, out_prefix='cohort.filtered.ldpruned')
        _archive = self.pigz_archive_bams()
