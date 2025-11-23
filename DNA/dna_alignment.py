import sys
import subprocess as sbp # enable activate and terminate program
from Utils.utils import call_log, time_stamp
from typing import List, Dict, Set, Tuple
from pathlib import Path

R1_PATTERNS = [
    "*_R1.cleaned.gz",
    "*_1.cleaned.gz",
    # If your suffix ever changes, add variants like:
    "*_R1.cleaned.fastq.gz",
    "*_1.cleaned.fastq.gz",
]

class DNAAlignment: 
    def __init__(self, args):
        self.args = args
        self.align_outdir = Path(args.align_outdir)
        self.align_outdir.mkdir(parents=True, exist_ok=True)
        self.align_log = self.align_outdir / "alignment.log"
        self.align_log.touch(exist_ok=True)
        self.seq_platform = args.seq_platform
        self.mark_duplicates = getattr(args, "mark_duplicates", True)
        # Allow explicit override; otherwise derive from platform with a safe default
        self.optical_distance = getattr(args, "optical_distance", None)
        if self.optical_distance is None:
            if self.seq_platform in ['NovaSeq']:
                self.optical_distance = 2500
                # ref: http://www.htslib.org/algorithms/duplicate.html
            elif self.seq_platform in ['HiSeq']:
                self.optical_distance = 100
            elif self.seq_platform.lower() == 'other':
                self.optical_distance = 100
                print('Sequencing platform set to "other"; defaulting optical distance to 100. Override via optical_distance if needed.')
            else:
                self.optical_distance = 100
                print('Sequencing platform is not recorded; defaulting optical distance to HiSeq (100).')
        
        self.processed_outdir = Path(args.fastp_outdir)

        self.threads = args.threads
        self.ref_fasta = args.ref_fasta
        self.population = args.population
        self.mem = self._normalize_mem(self.args.mem)
        self.skip_samples = set(getattr(args, "skip_align_samples", set()) or [])
        
    @staticmethod    
    def discover_pairs(wd: Path) -> List[Tuple[str, str, str, str]]:
        """
        Find paired-end files in wd.

        Returns list of (basename, style, r1_name, r2_name)
        style = "R" if filenames use _R1.cleaned/_R2.cleaned, else "1" for _1.cleaned/_2.cleaned
        """
        seen = set()
        pairs: List[Tuple[str, str, str, str]] = []
        for pat in R1_PATTERNS:
            for r1p in wd.glob(pat):
                fname = r1p.name
                if "_R1.cleaned" in fname:
                    base = fname.split("_R1.", 1)[0]
                    ext = fname.split("_R1.", 1)[1]
                    r2 = f"{base}_R2.{ext}"
                    style = "R"
                elif "_1.cleaned" in fname:
                    base = fname.split("_1.", 1)[0]
                    ext = fname.split("_1.", 1)[1]
                    r2 = f"{base}_2.{ext}"
                    style = "1"
                else:
                    continue
                r2p = wd / r2
                if r2p.exists():
                    key = (base, r1p.name, r2p.name)
                    if key not in seen:
                        seen.add(key)
                        pairs.append((base, style, r1p.name, r2p.name))
                else:
                    print(f"[info] Skip {fname}: missing pair {r2}", file=sys.stderr)
        return pairs
    
    @staticmethod
    def _normalize_mem(mem: str) -> str:
        s = str(mem).upper().replace(" ", "")
        return s.replace("GB", "G").replace("MB", "M")
        
    def run(self):
        #time_stamp("Starting DNA alignment", self.args.logfile)
        pairs = self.discover_pairs(self.processed_outdir)
        if not pairs:
            print(f"[warn] No cleaned read pairs found in {self.processed_outdir}", file=sys.stderr)
            return
        
        for base, style, r1, r2 in pairs:
            if base in self.skip_samples:
                time_stamp(f"[align] skip sample {base} (listed in skip_alignment_samples)")
                continue
            fastq1 = str(self.processed_outdir / r1)
            fastq2 = str(self.processed_outdir / r2)
            readgroup = base.split('.')[0]
            align_dir_str = str(self.align_outdir)

            cmd_parts = [
                "bwa-mem2 mem -t {0} -R '@RG\\tID:{7}\\tPL:illumina\\tLB:library\\tSM:{8}' {1} {2} {3}".format(
                    self.threads, self.ref_fasta, fastq1, fastq2, self.mem, self.align_outdir, base, self.population, readgroup
                )
            ]
            if self.mark_duplicates:
                cmd_parts.append("samtools collate -@ {0} -O -u -".format(self.threads))
                cmd_parts.append("samtools fixmate -m -@ {0} - -".format(self.threads))
            cmd_parts.append("samtools sort -m {0} -@ {1}".format(self.mem, self.threads))
            if self.mark_duplicates:
                cmd_parts.append(
                    "samtools markdup -r -s -d {0} -@ {1} -f {2}/{3}.dupstat - -".format(
                        self.optical_distance, self.threads, align_dir_str, base
                    )
                )
            cmd_parts.append(
                "samtools view -b -f 2 -@ {0} -F 2048 -o {1}/{2}.bam".format(
                    self.threads, align_dir_str, base
                )
            )
            # build single-line pipeline to avoid mangling with clean_cmd
            cmd = " | ".join(cmd_parts)
            cmd = f"set -o pipefail; {cmd}"
            # fixmate: fills in mate coordinates and insert size fields.
            # it remove duplicate during the process
            # -f 2:  require read mapped in proper pair (0x2)
            # -F 2048: remove supplementary alignment (0x800)

            # do not clean the command; keep pipes intact
            bam_path = self.align_outdir / f"{base}.bam"
            log_path = self.align_outdir / "bam_index.log"

            try:
                with open(self.align_log, "a") as log_fh:
                    log_fh.write(time_stamp(f"Running alignment for {base}") + "\n")
                    log_fh.write(cmd + "\n")
                    log_fh.flush()
                    sbp.run(
                        cmd,
                        stdout=log_fh,
                        stderr=sbp.STDOUT,
                        shell=True,
                        executable="/bin/bash",
                        check=True
                    )
            except sbp.CalledProcessError:
                call_log(self.align_outdir, 'alignment', cmd)
                sys.exit(1) #terminate the process and output the error messages

            # Create BAM index (.bai)
            index_cmd = "samtools index -@ {t} -b '{bam}' >> '{log}' 2>&1".format(
                t=self.threads,
                bam=str(bam_path),
                log=str(log_path),
            )
            try:
                sbp.run(index_cmd, shell=True, check=True, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)
            except sbp.CalledProcessError:
                call_log(self.align_outdir, 'bam_index', index_cmd)
                sys.exit(1)
        
        
