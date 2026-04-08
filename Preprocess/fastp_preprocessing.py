#!/usr/bin/env python3
"""
run_fastp.py 

Process paired-end FASTQ files with fastp.

- Detects pairs named either <base>_1 / <base>_2 or <base>_R1 / <base>_R2
  with extensions .fq/.fastq[.gz]
- Optional *front* trimming ONLY:
    * one global integer (e.g., --trim-front 10)
    * or a 2-column file mapping sample_key -> value, where sample_key ends
      with _1/_2 or _R1/_R2, e.g.:
         SRR26245751_1  8
         SRR26245751_2  6
         4202P-HPGonad_R1 12
         4202P-HPGonad_R2 12
- Writes cleaned reads to --outdir as <base>_<R1or1>.cleaned.<ext>
- Generates per-sample HTML/JSON reports

Notes:
- For PE data, adapter trimming is robust via overlap; you can add
  --detect-adapter-pe for ultra-cleaning.
- Avoid per-read cutting options when using --dedup; global front trim (-f/-F) is safe.
"""

import sys
import argparse
from pathlib import Path
from typing import Dict, Tuple, List, Optional
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
# Allow running this script from outside the repo.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from Utils.utils import discover_pairs_by_patterns

R1_PATTERNS = [
    "*_1.fastq", "*_1.fq", "*_1.fastq.gz", "*_1.fq.gz",
    "*_R1.fastq", "*_R1.fq", "*_R1.fastq.gz", "*_R1.fq.gz",
    "*_1_*.fastq", "*_1_*.fq", "*_1_*.fastq.gz", "*_1_*.fq.gz",
    "*_R1_*.fastq", "*_R1_*.fq", "*_R1_*.fastq.gz", "*_R1_*.fq.gz",
]
FASTP_MAX_THREADS = 64
FASTP_PREFERRED_THREADS = 16
FASTP_MIN_PARALLEL_THREADS = 8


def _distribute_fastp_threads(total_threads: int, n_jobs: int) -> List[int]:
    total_threads = max(1, int(total_threads))
    n_jobs = max(1, int(n_jobs))
    n_jobs = min(n_jobs, total_threads)
    base = total_threads // n_jobs
    extra = total_threads % n_jobs
    alloc = [base + (1 if idx < extra else 0) for idx in range(n_jobs)]
    return [max(1, min(FASTP_MAX_THREADS, val)) for val in alloc]


def _choose_fastp_job_count(total_threads: int, n_pairs: int) -> int:
    total_threads = max(1, int(total_threads))
    n_pairs = max(1, int(n_pairs))

    # Prefer more sample-level concurrency and keep per-process thread counts
    # in a range where fastp tends to scale better on large servers.
    preferred_jobs = max(1, (total_threads + FASTP_PREFERRED_THREADS - 1) // FASTP_PREFERRED_THREADS)
    max_jobs_by_floor = max(1, total_threads // FASTP_MIN_PARALLEL_THREADS)
    return max(1, min(n_pairs, preferred_jobs, max_jobs_by_floor))

def load_trim_file(trim_file: str) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Read a 2-column trim file:
        <sample_key> <trim_value>
    and split into R1 and R2 dicts based on *_1 / *_R1 / *_2 / *_R2 suffix.
    Values must be integers (bp).
    """
    r1: Dict[str, int] = {}
    r2: Dict[str, int] = {}
    with open(trim_file) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()  # split on any whitespace
            if len(parts) != 2:
                print(f"[warn] Bad trim line (skip): {line}", file=sys.stderr)
                continue
            key, val = parts
            try:
                iv = int(val)
            except ValueError:
                print(f"[warn] Non-integer trim value (skip): {line}", file=sys.stderr)
                continue
            if key.endswith(("_1", "_R1")):
                r1[key] = iv
            elif key.endswith(("_2", "_R2")):
                r2[key] = iv
            else:
                print(f"[warn] Key has no _1/_2 or _R1/_R2 suffix (skip): {key}", file=sys.stderr)
    return r1, r2

def discover_pairs(wd: Path) -> List[Tuple[str, str, str, str]]:
    return discover_pairs_by_patterns(wd, R1_PATTERNS)

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Batch run fastp on paired-end FASTQs (front-only trimming).",
        formatter_class=argparse.RawTextHelpFormatter
    )
    p.add_argument("-d", "--wd", required=True, help="Working directory containing input FASTQs")
    p.add_argument("-o", "--outdir", required=True, help="Output directory (must exist)")
    p.add_argument(
        "-w",
        "--threads",
        type=int,
        default=16,
        help="Threads for fastp (default: 16, max: 64)"
    )
    p.add_argument("-l", "--min-length", type=int, default=30, help="--length_required (default: 30)")
    p.add_argument("-q", "--qualified-phred", type=int, default=20, help="Qualified base threshold Q (default: 20)")
    p.add_argument("-e", "--average-qual", type=int, default=0, help="Drop reads with avg qual < this (0=off)")
    p.add_argument(
        "--lib-type",
        choices=["DNA", "RNA"],
        required=True,
        help="Library type tunes defaults. If RNA, enables --trim_poly_x"
    )
    p.add_argument(
        "--trim-front",
        help="Global integer (front trim) OR path to 2-col file mapping sample_key -> value",
        default=""
    )
    p.add_argument("--detect-adapter-pe", action="store_true",
                   help="Enable adapter detection for PE (adds --detect_adapter_for_pe)")
    p.add_argument(
        "--dedup",
        action="store_true",
        help="Enable fastp deduplication (--dedup). Not allowed with --lib-type RNA."
    )
    p.add_argument("--report-dir", default="", help="Directory for HTML/JSON reports (default: outdir)")
    return p

def run_fastp(
    wd: Path,
    outdir: Path,
    threads: int,
    min_length: int,
    qualified_phred: int,
    average_qual: int,
    lib_type: str,
    trim_front: str,
    detect_adapter_pe: bool,
    dedup: bool,
    report_dir: Optional[Path]
) -> None:

    if not wd.exists():
        raise FileNotFoundError(f"Working dir not found: {wd}")
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    if report_dir and not report_dir.exists():
        raise FileNotFoundError(f"Report dir not found: {report_dir}")

    if threads < 1:
        print(f"[warn] fastp threads must be >= 1; using 1 (requested {threads})", file=sys.stderr)
        threads = 1

    if lib_type == "RNA" and dedup:
        raise ValueError("--dedup is not allowed for RNA libraries. Remove --dedup or use --lib-type DNA.")
    if dedup:
        print(
            "[warn] fastp dedup is enabled. This mode is much heavier on large DNA FASTQs and can stall on some datasets; "
            "leave fastp_dedup=false unless you explicitly need fastp duplicate removal.",
            file=sys.stderr,
            flush=True,
        )

    # parse front trimming
    per_r1: Dict[str, int] = {}
    per_r2: Dict[str, int] = {}
    global_front: Optional[int] = None
    if trim_front:
        tf = Path(trim_front)
        if tf.exists():
            per_r1, per_r2 = load_trim_file(str(tf))
        else:
            try:
                global_front = int(trim_front)
            except ValueError:
                print(f"[warn] --trim-front must be int or file. Ignored: {trim_front}", file=sys.stderr)

    pairs = discover_pairs(wd)
    if not pairs:
        print(f"[info] No paired FASTQs found in {wd}", file=sys.stderr)
        return

    max_parallel_jobs = _choose_fastp_job_count(threads, len(pairs))
    thread_plan = _distribute_fastp_threads(threads, max_parallel_jobs)
    print(
        f"[info] fastp scheduling: total_threads={threads}, processes={max_parallel_jobs}, "
        f"per_process_threads={thread_plan}, target_threads_per_process={FASTP_PREFERRED_THREADS}",
        flush=True,
    )

    jobs = []
    for pair_idx, (base, style, r1, r2) in enumerate(pairs):
        # preserve original style in output names
        rtag1, rtag2 = ("R1", "R2") if style == "R" else ("1", "2")

        out1 = outdir / f"{base}_{rtag1}.cleaned.fastq.gz"
        out2 = outdir / f"{base}_{rtag2}.cleaned.fastq.gz"

        # sample-specific front trims (prefer exact key; fall back to alt style; else global)
        key_r1 = f"{base}_{rtag1}"
        key_r2 = f"{base}_{rtag2}"
        alt_r1 = f"{base}_{'1' if style=='R' else 'R1'}"
        alt_r2 = f"{base}_{'2' if style=='R' else 'R2'}"

        trim_r1 = per_r1.get(key_r1, per_r1.get(alt_r1, global_front))
        trim_r2 = per_r2.get(key_r2, per_r2.get(alt_r2, global_front))

        repdir = report_dir if report_dir else outdir
        html = repdir / f"{base}.fastp.html"
        jsn  = repdir / f"{base}.fastp.json"

        assigned_threads = thread_plan[pair_idx % len(thread_plan)]
        cmd = [
            "fastp",
            "-w", str(assigned_threads),
            "-i", str(wd / r1),
            "-I", str(wd / r2),
            "-o", str(out1),
            "-O", str(out2),
            "-l", str(min_length),
            "-q", str(qualified_phred),
            "-h", str(html),
            "-j", str(jsn)
        ]

        if average_qual and average_qual > 0:
            cmd += ["-e", str(average_qual)]

        if detect_adapter_pe:
            cmd += ["--detect_adapter_for_pe"]  # v1.0.2 supports this long flag

        if dedup:
            cmd += ["--dedup"]  # equivalent to -D

        # Library-type sensible defaults
        if lib_type == "RNA":
            # PolyX (e.g., polyA/T) tails common in mRNA-Seq; keep polyG auto
            cmd += ["--trim_poly_x"]

        # Apply *front* trims if set  <-- uses fastp -f/-F
        if isinstance(trim_r1, int):
            cmd += ["-f", str(trim_r1)]
        if isinstance(trim_r2, int):
            cmd += ["-F", str(trim_r2)]

        jobs.append((base, cmd))

    def _run_one(job: Tuple[str, List[str]]) -> str:
        base, cmd = job
        print("[info] Running:", " ".join(cmd), flush=True)
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[error] fastp failed for {base}: {e}", file=sys.stderr)
            raise
        return base

    if max_parallel_jobs == 1:
        for job in jobs:
            _run_one(job)
        return

    with ThreadPoolExecutor(max_workers=max_parallel_jobs) as ex:
        futures = [ex.submit(_run_one, job) for job in jobs]
        for fut in as_completed(futures):
            fut.result()

def main():
    ap = build_argparser()
    a = ap.parse_args()
    if a.lib_type == "RNA" and a.dedup:
        ap.error("--dedup cannot be used with --lib-type RNA.")

    run_fastp(
        wd=Path(a.wd).resolve(),
        outdir=Path(a.outdir).resolve(),
        threads=a.threads,
        min_length=a.min_length,
        qualified_phred=a.qualified_phred,
        average_qual=a.average_qual,
        lib_type=a.lib_type,
        trim_front=a.trim_front,
        detect_adapter_pe=a.detect_adapter_pe,
        dedup=a.dedup,
        report_dir=Path(a.report_dir).resolve() if a.report_dir else None
    )
    print("[info] All samples processed.", flush=True)

if __name__ == "__main__":
    main()
