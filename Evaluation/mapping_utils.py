#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import shutil
import subprocess as sbp
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.cbook import boxplot_stats

from Utils.utils import call_log, time_stamp

# --------------------------- FASTQ discovery ---------------------------

# Cleaned outputs from fastp (priority)
R1_PATTERNS_CLEANED = [
    "*_R1.cleaned.gz",
    "*_1.cleaned.gz",
    "*_R1.cleaned.fastq.gz",
    "*_1.cleaned.fastq.gz",
]

# Common raw patterns (fallback)
R1_PATTERNS_RAW = [
    "*_1.fastq", "*_1.fq", "*_1.fastq.gz", "*_1.fq.gz",
    "*_R1.fastq", "*_R1.fq", "*_R1.fastq.gz", "*_R1.fq.gz",
]


def _discover_pairs_in_dir(wd: Path) -> List[Dict]:
    """
    Return list of {sample, r1, r2|None}, preferring fastp-cleaned names,
    otherwise falling back to raw naming patterns.
    'sample' is the basename before _R1/_1 suffix + extension.
    """
    wd = Path(wd)
    pairs: List[Dict] = []

    # Try cleaned patterns first, then raw patterns
    seen = set()
    for patterns in (R1_PATTERNS_CLEANED, R1_PATTERNS_RAW):
        for pat in patterns:
            for r1p in sorted(wd.glob(pat)):
                name = r1p.name
                if "_R1." in name:
                    base, ext = name.split("_R1.", 1)
                    r2n = f"{base}_R2.{ext}"
                elif "_1." in name:
                    base, ext = name.split("_1.", 1)
                    r2n = f"{base}_2.{ext}"
                else:
                    # doesn't match paired-end naming; skip
                    continue

                r2p = wd / r2n
                sample = base
                key = (sample, r1p, r2p)
                if key in seen:
                    continue
                seen.add(key)

                pairs.append(
                    dict(
                        sample=sample,
                        r1=r1p,
                        r2=r2p if r2p.exists() else None,
                    )
                )

        # if we found anything with this tier of patterns, stop; don't fall back
        if pairs:
            break

    # If still nothing, consider single-end fastqs as singles
    if not pairs:
        all_fqs = (
            sorted(wd.glob("*.fastq"))
            + sorted(wd.glob("*.fq"))
            + sorted(wd.glob("*.fastq.gz"))
            + sorted(wd.glob("*.fq.gz"))
        )
        for p in all_fqs:
            pairs.append(dict(sample=p.stem, r1=p, r2=None))

    return pairs


# --------------------------- Kraken2 helpers ---------------------------

def _append_log_line(log_dir: Path, step: str, msg: str) -> Path:
    """
    Append a human-readable line to <log_dir>/<step>.log (creating dirs as needed).
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{step}.log"
    with open(log_path, "a", encoding="utf-8") as fh:
        fh.write(msg.rstrip() + "\n")
    return log_path


def _run_kraken2(
    r1: Path,
    r2: Optional[Path],
    db: Path,
    threads: int,
    out_report: Path,
    out_assign: Optional[Path] = None,
    log_dir: Optional[Path] = None,
    label: Optional[str] = None,
    gzip_input: bool = True,
) -> None:
    """
    Execute the built-in kraken2 command with the provided inputs.
    """
    out_report.parent.mkdir(parents=True, exist_ok=True)
    if out_assign:
        out_assign.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "kraken2",
        "--db", str(Path(db).resolve()),
        "--threads", str(threads),
        "--use-names",
        "--report", str(Path(out_report).resolve()),
    ]
    if gzip_input:
        cmd.append("--gzip-compressed")
    if out_assign:
        cmd += ["--output", str(Path(out_assign).resolve())]
    if r2:
        cmd.append("--paired")
    cmd.append(str(r1))
    if r2:
        cmd.append(str(r2))

    cmd = " ".join(a for a in cmd)
    time_stamp(f"[contam] {cmd}")
    stdout = None
    stderr = None
    log_handle = None
    if log_dir:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / "kraken2.log"
        log_handle = open(log_path, "ab")
        tag = label or "kraken2"
        header = f"\n[{tag}] {cmd}\n".encode()
        log_handle.write(header)
        stdout = log_handle
        stderr = sbp.STDOUT

    try:
        sbp.run(cmd, shell=True, check=True, stdout=stdout, stderr=stderr)
    finally:
        if log_handle:
            log_handle.flush()
            log_handle.close()


def _summarize_kreport_species(kreport: Path) -> pd.DataFrame:
    """
    Parse kreport -> species-level rows only (rank code starts with 'S').
    Returns DataFrame columns: [name, pct, reads_clade, reads_direct, taxid, rank]
    """
    rows = []
    with open(kreport) as fh:
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 6:  # pct, reads_clade, reads_direct, rank, taxid, name
                continue
            try:
                pct = float(parts[0])
                rc = int(parts[1])
                rd = int(parts[2])
            except (TypeError, ValueError):
                continue
            rank = parts[3].strip()
            taxid = parts[4].strip()
            name = parts[5].strip()
            if rank and rank[0].upper() == "S":  # species-level
                rows.append(dict(name=name, pct=pct, reads_clade=rc, reads_direct=rd, taxid=taxid, rank=rank))
    return pd.DataFrame(rows)


def _kreport_total_reads(kreport: Path) -> Optional[int]:
    """
    Return the total number of sequences Kraken2 saw (root reads_clade from kreport).
    """
    try:
        with open(kreport) as fh:
            for ln in fh:
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    total = int(parts[1])
                except ValueError:
                    continue
                return total
    except Exception as e:
        time_stamp(f"[contam] failed to read kreport {kreport}: {e}")
    return None


def _annotate_box_fliers(ax, x_values: Iterable, y_values: Iterable, labels: Iterable, rotation: int = 45) -> None:
    """
    Label boxplot fliers on the given axes.
    x_values / y_values / labels must be aligned and of equal length.
    """
    x_arr = np.asarray(list(x_values))
    y_arr = np.asarray(list(y_values), dtype=float)
    labels_arr = np.asarray(list(labels))

    if x_arr.size == 0 or y_arr.size == 0 or labels_arr.size == 0:
        return

    texts = []
    for x_val in np.unique(x_arr):
        mask = x_arr == x_val
        if not np.any(mask):
            continue
        x_vals = x_arr[mask]
        y_vals = y_arr[mask]
        labs = labels_arr[mask]
        bstats = boxplot_stats(y_vals, labels=None)[0]
        fliers = bstats.get("fliers", [])
        if not fliers:
            continue
        fliers_set = set(fliers)
        for xv, yv, lab in zip(x_vals, y_vals, labs):
            if yv in fliers_set:
                texts.append(ax.text(xv, yv, lab, rotation=rotation, ha="right", va="center", fontsize=7))

    if texts:
        try:
            adjust_text(
                texts,
                arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5),
                force_points=(0.1, 0.2),
                force_text=(0.05, 0.2),
                expand_points=(1.05, 1.1),
                lim=150,
            )
        except Exception:
            pass


# --------------------------- BAM / FASTQ counters ---------------------------

def _count_mapped_primary_reads(bam: Path, threads: int) -> int:
    """
    Count mapped, primary, non-supplementary reads (exclude unmapped/secondary/supplementary).
    """
    cmd = f"samtools view -c -F 0x904 -@ {threads} {bam}"
    try:
        out = sbp.check_output(cmd, shell=True, text=True).strip()
        return int(out)
    except Exception:
        call_log(bam.parent, "samtools_count_mapped_primary_reads", cmd)
        return 0


def _count_unmapped_primary_reads(bam: Path, threads: int) -> int:
    """
    Count unmapped primary reads (FLAG 0x4, excluding secondary/supplementary).
    """
    cmd = f"samtools view -c -f 0x4 -F 0x900 -@ {threads} {bam}"
    try:
        out = sbp.check_output(cmd, shell=True, text=True).strip()
        return int(out)
    except Exception:
        call_log(bam.parent, "samtools_count_unmapped_primary_reads", cmd)
        return 0


def _count_reads_in_fastq(path: Path) -> int:
    """Number of reads = lines/4 (handles .gz, uses pigz for speed)."""
    if str(path).endswith(".gz"):
        pigz = shutil.which("pigz")
        if not pigz:
            time_stamp("[contam] pigz not found on PATH; install pigz in the environment.")
            return 0
        cmd = f"{pigz} -cd {path} | wc -l"
    else:
        cmd = f"wc -l {path} | awk '{{print $1}}'"
    try:
        out = sbp.check_output(cmd, shell=True, text=True).strip()
        nlines = int(out.split()[0])
        return nlines // 4
    except Exception as e:
        time_stamp(f"[contam] failed to count reads in {path}: {e}")
        return 0


def _prepare_kraken_inputs(sample: str, r1: Path, r2: Optional[Path]):
    """
    Return (r1_path, r2_path, total_reads, gzip_input) for Kraken2.
    Returns None if no reads are present.
    """
    r1_reads = _count_reads_in_fastq(r1)
    r2_reads = _count_reads_in_fastq(r2) if r2 else 0
    if r2:
        if r1_reads != r2_reads:
            time_stamp(f"[contam] warning: read count mismatch for {sample} (R1={r1_reads}, R2={r2_reads}); using max")
        total_pairs = max(r1_reads, r2_reads)
    else:
        total_pairs = r1_reads

    if total_pairs == 0:
        return None
    gzip_input = str(r1).endswith(".gz") and (not r2 or str(r2).endswith(".gz"))

    return r1, r2, total_pairs, gzip_input
