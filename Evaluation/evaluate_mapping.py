#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import shutil
import subprocess as sbp
from pathlib import Path
from typing import Dict, List, Optional, Iterable, Tuple
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.cbook import boxplot_stats
from adjustText import adjust_text

from Utils.utils import time_stamp, clean_cmd, call_log

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

    cmd = clean_cmd(" ".join(a for a in cmd))
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


def _annotate_box_fliers(ax, x_values: Iterable, y_values: Iterable, labels: Iterable, rotation: int = 45, auto_adjust: bool = False) -> None:
    """
    Label boxplot fliers on the given axes.
    x_values / y_values / labels must be aligned and of equal length.
    If auto_adjust is True and adjustText is available, labels are nudged to reduce overlap.
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

    if auto_adjust and texts:
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
    cmd = clean_cmd(f"samtools view -c -F 0x904 -@ {threads} {bam}")
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
    cmd = clean_cmd(f"samtools view -c -f 0x4 -F 0x900 -@ {threads} {bam}")
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


# =========================== Evaluator ===========================

class MappingEvaluator:
    """
    Steps:
      0: FastQC + MultiQC
      1: Kraken2 on all reads; aggregate species Top-5 across samples (boxplot)
      2: Mapping yield ratio = unique (primary, non-supplementary) BAM reads / fastp reads; one boxplot
    """

    def _setup_local_env(self) -> None:
        """
        Keep temp/cache writes inside the configured output directory.
        """
        tmp_dir = self.outdir / ".tmp"
        mpl_dir = self.outdir / ".cache" / "matplotlib"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        mpl_dir.mkdir(parents=True, exist_ok=True)
        for key in ("TMPDIR", "TEMP", "TMP"):
            os.environ[key] = str(tmp_dir)
        os.environ["MPLCONFIGDIR"] = str(mpl_dir)
        os.environ.setdefault("MPLBACKEND", "Agg")
        self.tmp_dir = tmp_dir
        self.mpl_config_dir = mpl_dir

    def __init__(self, args):
        self.threads   = int(getattr(args, "threads", 1))
        self.outdir    = Path(getattr(args, "outdir", ".")).resolve()
        self.input_dir = Path(getattr(args, "input_dir", ".")).resolve()

        self._setup_local_env()

        self.fastp_outdir = Path(getattr(args, "fastp_outdir", self.input_dir)).resolve()
        self.align_dir    = Path(getattr(args, "align_outdir", self.outdir)).resolve()

        self.gff = Path(getattr(args, "gff", "")).expanduser().resolve() if getattr(args, "gff", None) else None

        raw_tag = getattr(args, "tag_outliers", True)
        if isinstance(raw_tag, str):
            self.tag_outliers = raw_tag.strip().lower() not in {"0", "false", "no", "off", ""}
        else:
            self.tag_outliers = bool(raw_tag)

        default_kraken_db = "/kraken_db"
        db_val = getattr(args, "kraken_db", default_kraken_db)
        if isinstance(db_val, str):
            db_val = (db_val or default_kraken_db).strip()
        else:
            db_val = str(db_val).strip() if db_val else default_kraken_db
        self.kraken_db = db_val

        self.report_dir = self.outdir / "evaluation"
        self.report_dir.mkdir(parents=True, exist_ok=True)
        (self.report_dir / "logs").mkdir(parents=True, exist_ok=True)

    # ---------------- step 0/1: QC ----------------
    def evaluate_qc(self, fastq_stage: str = "auto") -> None:
        """
        Run FastQC + MultiQC on FASTQs from a specified stage.
        fastq_stage:
          - 'raw': only consider the original input directory
          - 'processed': only consider fastp-cleaned outputs
          - 'auto' (default): prefer processed, fall back to raw if nothing found
        """
        stage = (fastq_stage or "auto").strip().lower()
        choices = {"raw", "processed", "auto"}
        if stage not in choices:
            stage = "auto"

        stage_dirs = []
        if stage in {"processed", "auto"}:
            stage_dirs.append(("processed", self.fastp_outdir))
        if stage in {"raw", "auto"}:
            stage_dirs.append(("raw", self.input_dir))

        used_stage = None
        pairs = []
        for label, wd in stage_dirs:
            pairs = _discover_pairs_in_dir(wd)
            if pairs:
                used_stage = label
                break

        fastqs: List[Path] = []
        for d in pairs:
            fastqs.append(d["r1"])
            if d["r2"]:
                fastqs.append(d["r2"])
        if not fastqs:
            time_stamp(f"[qc] no FASTQs found for stage '{stage}'; skip")
            return

        stage_dir = used_stage or stage
        fq_out = self.report_dir / "fastqc" / stage_dir
        mq_out = self.report_dir / "multiqc" / stage_dir
        fq_out.mkdir(parents=True, exist_ok=True)
        mq_out.mkdir(parents=True, exist_ok=True)

        files_str = " ".join(str(p) for p in fastqs)
        cmd1 = clean_cmd(f"fastqc -t {self.threads} -o {fq_out} {files_str}")
        time_stamp(f"[qc:{stage_dir}] {cmd1}")
        try:
            sbp.run(cmd1, shell=True, check=True)
        except sbp.CalledProcessError:
            call_log(self.report_dir / "logs", f"fastqc_{stage_dir}", cmd1)

        cmd2 = clean_cmd(f"multiqc -o {mq_out} {fq_out}")
        time_stamp(f"[qc:{stage_dir}] {cmd2}")
        try:
            sbp.run(cmd2, shell=True, check=True)
        except sbp.CalledProcessError:
            call_log(self.report_dir / "logs", f"multiqc_{stage_dir}", cmd2)

    def evaluate_contamination(self) -> Optional[pd.DataFrame]:
        if not self.kraken_db:
            time_stamp("[contam] skip: kraken_db not set")
            return None

        pairs = _discover_pairs_in_dir(self.fastp_outdir)
        if not pairs:
            pairs = _discover_pairs_in_dir(self.input_dir)
        if not pairs:
            time_stamp("[contam] no FASTQs for Kraken2")
            return None

        log_dir = self.report_dir / "logs"
        db_path = Path(self.kraken_db).expanduser()
        if not db_path.exists():
            msg = f"[contam] kraken_db path does not exist: {db_path}"
            time_stamp(msg)
            _append_log_line(log_dir, "kraken2", msg)
            return None
        db_path = db_path.resolve()
        required = ["hash.k2d", "opts.k2d", "taxo.k2d"]
        missing = [db_path / fname for fname in required if not (db_path / fname).exists()]
        if missing:
            missing_list = ", ".join(p.name for p in missing)
            msg = f"[contam] kraken_db missing required files ({missing_list}); skip Kraken2"
            time_stamp(msg)
            _append_log_line(log_dir, "kraken2", msg)
            return None

        out_dir = self.report_dir / "kraken"
        out_dir.mkdir(parents=True, exist_ok=True)

        all_rows = []

        for rec in pairs:
            sample, r1, r2 = rec["sample"], rec["r1"], rec["r2"]

            prepared = _prepare_kraken_inputs(sample, r1, r2)
            if not prepared:
                msg = f"[contam] sample {sample} has 0 reads; skip Kraken2"
                time_stamp(msg)
                _append_log_line(log_dir, "kraken2", msg)
                continue
            r1_for_kraken, r2_for_kraken, sample_total_records, gzip_input = prepared

            krep = out_dir / f"{sample}.kreport"
            kout = out_dir / f"{sample}.kout"
            try:
                _run_kraken2(
                    r1=r1_for_kraken, r2=r2_for_kraken,
                    db=db_path,
                    threads=self.threads,
                    out_report=krep,
                    out_assign=kout,
                    log_dir=log_dir,
                    label=sample,
                    gzip_input=gzip_input,
                )
                sample_total_records = _kreport_total_reads(krep) or sample_total_records
                if kout.exists():
                    try:
                        sbp.run(["pigz", "-f", str(kout)], check=True, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)
                    except Exception as e:
                        time_stamp(f"[contam] pigz failed for {kout}: {e}")
            except sbp.CalledProcessError as e:
                failed_cmd = e.cmd if isinstance(e.cmd, str) else "kraken2"
                call_log(
                    log_dir,
                    "kraken2",
                    f"[sample {sample}] {failed_cmd} (exit {e.returncode})",
                )
                continue

            df = _summarize_kreport_species(krep)
            if df is not None and not df.empty:
                df["sample"] = sample
                df["total_records"] = sample_total_records
                all_rows.append(df)

        if not all_rows:
            time_stamp("[contam] no species rows parsed")
            return None

        big = pd.concat(all_rows, ignore_index=True)

        big = big.rename(columns={"reads_clade": "n_records"})
        big["name_freq_pct"] = big["n_records"] * 100.0 / big["total_records"].clip(lower=1)

        med = big.groupby("name")["name_freq_pct"].median()
        top5_names = list(med.sort_values(ascending=False).head(5).index)
        plot_top = big[big["name"].isin(top5_names)]
        if plot_top.empty:
            time_stamp("[contam] skip plot: no rows after Top-5 filter")
            return big

        plt.figure(figsize=(7, 5))
        ax = sns.boxplot(
            data=plot_top,
            x="name",
            y="name_freq_pct",
            width=0.4,
        )

        x_values = plot_top["name"].map({name: i for i, name in enumerate(top5_names)})
        y_values = plot_top["name_freq_pct"]
        labels = plot_top["sample"]
        if self.tag_outliers:
            _annotate_box_fliers(ax, x_values, y_values, labels, rotation=90, auto_adjust=True)

        medians = plot_top.groupby("name")["name_freq_pct"].median()
        for i, name in enumerate(top5_names):
            med = medians.get(name)
            if med is None or pd.isna(med):
                continue
            ax.text(i, med, f"{med:.1f}%", ha="center", va="bottom", fontsize=8, color="red")

        ax.set_xlabel("Species (Top-5 by median frequency across samples)")
        ax.set_ylabel("Percentage of records per sample")
        ax.set_title(
            "Kraken2 contamination — Top-5 species\n"
            "by median per-sample frequency (scientific name)"
        )
        plt.tight_layout()
        out_png = out_dir / "kraken_top5_species_boxplot.png"
        plt.savefig(out_png, dpi=180)
        plt.close()
        time_stamp(f"[contam] wrote {out_png}")

        out_tsv = out_dir / "kraken_species_tidy.tsv"
        big.to_csv(out_tsv, sep="\t", index=False)
        time_stamp(f"[contam] wrote {out_tsv}")

        return big

    # ---------------- step 2: mapping yield (single whisker) ----------------
    def evaluate_alignment(self) -> pd.DataFrame:
        if not self.align_dir.exists():
            time_stamp(f"[kept] align_dir not found: {self.align_dir}")
            return pd.DataFrame()

        bams = sorted(self.align_dir.glob("*.bam")) + sorted(self.align_dir.glob("*.cram"))
        if not bams:
            time_stamp(f"[kept] no BAM/CRAM in {self.align_dir}")
            return pd.DataFrame()

        pairs = _discover_pairs_in_dir(self.fastp_outdir) or _discover_pairs_in_dir(self.input_dir)
        by_sample = {d["sample"]: d for d in pairs}

        rows = []
        for bam in bams:
            sample = bam.stem
            num_mapped = _count_mapped_primary_reads(bam, self.threads)
            num_unmapped = _count_unmapped_primary_reads(bam, self.threads)

            denom_reads = None
            if sample in by_sample:
                r1 = by_sample[sample]["r1"]
                r2 = by_sample[sample]["r2"]
                denom_reads = _count_reads_in_fastq(r1)
                if r2:
                    denom_reads += _count_reads_in_fastq(r2)

            ratio = (100.0 * num_mapped / denom_reads) if denom_reads else None
            unmapped_ratio = (100.0 * num_unmapped / denom_reads) if denom_reads else None
            rows.append(dict(sample=sample,
                             unique_reads_in_bam=num_mapped,
                             unmapped_reads_in_bam=num_unmapped,
                             fastp_reads=denom_reads,
                             ratio_unique_pct=ratio,
                             ratio_unmapped_pct=unmapped_ratio))

        df = pd.DataFrame(rows)
        out_csv = self.report_dir / "unique_reads_in_bam_vs_fastp.csv"
        df.to_csv(out_csv, index=False)
        time_stamp(f"[kept] wrote {out_csv}")

        melt_cols = []
        if df["ratio_unique_pct"].notna().any():
            melt_cols.append("ratio_unique_pct")
        if "ratio_unmapped_pct" in df and df["ratio_unmapped_pct"].notna().any():
            melt_cols.append("ratio_unmapped_pct")

        if melt_cols:
            plot_df = df[["sample"] + melt_cols].dropna().melt(id_vars="sample", var_name="metric", value_name="pct")
            metric_labels = {
                "ratio_unique_pct": "Unique / fastp (%)",
                "ratio_unmapped_pct": "Unmapped / fastp (%)",
            }
            plot_df["metric_label"] = plot_df["metric"].map(metric_labels)

            plt.figure(figsize=(5, 5))
            ax = sns.boxplot(data=plot_df,
                             x="metric_label",
                             y="pct",
                             width=0.4)

            x_map = {lab: i for i, lab in enumerate(plot_df["metric_label"].unique())}
            if self.tag_outliers:
                _annotate_box_fliers(
                    ax,
                    x_values=plot_df["metric_label"].map(x_map),
                    y_values=plot_df["pct"],
                    labels=plot_df["sample"],
                    auto_adjust=True,
                )

            medians = plot_df.groupby("metric_label")["pct"].median()
            for label, x_pos in x_map.items():
                med = medians.get(label)
                if med is None or pd.isna(med):
                    continue
                ax.text(x_pos, med, f"{med:.1f}%", ha="center", va="bottom", fontsize=8, color="red")

            ax.set_xlabel("")
            ax.set_ylabel("Percentage of fastp reads")
            ax.set_title("Mapping summary (unique vs. unmapped)")
            plt.tight_layout()
            out_png = self.report_dir / "mapping_yield_boxplots.png"
            plt.savefig(out_png, dpi=180)
            plt.close()
            time_stamp(f"[kept] wrote {out_png}")

        return df
