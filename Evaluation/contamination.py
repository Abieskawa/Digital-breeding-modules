#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess as sbp
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from Utils.utils import (
    _annotate_box_fliers,
    _append_log_line,
    _count_reads_in_fastq,
    _discover_pairs_in_dir,
    _kreport_total_reads,
    call_log,
    time_stamp,
)
from Evaluation.mapping_base import MappingEvalBase


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
    rows = []
    with open(kreport) as fh:
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 6:
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
            if rank and rank[0].upper() == "S":
                rows.append(
                    dict(
                        name=name,
                        pct=pct,
                        reads_clade=rc,
                        reads_direct=rd,
                        taxid=taxid,
                        rank=rank,
                    )
                )
    return pd.DataFrame(rows)


def _prepare_kraken_inputs(sample: str, r1: Path, r2: Optional[Path]):
    r1_reads = _count_reads_in_fastq(r1)
    r2_reads = _count_reads_in_fastq(r2) if r2 else 0
    if r2:
        if r1_reads != r2_reads:
            time_stamp(
                f"[contam] warning: read count mismatch for {sample} (R1={r1_reads}, R2={r2_reads}); using max"
            )
        total_pairs = max(r1_reads, r2_reads)
    else:
        total_pairs = r1_reads

    if total_pairs == 0:
        return None
    gzip_input = str(r1).endswith(".gz") and (not r2 or str(r2).endswith(".gz"))

    return r1, r2, total_pairs, gzip_input


class ContaminationEvaluator(MappingEvalBase):
    """
    Step 1: Kraken2 on all reads; aggregate species Top-5 across samples (boxplot).
    """

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
            _annotate_box_fliers(ax, x_values, y_values, labels, rotation=90)

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
