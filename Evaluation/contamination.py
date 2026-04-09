#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import subprocess as sbp
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from Utils.utils import (
    _annotate_box_fliers,
    _append_log_line,
    _count_reads_in_fastq,
    _discover_pairs_in_dir,
    _kreport_total_reads,
    time_stamp,
)
from Evaluation.mapping_base import MappingEvalBase


def _normalize_species_hint(raw: Optional[str]) -> str:
    if not raw:
        return ""
    return " ".join(str(raw).replace("_", " ").split()).strip().lower()


def _resolve_kraken_db_path(raw_db: str) -> Path:
    raw_db = str(raw_db or "").strip()
    if not raw_db:
        return Path("")

    candidates = [Path(raw_db).expanduser()]
    raw_path = Path(raw_db)
    if not raw_path.is_absolute():
        candidates.append(Path("/") / raw_db)

    seen = set()
    for cand in candidates:
        cand_str = str(cand)
        if cand_str in seen:
            continue
        seen.add(cand_str)
        if cand.exists():
            return cand.resolve()
    return candidates[0].resolve(strict=False)


def _build_non_target_summaries(
    species_df: pd.DataFrame,
    target_species: Optional[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if species_df is None or species_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    work = species_df.copy()
    work["name_norm"] = (
        work["name"]
        .astype(str)
        .str.lower()
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )
    target_norm = _normalize_species_hint(target_species)
    if target_norm:
        work = work[~work["name_norm"].str.contains(re.escape(target_norm), regex=True, na=False)]
    if work.empty:
        return pd.DataFrame(), pd.DataFrame()

    per_sample = (
        work.sort_values(["sample", "name_freq_pct", "n_records"], ascending=[True, False, False])
        .groupby("sample", as_index=False)
        .first()
        .rename(
            columns={
                "name": "top_non_target_species",
                "name_freq_pct": "top_non_target_pct",
                "n_records": "top_non_target_records",
            }
        )
    )
    per_sample = per_sample[
        ["sample", "top_non_target_species", "top_non_target_pct", "top_non_target_records", "total_records"]
    ]

    cohort = (
        work.groupby("name", as_index=False)
        .agg(
            median_freq_pct=("name_freq_pct", "median"),
            max_freq_pct=("name_freq_pct", "max"),
            samples_detected=("sample", "nunique"),
        )
        .sort_values(["median_freq_pct", "max_freq_pct", "samples_detected"], ascending=[False, False, False])
    )
    return per_sample, cohort


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


def _run_kraken_job(job: Dict[str, object]) -> Tuple[str, Optional[str]]:
    sample = str(job["sample"])
    try:
        _run_kraken2(
            r1=job["r1"],
            r2=job["r2"],
            db=job["db"],
            threads=int(job["threads"]),
            out_report=job["kreport"],
            out_assign=job["kout"],
            log_dir=job.get("log_dir"),
            label=sample,
            gzip_input=bool(job["gzip_input"]),
        )
        kout = job["kout"]
        if isinstance(kout, Path) and kout.exists():
            try:
                sbp.run(["pigz", "-f", str(kout)], check=True, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)
            except Exception as exc:
                return sample, f"[contam] pigz failed for {kout}: {exc}"
        return sample, None
    except sbp.CalledProcessError as exc:
        failed_cmd = exc.cmd if isinstance(exc.cmd, str) else "kraken2"
        return sample, f"[sample {sample}] {failed_cmd} (exit {exc.returncode})"
    except Exception as exc:
        return sample, f"[sample {sample}] {exc}"


class ContaminationEvaluator(MappingEvalBase):
    """
    Step 1: Kraken2 on all reads; aggregate species Top-5 across samples (boxplot).
    """

    def evaluate_contamination(self) -> Optional[pd.DataFrame]:
        pairs = _discover_pairs_in_dir(self.fastp_outdir)
        if not pairs:
            pairs = _discover_pairs_in_dir(self.input_dir)
        if not pairs:
            time_stamp("[contam] no FASTQs for Kraken2")
            return None

        log_dir = self.report_dir / "logs"
        db_path = _resolve_kraken_db_path(self.kraken_db)
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
        jobs: List[Dict[str, object]] = []

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
            jobs.append(
                dict(
                    sample=sample,
                    r1=r1_for_kraken,
                    r2=r2_for_kraken,
                    db=db_path,
                    kreport=krep,
                    kout=kout,
                    log_dir=log_dir,
                    gzip_input=gzip_input,
                    total_records=sample_total_records,
                )
            )

        if not jobs:
            time_stamp("[contam] no Kraken2 jobs to run")
            return None

        total_threads = max(1, int(self.threads))
        if total_threads < 16:
            workers = 1
        else:
            workers = max(1, total_threads // 16)
        workers = min(workers, len(jobs))
        threads_per_job = max(1, total_threads // max(1, workers))

        time_stamp(f"[contam] Kraken2 workers={workers}, threads/job={threads_per_job}")
        for job in jobs:
            job["threads"] = threads_per_job

        failed_samples = set()
        if workers > 1:
            with ThreadPoolExecutor(max_workers=workers) as ex:
                futs = [ex.submit(_run_kraken_job, job) for job in jobs]
                for fut in as_completed(futs):
                    sample, err = fut.result()
                    if err:
                        failed_samples.add(sample)
                        _append_log_line(log_dir, "kraken2", err)
        else:
            for job in jobs:
                sample, err = _run_kraken_job(job)
                if err:
                    failed_samples.add(sample)
                    _append_log_line(log_dir, "kraken2", err)
                    continue

        for job in jobs:
            sample = job["sample"]
            if sample in failed_samples:
                continue
            krep = job["kreport"]
            if not krep.exists():
                continue
            sample_total_records = _kreport_total_reads(krep) or int(job["total_records"])
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

        sample_summary, cohort_summary = _build_non_target_summaries(big, self.species_name)
        if not sample_summary.empty:
            sample_summary_path = out_dir / "kraken_non_target_top_hits.tsv"
            sample_summary.to_csv(sample_summary_path, sep="\t", index=False)
            time_stamp(f"[contam] wrote {sample_summary_path}")
        else:
            time_stamp("[contam] no non-target species summary rows after applying target-species hint")

        if not cohort_summary.empty:
            cohort_summary_path = out_dir / "kraken_non_target_median_summary.tsv"
            cohort_summary.to_csv(cohort_summary_path, sep="\t", index=False)
            time_stamp(f"[contam] wrote {cohort_summary_path}")
            hint_label = _normalize_species_hint(self.species_name) or "none"
            time_stamp(f"[contam] target-species hint: {hint_label}")
            for row in cohort_summary.head(5).itertuples(index=False):
                time_stamp(
                    "[contam] top non-target median: "
                    f"{row.name} median={row.median_freq_pct:.3f}% "
                    f"max={row.max_freq_pct:.3f}% samples={int(row.samples_detected)}"
                )

        return big
