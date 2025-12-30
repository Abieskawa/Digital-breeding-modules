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
    _discover_pairs_in_dir,
    _kreport_total_reads,
    _prepare_kraken_inputs,
    _run_kraken2,
    _summarize_kreport_species,
    call_log,
    time_stamp,
)
from Evaluation.mapping_base import MappingEvalBase


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
