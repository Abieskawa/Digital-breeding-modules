#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from Utils.utils import time_stamp
from Evaluation.mapping_base import MappingEvalBase
from Evaluation.mapping_utils import (
    _annotate_box_fliers,
    _count_mapped_primary_reads,
    _count_reads_in_fastq,
    _count_unmapped_primary_reads,
    _discover_pairs_in_dir,
)


class MappingYieldEvaluator(MappingEvalBase):
    """
    Step 2: Mapping yield ratio = unique (primary, non-supplementary) BAM reads / fastp reads.
    """

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
