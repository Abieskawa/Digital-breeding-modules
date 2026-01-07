#!/usr/bin/env python3
import os
from typing import Dict, List, Optional

import pandas as pd
import matplotlib.pyplot as plt
from qmplot import manhattanplot, qqplot

from Utils.utils import _resolve_outdir


def _gwas_plot_dirs(config: Dict[str, str]) -> Dict[str, object]:
    output_dir = _resolve_outdir(
        config,
        key="prediction_output_dir",
        default="results",
        ensure_dir=True,
    )
    gwas_output_dir = _resolve_outdir(
        base_outdir=output_dir,
        subdir="GWAS_dir",
        ensure_dir=True,
    )
    return {
        "output_dir": output_dir,
        "gwas_output_dir": gwas_output_dir,
    }


def _load_chromosome_recode(chromosome_csv_path: str) -> Dict[str, str]:
    if not chromosome_csv_path:
        raise ValueError("chromosome_csv_path is required for Manhattan/QQ plots")
    chromosome_mapping = pd.read_csv(chromosome_csv_path)
    recode_dict = dict(
        zip(chromosome_mapping["OriginalName"], chromosome_mapping["Label"])
    )
    return recode_dict


class GWASPCAEvaluator:
    def __init__(self, config: Dict[str, str]):
        self.config = config
        self.species_name = self.config.get("species_name")
        self.phenotype_column = self.config.get("phenotype_column", "Phenotype")
        self.chromosome_csv_path = self.config.get("chromosome_csv_path")
        self.gwas_output_dir = _gwas_plot_dirs(self.config)["gwas_output_dir"]

    def collect_targets(self) -> List[str]:  
        targets: List[str] = []
        if not os.path.isdir(self.gwas_output_dir):
            return targets
        for name in sorted(os.listdir(self.gwas_output_dir)):
            current_dir = os.path.join(self.gwas_output_dir, name)
            if not os.path.isdir(current_dir):
                continue
            p_value_file = os.path.join(
                current_dir, f"train_{name}_{self.phenotype_column}_GWAS_result.txt"
            )
            if os.path.isfile(p_value_file):
                targets.append(name)
        return targets

    def plot_target(self, file_fold_npc_index: str) -> None:
        recode_dict = _load_chromosome_recode(self.chromosome_csv_path)

        current_dir = os.path.join(self.gwas_output_dir, file_fold_npc_index)
        os.makedirs(current_dir, exist_ok=True)

        p_value_file = os.path.join(
            current_dir, f"train_{file_fold_npc_index}_{self.phenotype_column}_GWAS_result.txt"
        )
        df = pd.read_table(p_value_file, sep="\t")

        df = df.sort_values(by=["chr", "pos"])
        df = df.rename(columns={"chr": "#CHROM", "pos": "POS", "p_value": "P", "taxa": "ID"})
        df["#CHROM"] = df["#CHROM"].map(recode_dict)
        df = df[df["P"] != 0]

        fig, ax = plt.subplots(figsize=(12, 6))
        manhattanplot(
            data=df,
            sign_line_cols=["#D62728", "#2CA02C"],
            hline_kws={"linestyle": "--", "lw": 1.3},
            xticklabel_kws={"rotation": 45},
            sign_marker_p=0.05 / df.shape[0],
            sign_marker_color="r",
            title=f"Manhattan Plot - Fold {file_fold_npc_index} of {self.phenotype_column} in {self.species_name}",
            xlabel="Chromosome",
            ylabel=r"$-log_{10}{(P)}$",
            ax=ax,
        )
        output_manplot_path = os.path.join(current_dir, f"manhattan_plot_{file_fold_npc_index}.png")
        plt.savefig(output_manplot_path)
        plt.close()

        plt.figure(figsize=(6, 6))
        qqplot(data=df["P"])
        plt.title(
            f"QQ Plot - Fold {file_fold_npc_index} of {self.phenotype_column} in {self.species_name}"
        )
        output_qqplot_path = os.path.join(current_dir, f"qq_plot_shuffle-fold_{file_fold_npc_index}.png")
        plt.savefig(output_qqplot_path)
        plt.close()
