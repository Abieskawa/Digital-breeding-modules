#!/usr/bin/env python3
import os
from typing import Dict, List, Optional

import pandas as pd
import matplotlib.pyplot as plt
from qmplot import manhattanplot, qqplot

from Utils.utils import _resolve_outdir, parallel_process


def _num_threads(config: Dict[str, str]) -> int:
    threads_raw = (config.get("num_threads") or config.get("threads") or "4").strip()
    try:
        return int(threads_raw)
    except ValueError:
        return 4


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


def _collect_plot_targets(gwas_dir: str, phenotype_column: str) -> List[str]:
    targets: List[str] = []
    if not os.path.isdir(gwas_dir):
        return targets
    for name in sorted(os.listdir(gwas_dir)):
        current_dir = os.path.join(gwas_dir, name)
        if not os.path.isdir(current_dir):
            continue
        p_value_file = os.path.join(
            current_dir, f"train_{name}_{phenotype_column}_GWAS_result.txt"
        )
        if os.path.isfile(p_value_file):
            targets.append(name)
    return targets


def _load_chromosome_recode(chromosome_csv_path: str) -> Dict[str, str]:
    if not chromosome_csv_path:
        raise ValueError("chromosome_csv_path is required for Manhattan/QQ plots")
    chromosome_mapping = pd.read_csv(chromosome_csv_path)
    recode_dict = dict(
        zip(chromosome_mapping["OriginalName"], chromosome_mapping["Label"])
    )
    return recode_dict


def _plot_manhattan_and_qq(job: Dict[str, object]) -> str:
    config = job["config"]
    file_fold_npc_index = str(job["file_fold_npc_index"])

    species_name = config.get("species_name")
    phenotype_column = config.get("phenotype_column", "Phenotype")
    chromosome_csv_path = config.get("chromosome_csv_path")
    recode_dict = _load_chromosome_recode(chromosome_csv_path)

    gwas_output_dir = _gwas_plot_dirs(config)["gwas_output_dir"]
    current_dir = os.path.join(gwas_output_dir, file_fold_npc_index)
    os.makedirs(current_dir, exist_ok=True)

    p_value_file = os.path.join(
        current_dir, f"train_{file_fold_npc_index}_{phenotype_column}_GWAS_result.txt"
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
        title=f"Manhattan Plot - Fold {file_fold_npc_index} of {phenotype_column} in {species_name}",
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
        f"QQ Plot - Fold {file_fold_npc_index} of {phenotype_column} in {species_name}"
    )
    output_qqplot_path = os.path.join(current_dir, f"qq_plot_shuffle-fold_{file_fold_npc_index}.png")
    plt.savefig(output_qqplot_path)
    plt.close()
    return file_fold_npc_index


class GWASPCAEvaluator:
    def __init__(self, config: Dict[str, str]):
        self.config = config
        self.species_name = self.config.get("species_name")
        self.phenotype_column = self.config.get("phenotype_column", "Phenotype")
        self.chromosome_csv_path = self.config.get("chromosome_csv_path")
        self.num_threads = _num_threads(self.config)
        self.gwas_output_dir = _gwas_plot_dirs(self.config)["gwas_output_dir"]

    def run(self) -> None:
        targets = _collect_plot_targets(self.gwas_output_dir, self.phenotype_column)
        if not targets:
            print(f"No GWAS results found under {self.gwas_output_dir}")
            return
        jobs = [{"config": self.config, "file_fold_npc_index": target} for target in targets]
        parallel_process(
            jobs,
            _plot_manhattan_and_qq,
            max_workers=self.num_threads,
        )
