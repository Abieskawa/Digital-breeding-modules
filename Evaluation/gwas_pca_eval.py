#!/usr/bin/env python3
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Utils.utils import _load_chromosome_recode, _resolve_outdir


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
        from qmplot import manhattanplot, qqplot

        recode_dict = _load_chromosome_recode(self.chromosome_csv_path, required=True)

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


class PCAEvaluator:
    def __init__(self, config: Dict[str, str]):
        self.config = config
        self.species_name = self.config.get("species_name")
        self.phenotype_csv_path = self.config.get("phenotype_csv_path")
        self.phenotype_column = self.config.get("phenotype_column", "Phenotype")
        self.data_type = self.config.get("data_type", "binary")

        self.output_dir = _resolve_outdir(
            self.config,
            key="prediction_output_dir",
            default="results",
            ensure_dir=True,
        )
        self.pca_scree_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="evaluation/PCA_Scree_Plot",
            ensure_dir=True,
        )
        self.pca_eval_output_dir = self.pca_scree_output_dir

    @staticmethod
    def _format_fold_tag(file_fold_index: str) -> str:
        parts = file_fold_index.split("-", 1)
        if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
            return f"fold{parts[1]}"
        if file_fold_index.isdigit():
            return f"fold{file_fold_index}"
        safe = "".join(c if c.isalnum() else "_" for c in file_fold_index)
        return f"fold{safe}"

    def _load_phenotypes(self) -> pd.DataFrame:
        phenotype_df = pd.read_csv(self.phenotype_csv_path, sep=None, engine="python")
        if "taxa" not in phenotype_df.columns:
            raise ValueError("phenotype_csv_path must contain a 'taxa' column.")
        if self.phenotype_column not in phenotype_df.columns:
            raise ValueError(f"phenotype column not found: {self.phenotype_column}")
        phenotype_df = phenotype_df[["taxa", self.phenotype_column]].dropna()
        phenotype_df["taxa"] = phenotype_df["taxa"].astype(str)
        return phenotype_df.set_index("taxa")

    def _prepare_pca_plot_df(
        self,
        principal_df: pd.DataFrame,
        phenotype_df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        if phenotype_df is None:
            phenotype_df = self._load_phenotypes()
        plot_df = principal_df.copy()
        plot_df["taxa"] = plot_df["taxa"].astype(str)
        plot_df = plot_df.set_index("taxa")
        plot_df = plot_df.join(phenotype_df, how="inner")
        return plot_df

    @staticmethod
    def _axis_label(pc_idx: int, explained_variance_ratio: np.ndarray) -> str:
        if 0 <= pc_idx - 1 < len(explained_variance_ratio):
            pct = explained_variance_ratio[pc_idx - 1] * 100
            return f"Principal component {pc_idx} ({pct:.1f}%)"
        return f"Principal component {pc_idx}"

    def _plot_pca_scatter(
        self,
        plot_df: pd.DataFrame,
        explained_variance_ratio: np.ndarray,
        out_path: str,
        title: str,
    ) -> None:
        if plot_df.empty:
            return
        if not {"PC1", "PC2", "PC3"}.issubset(plot_df.columns):
            return

        fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
        pairs = [(1, 2), (2, 3), (1, 3)]
        phenos = plot_df[self.phenotype_column]
        is_numeric = pd.api.types.is_numeric_dtype(phenos)
        unique_count = phenos.nunique(dropna=True)
        continuous = self.data_type == "regression" or (is_numeric and unique_count > 10)

        if continuous:
            scatter = None
            for ax, (x_idx, y_idx) in zip(axes, pairs):
                scatter = ax.scatter(
                    plot_df[f"PC{x_idx}"],
                    plot_df[f"PC{y_idx}"],
                    c=phenos.values,
                    cmap="viridis",
                    s=12,
                    alpha=0.75,
                )
                ax.set_xlabel(self._axis_label(x_idx, explained_variance_ratio))
                ax.set_ylabel(self._axis_label(y_idx, explained_variance_ratio))
            if scatter is not None:
                fig.colorbar(scatter, ax=axes, label=self.phenotype_column, shrink=0.85)
        else:
            pheno_labels = phenos.astype(str)
            groups = sorted(pheno_labels.unique())
            base_cmap = plt.get_cmap("tab20", max(len(groups) + 2, 1))
            color_map = {}
            for idx, group in enumerate(groups):
                if idx == 0:
                    color_map[group] = "skyblue"
                elif idx == 1:
                    color_map[group] = "salmon"
                else:
                    color_map[group] = base_cmap(idx + 2)
            for ax, (x_idx, y_idx) in zip(axes, pairs):
                for group in groups:
                    subset = plot_df[pheno_labels == group]
                    ax.scatter(
                        subset[f"PC{x_idx}"],
                        subset[f"PC{y_idx}"],
                        label=group,
                        color=color_map[group],
                        s=12,
                        alpha=0.75,
                    )
                ax.set_xlabel(self._axis_label(x_idx, explained_variance_ratio))
                ax.set_ylabel(self._axis_label(y_idx, explained_variance_ratio))
            axes[-1].legend(title=self.phenotype_column, fontsize="small", loc="best")

        fig.suptitle(title)
        fig.savefig(out_path, dpi=300)
        plt.close(fig)

    def _plot_scree(
        self,
        explained_variance_ratio: np.ndarray,
        out_path: str,
        title: str,
    ) -> None:
        n_components = min(15, len(explained_variance_ratio))
        if n_components < 1:
            return
        explained = explained_variance_ratio[:n_components] * 100.0
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(range(1, n_components + 1), explained, marker="o")
        ax.set_xticks(range(1, n_components + 1))
        ax.set_ylim(0, 100)
        ax.set_xlabel("Principal component")
        ax.set_ylabel("Explained variance (%)")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        fig.savefig(out_path, dpi=300)
        plt.close(fig)

    def _plot_fold_pca(
        self,
        file_fold_index: str,
        principal_df: pd.DataFrame,
        explained_variance_ratio: np.ndarray,
    ) -> None:
        tag = self._format_fold_tag(file_fold_index)
        plot_df = self._prepare_pca_plot_df(principal_df)
        scatter_path = os.path.join(self.pca_eval_output_dir, f"pca_scatter_{tag}.png")
        scree_path = os.path.join(self.pca_eval_output_dir, f"scree_{tag}.png")

        title_base = f"PCA (CV Fold {file_fold_index})"
        if self.species_name:
            title_base = f"{title_base} - {self.species_name}"

        self._plot_pca_scatter(plot_df, explained_variance_ratio, scatter_path, title_base)
        self._plot_scree(
            explained_variance_ratio,
            scree_path,
            f"Scree Plot (CV Fold {file_fold_index})",
        )

    def collect_pca_targets(self) -> List[str]:
        if not os.path.isdir(self.pca_eval_output_dir) and not os.path.isdir(self.pca_scree_output_dir):
            return []
        targets: List[str] = []
        prefix = "pca_components_"
        suffix = ".csv"
        search_dirs = []
        if os.path.isdir(self.pca_eval_output_dir):
            search_dirs.append(self.pca_eval_output_dir)
        if os.path.isdir(self.pca_scree_output_dir):
            search_dirs.append(self.pca_scree_output_dir)
        for directory in search_dirs:
            for name in sorted(os.listdir(directory)):
                if not (name.startswith(prefix) and name.endswith(suffix)):
                    continue
                fold_id = name[len(prefix):-len(suffix)]
                if fold_id:
                    targets.append(fold_id)
        return targets

    def plot_fold_pca_from_files(self, file_fold_index: str) -> bool:
        candidates = [
            (
                os.path.join(self.pca_eval_output_dir, f"pca_components_{file_fold_index}.csv"),
                os.path.join(self.pca_eval_output_dir, f"explained_variance_{file_fold_index}.csv"),
            ),
            (
                os.path.join(self.pca_scree_output_dir, f"pca_components_{file_fold_index}.csv"),
                os.path.join(self.pca_scree_output_dir, f"explained_variance_{file_fold_index}.csv"),
            ),
        ]
        components_path = None
        variance_path = None
        for comp_path, var_path in candidates:
            if os.path.isfile(comp_path) and os.path.isfile(var_path):
                components_path = comp_path
                variance_path = var_path
                break
        if components_path is None or variance_path is None:
            return False

        principal_df = pd.read_csv(components_path)
        variance_df = pd.read_csv(variance_path)
        if "explained_variance_ratio" in variance_df.columns:
            explained = variance_df["explained_variance_ratio"].values
        else:
            explained = variance_df.iloc[:, 0].values

        self._plot_fold_pca(file_fold_index, principal_df, explained)
        return True

    def plot_global_qc(self, principal_df: pd.DataFrame, explained_variance_ratio: np.ndarray) -> None:
        plot_df = self._prepare_pca_plot_df(principal_df)
        scatter_path = os.path.join(self.pca_eval_output_dir, "pca_scatter_global_qc.png")
        scree_path = os.path.join(self.pca_eval_output_dir, "scree_global_qc.png")
        title_base = "Global QC PCA"
        if self.species_name:
            title_base = f"{title_base} - {self.species_name}"

        self._plot_pca_scatter(plot_df, explained_variance_ratio, scatter_path, title_base)
        self._plot_scree(explained_variance_ratio, scree_path, "Global QC Scree Plot")
