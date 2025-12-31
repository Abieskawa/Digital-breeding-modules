#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from qmplot import manhattanplot, qqplot
from sklearn.metrics import roc_curve, auc

from Utils.utils import _resolve_outdir


class SNPPredictionModelTool:
    def __init__(self, config):
        """
        Lightweight helper for evaluation/plotting utilities.
        """
        self.config = config
        self.species_name = self.config.get("species_name")
        self.phenotype_column = self.config.get("phenotype_column", "Phenotype")
        self.chromosome_csv_path = self.config.get("chromosome_csv_path")

        self.output_dir = _resolve_outdir(
            self.config,
            key="prediction_output_dir",
            default="results",
            ensure_dir=True,
        )
        threads_raw = (self.config.get("num_threads") or self.config.get("threads") or "4").strip()
        try:
            self.num_threads = int(threads_raw)
        except ValueError:
            self.num_threads = 4

        self._create_directories()
        self._chromosome_loaded = False
        if self.chromosome_csv_path:
            self.load_chromosome_recode()
            self._chromosome_loaded = True

    def _create_directories(self):
        self.gwas_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="GWAS_dir",
            ensure_dir=True,
        )
        self.model_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="Model_result",
            ensure_dir=True,
        )
        self.auc_plot_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="AUC_dir",
            ensure_dir=True,
        )

    def load_chromosome_recode(self):
        if not self.chromosome_csv_path:
            raise ValueError("chromosome_csv_path is required for Manhattan/QQ plots")
        self.chromosome_mapping = pd.read_csv(self.chromosome_csv_path)
        self.chr_name_to_recoded = dict(
            zip(self.chromosome_mapping["OriginalName"], self.chromosome_mapping["RecodedNumber"])
        )
        self.recoded_to_label = dict(
            zip(self.chromosome_mapping["RecodedNumber"].astype(str), self.chromosome_mapping["Label"])
        )

    def plot_performance_comparison(self):
        """
        Generate comparison plots for n_snps vs. models and n_pcs vs. models in each fold,
        reading data from saved .npy files for probabilities and true labels.
        """
        probabilities_files = [f for f in os.listdir(self.model_output_dir) if f.startswith("Probabilities")]
        testy_files = [f for f in os.listdir(self.model_output_dir) if f.startswith("Testy")]

        probabilities_files.sort()
        testy_files.sort()

        auc_scores = []

        for prob_file, testy_file in zip(probabilities_files, testy_files):
            meta = prob_file.split("_")
            model_name = meta[1]
            fold = next((frag.replace("Fold", "") for frag in meta if frag.startswith("Fold")), None)
            n_pcs = int(next((frag.replace("PC", "") for frag in meta if "PC" in frag), None))
            n_snps = int(next((frag.replace("SNPs.npy", "") for frag in meta if "SNPs.npy" in frag), None))

            y_scores = np.load(os.path.join(self.model_output_dir, prob_file))
            test_y = np.load(os.path.join(self.model_output_dir, testy_file))

            fpr, tpr, _ = roc_curve(test_y, y_scores)

            auc_scores.append(
                {
                    "model_name": model_name,
                    "fold": fold,
                    "n_pcs": n_pcs,
                    "n_snps": n_snps,
                    "fpr": fpr,
                    "tpr": tpr,
                }
            )

        auc_df = pd.DataFrame(auc_scores)
        auc_df["roc_auc"] = auc_df.apply(
            lambda row: auc(row["fpr"], row["tpr"]) if isinstance(row["fpr"], (list, np.ndarray)) and isinstance(row["tpr"], (list, np.ndarray)) and len(row["fpr"]) == len(row["tpr"]) else np.nan,
            axis=1,
        )

        model_colors = sns.color_palette("tab10", n_colors=len(auc_df["model_name"].unique()))
        model_color_map = {model: model_colors[i] for i, model in enumerate(auc_df["model_name"].unique())}

        line_styles = [
            "solid",
            "dashed",
            "dotted",
            "dashdot",
            (0, (1, 1)),
            (0, (5, 5)),
            (0, (3, 5, 1, 5)),
        ]

        for fold in auc_df["fold"].unique():
            for n_snp in auc_df["n_snps"].unique():
                plt.figure(figsize=(6, 6))
                fold_nsnp_data = auc_df[(auc_df["fold"] == fold) & (auc_df["n_snps"] == n_snp)]
                sorted_data = fold_nsnp_data.sort_values(by=["model_name", "n_pcs"])
                plt.plot([0, 1], [0, 1], "k--", label="Random (y=x)")
                for _, row in sorted_data.iterrows():
                    model = row["model_name"]
                    fpr = row["fpr"]
                    tpr = row["tpr"]
                    npc_or_nsnps_styles = {n: line_styles[i % len(line_styles)] for i, n in enumerate(sorted(auc_df["n_pcs"].unique()))}

                    plt.plot(
                        fpr,
                        tpr,
                        label=f"{model}, {row['n_pcs']}PCs (AUC = {row['roc_auc']:.2f})",
                        color=model_color_map[model],
                        linestyle=npc_or_nsnps_styles[row["n_pcs"]],
                    )

                plt.title(
                    f"ROC curve (Fold {fold}, {n_snp}SNPs): PCs vs. Models of {self.phenotype_column} in {self.species_name}"
                )
                plt.xlabel("False Positive Rate")
                plt.ylabel("True Positive Rate")
                plt.grid()
                plt.legend(
                    loc="upper left",
                    bbox_to_anchor=(1, 1),
                    title="Legend",
                    fontsize="small",
                    markerscale=0,
                    handletextpad=0.8,
                )
                output_path = os.path.join(self.auc_plot_output_dir, f"auc_comparison_fold{fold}_{n_snp}SNPs.png")
                plt.savefig(output_path, bbox_inches="tight")
                plt.close()
                print(f"Plot saved: {output_path}")

            for n_pc in auc_df["n_pcs"].unique():
                plt.figure(figsize=(6, 6))
                fold_npc_data = auc_df[(auc_df["fold"] == fold) & (auc_df["n_pcs"] == n_pc)]

                plt.plot([0, 1], [0, 1], "k--", label="Random (y=x)")

                sorted_data = fold_npc_data.sort_values(by=["model_name", "n_snps"])
                for _, row in sorted_data.iterrows():
                    model = row["model_name"]
                    fpr = row["fpr"]
                    tpr = row["tpr"]

                    npc_or_nsnps_styles = {n: line_styles[i % len(line_styles)] for i, n in enumerate(sorted(auc_df["n_snps"].unique()))}
                    plt.plot(
                        fpr,
                        tpr,
                        label=f"{model}, {row['n_snps']}SNPs (AUC = {row['roc_auc']:.2f})",
                        color=model_color_map[model],
                        linestyle=npc_or_nsnps_styles[row["n_snps"]],
                    )

                plt.title(
                    f"ROC curve (Fold {fold}, {n_pc}PCs): SNPs vs. Models of {self.phenotype_column} in {self.species_name}"
                )
                plt.xlabel("False Positive Rate")
                plt.ylabel("True Positive Rate")
                plt.grid()
                plt.legend(
                    loc="upper left",
                    bbox_to_anchor=(1, 1),
                    title="Legend",
                    fontsize="small",
                    markerscale=0,
                    handletextpad=0.8,
                )
                output_path = os.path.join(self.auc_plot_output_dir, f"auc_comparison_fold{fold}_{n_pc}PCs.png")
                plt.savefig(output_path, bbox_inches="tight")
                plt.close()
                print(f"Plot saved: {output_path}")

        auc_df[["fpr", "tpr"]] = auc_df[["fpr", "tpr"]].applymap(lambda x: f"{x}")
        summary_csv = os.path.join(self.auc_plot_output_dir, "auc_summary.csv")
        auc_df.to_csv(summary_csv, index=False)
        print(f"Summary CSV saved: {summary_csv}")

    def plot_manhattan_and_qq(self, file_fold_npc_index):
        """
        Generate Manhattan and QQ plots for GWAS results.
        """
        if not self._chromosome_loaded:
            self.load_chromosome_recode()
            self._chromosome_loaded = True
        current_dir = os.path.join(self.gwas_output_dir, file_fold_npc_index)
        os.makedirs(current_dir, exist_ok=True)

        p_value_file = os.path.join(
            current_dir, f"train_{file_fold_npc_index}_{self.phenotype_column}_GWAS_result.txt"
        )
        df = pd.read_table(p_value_file, sep="\t")

        df = df.sort_values(by=["chr", "pos"])
        df = df.rename(columns={"chr": "#CHROM", "pos": "POS", "p_value": "P", "taxa": "ID"})

        recode_dict = dict(
            zip(self.chromosome_mapping["OriginalName"], self.chromosome_mapping["Label"])
        )
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
