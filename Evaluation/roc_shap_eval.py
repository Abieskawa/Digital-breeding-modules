#!/usr/bin/env python3
import os
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap
import torch
from sklearn.metrics import roc_curve, auc
from xgboost import XGBClassifier

from Utils.utils import _resolve_outdir
from Prediction_model.model_construction import nnModel


def _prediction_plot_dirs(config: Dict[str, str]) -> Dict[str, object]:
    output_dir = _resolve_outdir(
        config,
        key="prediction_output_dir",
        default="results",
        ensure_dir=True,
    )
    model_output_dir = _resolve_outdir(
        base_outdir=output_dir,
        subdir="Model_result",
        ensure_dir=True,
    )
    auc_plot_output_dir = _resolve_outdir(
        base_outdir=output_dir,
        subdir="AUC_dir",
        ensure_dir=True,
    )
    shap_output_dir = _resolve_outdir(
        base_outdir=output_dir,
        subdir="Shap_dir",
        ensure_dir=True,
    )
    saved_models_dir = _resolve_outdir(
        base_outdir=model_output_dir,
        subdir="Saved_models",
        ensure_dir=True,
    )
    return {
        "output_dir": output_dir,
        "model_output_dir": model_output_dir,
        "auc_plot_output_dir": auc_plot_output_dir,
        "shap_output_dir": shap_output_dir,
        "saved_models_dir": saved_models_dir,
    }


def _plot_performance_comparison(
    *,
    model_output_dir: str,
    auc_plot_output_dir: str,
    phenotype_column: str,
    species_name: Optional[str],
) -> None:
    probabilities_files = [f for f in os.listdir(model_output_dir) if f.startswith("Probabilities")]
    testy_files = [f for f in os.listdir(model_output_dir) if f.startswith("Testy")]

    probabilities_files.sort()
    testy_files.sort()

    auc_scores = []

    for prob_file, testy_file in zip(probabilities_files, testy_files):
        meta = prob_file.split("_")
        model_name = meta[1]
        fold = next((frag.replace("Fold", "") for frag in meta if frag.startswith("Fold")), None)
        n_pcs = int(next((frag.replace("PC", "") for frag in meta if "PC" in frag), None))
        n_snps = int(next((frag.replace("SNPs.npy", "") for frag in meta if "SNPs.npy" in frag), None))

        y_scores = np.load(os.path.join(model_output_dir, prob_file))
        test_y = np.load(os.path.join(model_output_dir, testy_file))

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
        lambda row: auc(row["fpr"], row["tpr"])
        if isinstance(row["fpr"], (list, np.ndarray))
        and isinstance(row["tpr"], (list, np.ndarray))
        and len(row["fpr"]) == len(row["tpr"])
        else np.nan,
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
            npc_or_nsnps_styles = {
                n: line_styles[i % len(line_styles)] for i, n in enumerate(sorted(auc_df["n_pcs"].unique()))
            }
            for _, row in sorted_data.iterrows():
                model = row["model_name"]
                fpr = row["fpr"]
                tpr = row["tpr"]

                plt.plot(
                    fpr,
                    tpr,
                    label=f"{model}, {row['n_pcs']}PCs (AUC = {row['roc_auc']:.2f})",
                    color=model_color_map[model],
                    linestyle=npc_or_nsnps_styles[row["n_pcs"]],
                )

            plt.title(
                f"ROC curve (Fold {fold}, {n_snp}SNPs): PCs vs. Models of {phenotype_column} in {species_name}"
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
            output_path = os.path.join(auc_plot_output_dir, f"auc_comparison_fold{fold}_{n_snp}SNPs.png")
            plt.savefig(output_path, bbox_inches="tight")
            plt.close()
            print(f"Plot saved: {output_path}")

        for n_pc in auc_df["n_pcs"].unique():
            plt.figure(figsize=(6, 6))
            fold_npc_data = auc_df[(auc_df["fold"] == fold) & (auc_df["n_pcs"] == n_pc)]

            plt.plot([0, 1], [0, 1], "k--", label="Random (y=x)")

            sorted_data = fold_npc_data.sort_values(by=["model_name", "n_snps"])
            npc_or_nsnps_styles = {
                n: line_styles[i % len(line_styles)] for i, n in enumerate(sorted(auc_df["n_snps"].unique()))
            }
            for _, row in sorted_data.iterrows():
                model = row["model_name"]
                fpr = row["fpr"]
                tpr = row["tpr"]

                plt.plot(
                    fpr,
                    tpr,
                    label=f"{model}, {row['n_snps']}SNPs (AUC = {row['roc_auc']:.2f})",
                    color=model_color_map[model],
                    linestyle=npc_or_nsnps_styles[row["n_snps"]],
                )

            plt.title(
                f"ROC curve (Fold {fold}, {n_pc}PCs): SNPs vs. Models of {phenotype_column} in {species_name}"
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
            output_path = os.path.join(auc_plot_output_dir, f"auc_comparison_fold{fold}_{n_pc}PCs.png")
            plt.savefig(output_path, bbox_inches="tight")
            plt.close()
            print(f"Plot saved: {output_path}")

    auc_df[["fpr", "tpr"]] = auc_df[["fpr", "tpr"]].applymap(lambda x: f"{x}")
    summary_csv = os.path.join(auc_plot_output_dir, "auc_summary.csv")
    auc_df.to_csv(summary_csv, index=False)
    print(f"Summary CSV saved: {summary_csv}")


def _parse_model_artifact(name: str) -> Optional[Dict[str, str]]:
    if "_Fold_" not in name:
        return None
    base, ext = os.path.splitext(name)
    if not ext:
        return None
    if ext not in {".pt", ".json", ".joblib", ".pkl"}:
        return None
    if "_PC" not in base or "_Fold_" not in base or "SNPs" not in base:
        return None
    model_name, rest = base.split("_Fold_", 1)
    if "_PC" not in rest:
        return None
    file_index, rest = rest.split("_PC", 1)
    if "_" not in rest:
        return None
    pcs_raw, snps_raw = rest.split("_", 1)
    if not snps_raw.endswith("SNPs"):
        return None
    n_snps_raw = snps_raw.replace("SNPs", "")
    return {
        "model_name": model_name,
        "file_index": file_index,
        "n_pcs": pcs_raw,
        "n_snps": n_snps_raw,
        "ext": ext.lstrip("."),
        "base": base,
    }


def _load_feature_list(saved_models_dir: str, base: str) -> List[str]:
    feature_path = os.path.join(saved_models_dir, f"{base}.features.txt")
    with open(feature_path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def _load_train_test_arrays(
    model_output_dir: str,
    file_index: str,
    n_pcs: str,
    n_snps: str,
    feature_list: List[str],
    phenotype_column: str,
) -> Dict[str, np.ndarray]:
    fold_tag = f"{file_index}_{n_pcs}PCs_{n_snps}SNPs"
    current_dir = os.path.join(model_output_dir, fold_tag)
    train_csv_path = os.path.join(current_dir, f"train_selected_snps_{fold_tag}.csv")
    test_csv_path = os.path.join(current_dir, f"test_selected_snps_{fold_tag}.csv")

    train_data = pd.read_csv(train_csv_path, index_col="taxa")
    test_data = pd.read_csv(test_csv_path, index_col="taxa")

    missing = [col for col in feature_list if col not in train_data.columns]
    if missing:
        raise ValueError(f"Missing features in training data: {missing}")

    train_y = train_data[phenotype_column].values.astype(int)
    test_y = test_data[phenotype_column].values.astype(int)
    train_X = train_data[feature_list].values
    test_X = test_data[feature_list].values
    return {
        "train_X": train_X,
        "test_X": test_X,
        "train_y": train_y,
        "test_y": test_y,
    }


def _load_model(
    *,
    ext: str,
    path: str,
    feature_list: List[str],
    data_type: str,
) -> object:
    if ext == "pt":
        num_o = 1 if data_type == "binary" else 2
        model = nnModel(num_i=len(feature_list), num_h=64, num_o=num_o)
        payload = torch.load(path, map_location="cpu")
        state = payload["state_dict"] if isinstance(payload, dict) and "state_dict" in payload else payload
        model.load_state_dict(state)
        model.eval()
        return model
    if ext == "json":
        model = XGBClassifier()
        model.load_model(path)
        return model
    if ext == "joblib":
        import joblib
        payload = joblib.load(path)
        return payload["model"] if isinstance(payload, dict) and "model" in payload else payload
    if ext == "pkl":
        import pickle
        with open(path, "rb") as f:
            payload = pickle.load(f)
        return payload["model"] if isinstance(payload, dict) and "model" in payload else payload
    raise ValueError(f"Unsupported model artifact: {path}")


def _generate_shap_outputs(
    *,
    model_output_dir: str,
    saved_models_dir: str,
    shap_output_dir: str,
    phenotype_column: str,
    data_type: str,
) -> None:
    if not os.path.isdir(saved_models_dir):
        print(f"No saved models found under {saved_models_dir}")
        return
    artifacts = sorted(os.listdir(saved_models_dir))
    for artifact in artifacts:
        parsed = _parse_model_artifact(artifact)
        if not parsed:
            continue
        model_name = parsed["model_name"]
        file_index = parsed["file_index"]
        n_pcs = parsed["n_pcs"]
        n_snps = parsed["n_snps"]
        base = parsed["base"]
        ext = parsed["ext"]
        model_path = os.path.join(saved_models_dir, artifact)

        feature_list = _load_feature_list(saved_models_dir, base)
        arrays = _load_train_test_arrays(
            model_output_dir,
            file_index,
            n_pcs,
            n_snps,
            feature_list,
            phenotype_column,
        )
        train_X = arrays["train_X"]
        test_X = arrays["test_X"]

        model = _load_model(
            ext=ext,
            path=model_path,
            feature_list=feature_list,
            data_type=data_type,
        )

        if model_name in ["RandomForest", "XGBoost"]:
            explainer = shap.TreeExplainer(model)
            shap_values = explainer.shap_values(test_X)
            if model_name == "RandomForest":
                shap_values_to_plot = shap_values[:, :, 1]
            else:
                shap_values_to_plot = shap_values
        elif model_name == "MLP":
            X_train_cpu = torch.tensor(train_X, dtype=torch.float32)
            X_test_cpu = torch.tensor(test_X, dtype=torch.float32)
            explainer = shap.GradientExplainer(model, X_train_cpu)
            shap_values = explainer.shap_values(X_test_cpu)
            shap_values_to_plot = np.squeeze(shap_values)
        else:
            if not hasattr(model, "predict_proba"):
                print(f"Skipping SHAP for {model_name}: predict_proba unavailable")
                continue
            explainer = shap.KernelExplainer(model.predict_proba, train_X)
            shap_values = explainer.shap_values(test_X)
            shap_values_to_plot = shap_values[:, :, 1]

        df_shap_values = pd.DataFrame(shap_values_to_plot, columns=feature_list)
        df_feature_importance = pd.DataFrame(columns=["feature", "importance"])
        for col in df_shap_values.columns:
            importance = df_shap_values[col].abs().mean()
            df_feature_importance.loc[len(df_feature_importance)] = [col, importance]
            df_feature_importance = df_feature_importance.sort_values("importance", ascending=False)

        df_feature_importance_file = os.path.join(
            model_output_dir,
            f"ShapTop10_{model_name}_Fold_{file_index}_PC{n_pcs}_{n_snps}SNPs.csv",
        )
        df_feature_importance.to_csv(df_feature_importance_file, index=False)

        plt.figure()
        shap.summary_plot(shap_values_to_plot, test_X, feature_names=feature_list, show=False)
        plt.title(f"{model_name} SHAP Summary Plot (Fold {file_index},PC{n_pcs},{n_snps}SNPs)")
        shap_plot = os.path.join(
            shap_output_dir,
            f"{model_name}_SHAP_Fold_{file_index}_PC{n_pcs}_{n_snps}SNPs.png",
        )
        plt.savefig(shap_plot)
        plt.close()


class RocShapEvaluator:
    def __init__(self, config: Dict[str, str]):
        self.config = config
        self.species_name = self.config.get("species_name")
        self.phenotype_column = self.config.get("phenotype_column", "Phenotype")
        self.data_type = self.config.get("data_type", "binary")
        dirs = _prediction_plot_dirs(self.config)
        self.model_output_dir = dirs["model_output_dir"]
        self.auc_plot_output_dir = dirs["auc_plot_output_dir"]
        self.shap_output_dir = dirs["shap_output_dir"]
        self.saved_models_dir = dirs["saved_models_dir"]

    def run(self) -> None:
        probabilities_files = [f for f in os.listdir(self.model_output_dir) if f.startswith("Probabilities")]
        testy_files = [f for f in os.listdir(self.model_output_dir) if f.startswith("Testy")]
        if not probabilities_files or not testy_files:
            print(f"No prediction probability outputs found under {self.model_output_dir}")
            return
        _plot_performance_comparison(
            model_output_dir=str(self.model_output_dir),
            auc_plot_output_dir=str(self.auc_plot_output_dir),
            phenotype_column=self.phenotype_column,
            species_name=self.species_name,
        )
        _generate_shap_outputs(
            model_output_dir=str(self.model_output_dir),
            saved_models_dir=str(self.saved_models_dir),
            shap_output_dir=str(self.shap_output_dir),
            phenotype_column=self.phenotype_column,
            data_type=self.data_type,
        )
