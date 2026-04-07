#!/usr/bin/env python3
"""
Step 5: Polygenic Risk Score (PRS) calculation (post-GWAS).

This script derives SNP weights from GWAS summary statistics and applies them
to genotype matrices to produce PRS scores per sample.
"""
import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from Prediction_model.snp_numeric_transformer import SNP_numerical
from Utils.cv_preparation import CV_preparation
from Utils.utils import (
    _resolve_outdir,
    coerce_bool,
    encode_phenotype_series,
    parse_int_list,
    read_config,
    time_stamp,
)


def _dedupe(seq: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for item in seq:
        if item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out


def _parse_prs_top_specs(raw: str) -> List[Dict[str, object]]:
    """
    Parse PRS selection specs from config.
    Accepts comma/semicolon-separated range tokens only:
      - "1-10,11-20,21-50"
    Returns list of dicts with start/end ranks and output labels.
    """
    if raw is None:
        return []
    cleaned = str(raw).strip()
    if not cleaned:
        return []
    tokens = [tok.strip() for tok in cleaned.replace(";", ",").split(",") if tok.strip()]
    specs: List[Dict[str, object]] = []
    seen = set()
    for tok in tokens:
        if "-" not in tok:
            raise ValueError(
                f"Invalid PRS token '{tok}'. Expected range format like '1-10,11-20,21-50'."
            )
        left, right = [part.strip() for part in tok.split("-", 1)]
        if not (left.isdigit() and right.isdigit()):
            raise ValueError(f"Invalid PRS range token: {tok}")
        start = int(left)
        end = int(right)
        if start <= 0 or end <= 0:
            raise ValueError(f"PRS ranges must be positive integers: {tok}")
        if start > end:
            start, end = end, start
        label = f"prs_bin{start}_{end}"
        if label in seen:
            continue
        specs.append({"start": start, "end": end, "label": label, "raw": tok})
        seen.add(label)
    return specs


class PRSCalculator:
    def __init__(self, config: Dict[str, str]):
        self.config = config
        self.phenotype_csv_path = self.config.get("phenotype_csv_path")
        self.vcf_file_path = self.config.get("vcf_file_path")
        self.phenotype_column = self.config.get("phenotype_column", "Phenotype")

        self.output_dir = _resolve_outdir(
            self.config,
            key="prediction_output_dir",
            default="results",
            ensure_dir=True,
        )
        self.gwas_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="GWAS_dir",
            ensure_dir=True,
        )
        self.prs_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="PRS_dir",
            ensure_dir=True,
        )

        raw_prs_top = (self.config.get("prs_top_n_snps_list", "") or "").strip()
        if not raw_prs_top:
            raw_prs_top = (self.config.get("top_n_snps_list", "") or "").strip()
        if not raw_prs_top:
            raw_prs_top = "1-10"
        self.prs_top_specs = _parse_prs_top_specs(raw_prs_top)

        self.weight_mode = (self.config.get("prs_weight_mode", "logp") or "logp").strip().lower()
        self.min_pvalue = float(self.config.get("prs_min_pvalue", "0") or 0)
        self.include_phenotype = coerce_bool(self.config.get("prs_include_phenotype", "true"), True)

        self.snp_tool = SNP_numerical(self.config)

    def _genotypes_to_numeric_no_impute(self, genotype_df: pd.DataFrame) -> pd.DataFrame:
        genotype_numeric_df = genotype_df.apply(
            lambda col: col.map(lambda x: self.snp_tool._map_gt(self.snp_tool._to_gt(x)))
        )
        return genotype_numeric_df.astype(float)

    def _prepare_train_test_numeric(
        self,
        train_geno_raw: pd.DataFrame,
        test_geno_raw: pd.DataFrame,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        train_num = self._genotypes_to_numeric_no_impute(train_geno_raw)
        test_num = self._genotypes_to_numeric_no_impute(test_geno_raw) if not test_geno_raw.empty else pd.DataFrame(index=test_geno_raw.index)

        all_missing = train_num.isna().all(axis=0)
        if all_missing.any():
            train_num = train_num.loc[:, ~all_missing]
            test_num = test_num.reindex(columns=train_num.columns)

        if train_num.shape[1] == 0:
            raise ValueError("All selected PRS variants are missing in the train fold.")

        medians = train_num.median(axis=0, skipna=True)
        valid_cols = medians[medians.notna()].index.tolist()
        if not valid_cols:
            raise ValueError("Unable to derive PRS train-fold medians for imputation.")

        train_num = train_num[valid_cols].fillna(medians.loc[valid_cols])
        test_num = test_num.reindex(columns=valid_cols).fillna(medians.loc[valid_cols])
        return train_num, test_num

    def _load_phenotypes(self) -> pd.DataFrame:
        if not self.phenotype_csv_path:
            raise ValueError("phenotype_csv_path is required for PRS.")
        pheno_df = pd.read_csv(self.phenotype_csv_path, sep=None, engine="python")
        if "taxa" not in pheno_df.columns:
            raise ValueError("phenotype_csv_path must contain a 'taxa' column.")
        if self.phenotype_column not in pheno_df.columns:
            raise ValueError(f"phenotype column not found: {self.phenotype_column}")
        pheno_df = pheno_df[["taxa", self.phenotype_column]].dropna()
        pheno_df[self.phenotype_column], _ = encode_phenotype_series(
            pheno_df[self.phenotype_column],
            log_prefix="PRS phenotype",
        )
        pheno_df = pheno_df.set_index("taxa")
        return pheno_df

    def _load_gwas_results(self, path: Path) -> pd.DataFrame:
        df = pd.read_csv(path, sep="\t")
        snp_col = None
        for cand in ("taxa", "SNP", "Marker", "ID"):
            if cand in df.columns:
                snp_col = cand
                break
        if snp_col is None:
            raise ValueError(f"GWAS file missing SNP column: {path}")
        pval_col = None
        for cand in ("p_value", "pvalue", "P"):
            if cand in df.columns:
                pval_col = cand
                break
        if pval_col is None:
            raise ValueError(f"GWAS file missing p-value column: {path}")

        df = df[[snp_col, pval_col]].dropna()
        df[pval_col] = pd.to_numeric(df[pval_col], errors="coerce")
        df = df[(df[pval_col] > 0) & df[pval_col].notna()]
        if self.min_pvalue > 0:
            df = df[df[pval_col] >= self.min_pvalue]
        if df.empty:
            return df

        df = df.rename(columns={snp_col: "snp_id", pval_col: "p_value"})
        df["logp"] = -np.log10(df["p_value"])
        df = df.sort_values(by="logp", ascending=False)
        df = df.drop_duplicates(subset=["snp_id"], keep="first")
        return df

    @staticmethod
    def _select_snps_by_rank(gwas_df: pd.DataFrame, start_rank: int, end_rank: int) -> pd.DataFrame:
        if gwas_df.empty:
            return gwas_df
        start_rank = max(1, int(start_rank))
        end_rank = max(start_rank, int(end_rank))
        unique_vals = gwas_df["logp"].unique()
        if unique_vals.size == 0:
            return gwas_df.iloc[0:0]
        start_idx = start_rank - 1
        if start_idx >= len(unique_vals):
            return gwas_df.iloc[0:0]
        end_idx = min(end_rank, len(unique_vals))
        selected_vals = unique_vals[start_idx:end_idx]
        selected = gwas_df[gwas_df["logp"].isin(selected_vals)].copy()
        selected = selected.sort_values(by="logp", ascending=False)
        return selected

    @staticmethod
    def _filter_samples(geno_df: pd.DataFrame, samples: Sequence[str]) -> Tuple[pd.DataFrame, List[str]]:
        if not samples:
            return geno_df, list(geno_df.index)
        present = [s for s in samples if s in geno_df.index]
        return geno_df.loc[present], present

    @staticmethod
    def _filter_snps(geno_df: pd.DataFrame, snps: Sequence[str]) -> Tuple[pd.DataFrame, List[str]]:
        if not snps:
            return geno_df, []
        present = [s for s in snps if s in geno_df.columns]
        return geno_df[present], present

    @staticmethod
    def _corr_sign(geno_df: pd.DataFrame, pheno: pd.Series) -> pd.Series:
        common = geno_df.index.intersection(pheno.index)
        if common.empty:
            return pd.Series(1.0, index=geno_df.columns)
        geno = geno_df.loc[common]
        pheno = pheno.loc[common]
        corr = geno.apply(lambda col: col.corr(pheno))
        sign = np.sign(corr).replace({0: 1}).fillna(1.0)
        return sign

    @staticmethod
    def _score(geno_df: pd.DataFrame, weights: pd.Series) -> pd.Series:
        if geno_df.empty or weights.empty:
            return pd.Series(dtype=float)
        ordered = [s for s in weights.index if s in geno_df.columns]
        if not ordered:
            return pd.Series(dtype=float)
        mat = geno_df[ordered].to_numpy(dtype=float, copy=False)
        w = weights.loc[ordered].to_numpy(dtype=float, copy=False)
        scores = mat.dot(w)
        return pd.Series(scores, index=geno_df.index)

    def _build_weights(
        self,
        gwas_df: pd.DataFrame,
        top_df: pd.DataFrame,
        corr_sign: Optional[pd.Series],
    ) -> pd.Series:
        snps = top_df["snp_id"].tolist()
        if self.weight_mode == "uniform":
            return pd.Series(1.0, index=snps)
        weights = top_df.set_index("snp_id")["logp"]
        if self.weight_mode == "logp":
            return weights
        if self.weight_mode == "corr_logp":
            if corr_sign is None:
                return weights
            signs = corr_sign.reindex(weights.index).fillna(1.0)
            return weights * signs
        raise ValueError(f"Unsupported prs_weight_mode: {self.weight_mode}")

    def _assemble_scores(
        self,
        split: str,
        sample_order: Sequence[str],
        score_map: Dict[str, pd.Series],
        phenotype: Optional[pd.Series],
    ) -> pd.DataFrame:
        if not score_map:
            return pd.DataFrame()
        df = pd.DataFrame({"taxa": list(sample_order)})
        df["split"] = split
        for label, scores in score_map.items():
            df[label] = scores.reindex(sample_order).to_numpy(dtype=float, copy=False)
        if phenotype is not None:
            df[self.phenotype_column] = phenotype.reindex(sample_order).to_numpy(dtype=float, copy=False)
        return df

    def run_fold(
        self,
        file_fold_index: str,
        gwas_pcs_list: Sequence[int],
        *,
        train_samples: Sequence[str],
        test_samples: Sequence[str],
    ) -> None:
        gwas_pcs = [int(pc) for pc in gwas_pcs_list] if gwas_pcs_list else [0]
        pheno_df = None
        if self.include_phenotype or self.weight_mode == "corr_logp":
            pheno_df = self._load_phenotypes()
        if not self.vcf_file_path:
            raise ValueError("vcf_file_path is required for PRS scoring.")
        cohort_geno_raw = self.snp_tool.vcf_to_dataframe(self.vcf_file_path)

        for pc in gwas_pcs:
            file_fold_npc_index = f"{file_fold_index}_{pc}PCs"
            gwas_dir = Path(self.gwas_output_dir) / file_fold_npc_index
            gwas_file = gwas_dir / f"train_{file_fold_npc_index}_{self.phenotype_column}_GWAS_result.txt"

            if not gwas_file.exists():
                time_stamp(f"[PRS] GWAS results missing: {gwas_file}")
                continue

            gwas_df = self._load_gwas_results(gwas_file)
            if gwas_df.empty:
                time_stamp(f"[PRS] No GWAS rows after filtering: {gwas_file}")
                continue

            top_dfs = {}
            union_snps: List[str] = []
            for spec in self.prs_top_specs:
                top_df = self._select_snps_by_rank(gwas_df, spec["start"], spec["end"])
                if top_df.empty:
                    continue
                top_dfs[spec["label"]] = top_df
                union_snps.extend(top_df["snp_id"].tolist())
            union_snps = _dedupe(union_snps)
            if not union_snps:
                time_stamp(f"[PRS] No SNPs selected for {file_fold_npc_index}")
                continue

            train_geno_raw, train_present = self._filter_samples(cohort_geno_raw, train_samples)
            if not train_present or train_geno_raw.empty:
                time_stamp(f"[PRS] No train samples available for {file_fold_npc_index}")
                continue
            test_geno_raw, test_present = self._filter_samples(cohort_geno_raw, test_samples)
            if not test_present or test_geno_raw.empty:
                time_stamp(f"[PRS] No test samples available for {file_fold_npc_index}")
                continue

            train_geno_raw, available_snps = self._filter_snps(train_geno_raw, union_snps)
            if not available_snps:
                time_stamp(f"[PRS] No selected SNPs found in cohort VCF for {file_fold_npc_index}")
                continue
            test_geno_raw, _ = self._filter_snps(test_geno_raw, available_snps)
            train_geno, test_geno = self._prepare_train_test_numeric(train_geno_raw, test_geno_raw)

            corr_sign = None
            if self.weight_mode == "corr_logp" and pheno_df is not None:
                pheno_series = pheno_df[self.phenotype_column]
                corr_sign = self._corr_sign(train_geno, pheno_series)

            train_scores: Dict[str, pd.Series] = {}
            test_scores: Dict[str, pd.Series] = {}
            for label, top_df in top_dfs.items():
                top_snps = [s for s in top_df["snp_id"].tolist() if s in available_snps]
                if not top_snps:
                    continue
                filtered_top_df = top_df[top_df["snp_id"].isin(top_snps)]
                weights = self._build_weights(gwas_df, filtered_top_df, corr_sign)
                weights = weights.reindex(top_snps)
                train_scores[label] = self._score(train_geno, weights)
                test_scores[label] = self._score(test_geno, weights)
            if not train_scores:
                time_stamp(f"[PRS] No usable PRS score columns for {file_fold_npc_index}")
                continue

            out_dir = Path(self.prs_output_dir) / file_fold_npc_index
            out_dir.mkdir(parents=True, exist_ok=True)

            pheno_series = None
            if self.include_phenotype and pheno_df is not None:
                pheno_series = pheno_df[self.phenotype_column]

            train_df = self._assemble_scores("train", train_present, train_scores, pheno_series)
            test_df = self._assemble_scores("test", test_present, test_scores, pheno_series)
            out_scores = pd.concat([train_df, test_df], ignore_index=True)
            out_path = out_dir / f"prs_scores_{file_fold_npc_index}.csv"
            out_scores.to_csv(out_path, index=False)
            time_stamp(f"[PRS] Wrote {out_path}")

            for label, top_df in top_dfs.items():
                top_snps = [s for s in top_df["snp_id"].tolist() if s in available_snps]
                if not top_snps:
                    continue
                filtered_top_df = top_df[top_df["snp_id"].isin(top_snps)]
                weights = self._build_weights(gwas_df, filtered_top_df, corr_sign)
                weights = weights.reindex(top_snps)
                weight_df = filtered_top_df.set_index("snp_id")[["p_value", "logp"]]
                weight_df["weight"] = weights
                label_suffix = str(label).replace("prs_", "")
                weight_path = out_dir / f"prs_weights_{file_fold_npc_index}_{label_suffix}.csv"
                weight_df.reset_index().to_csv(weight_path, index=False)


def run_prs_from_config(
    config_file: str,
    *,
    cv_config_dir: Optional[str] = None,
    folds: Optional[Sequence[str]] = None,
    gwas_pcs_override: str = "",
) -> None:
    config = read_config(config_file)
    prs_tool = PRSCalculator(config)

    output_dir = _resolve_outdir(
        config,
        key="prediction_output_dir",
        default="results",
        ensure_dir=True,
    )
    cv_dir = Path(cv_config_dir) if cv_config_dir else (Path(output_dir) / "cv_configs")
    if not cv_dir.exists():
        raise FileNotFoundError(f"cv_configs not found: {cv_dir}")

    fold_cfg_paths = sorted(cv_dir.glob("fold_*_config.txt"))
    fold_set = set(str(f) for f in folds) if folds else set()

    for fold_cfg_path in fold_cfg_paths:
        fold_cfg = CV_preparation.read_fold_config(fold_cfg_path)
        file_fold_index = str(fold_cfg["file_fold_index"])
        if fold_set and file_fold_index not in fold_set:
            continue
        gwas_pcs_raw = gwas_pcs_override or str(fold_cfg.get("gwas_pcs", "") or "")
        pcs = parse_int_list(gwas_pcs_raw) if gwas_pcs_raw.strip() else [0]

        train_samples = CV_preparation.read_samples_file(Path(fold_cfg["train_samples_file"]))
        test_samples = CV_preparation.read_samples_file(Path(fold_cfg["test_samples_file"]))

        prs_tool.run_fold(
            file_fold_index,
            pcs,
            train_samples=train_samples,
            test_samples=test_samples,
        )


def main():
    parser = argparse.ArgumentParser(description="PRS calculation (post-GWAS).")
    parser.add_argument("--config_file", required=True, help="Path to pipeline config file.")
    parser.add_argument("--cv_config_dir", default="", help="Optional CV config directory.")
    parser.add_argument("--folds", default="", help="Comma-separated fold ids to run.")
    parser.add_argument("--gwas_pcs", default="", help="Override GWAS PCs list (comma-separated).")
    args = parser.parse_args()

    folds = [f.strip() for f in args.folds.split(",") if f.strip()] if args.folds else None
    run_prs_from_config(
        args.config_file,
        cv_config_dir=args.cv_config_dir or None,
        folds=folds,
        gwas_pcs_override=args.gwas_pcs or "",
    )


if __name__ == "__main__":
    main()
