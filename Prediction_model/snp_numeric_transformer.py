#!/usr/bin/env python3
"""
Step 1: SNP convert to numerical (and model input preparation)

This script is split out from the original monolithic prediction_model.py logic.

It supports TWO modes that match your intended pipeline flow:

1) pca_numeric
   Input : VCF + train sample list
   Output: numeric genotype matrix saved as NPZ (for PCA/GWAS)

2) model_prep
   Input : GWAS outputs (from step 2) + test sample list
   Output: model-ready train/test CSVs (selected top SNPs from GWAS)

Unrelated parts (PCA, GWAS running, model training, CV loops) are NOT implemented here.
"""
import gzip
import os
from typing import List

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from Utils.utils import _resolve_outdir, parse_int_list


def _open_text_maybe_gz(path: str):
    """
    Open plain text or bgzip/gzip VCF transparently.
    """
    with open(path, "rb") as f:
        sig = f.read(2)
    if sig == b"\x1f\x8b":  # gzip magic
        return gzip.open(path, "rt")
    return open(path, "rt")


class SNP_numerical:
    def __init__(self, config: dict):
        self.config = config
        self.species_name = self.config.get('species_name')
        self.phenotype_csv_path = self.config.get('phenotype_csv_path')
        self.vcf_file_path = self.config.get('vcf_file_path')
        self.phenotype_column = self.config.get('phenotype_column', 'Phenotype')

        self.output_dir = _resolve_outdir(
            self.config,
            key="prediction_output_dir",
            default="results",
            ensure_dir=True,
        )

        # Keep same directory naming as the reference script
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

        self.top_n_snps_list = parse_int_list(self.config.get('top_n_snps_list', '10'))

        # GWAS PC list (can be independent from PCA extraction count)
        self.gwas_n_pcs_list = parse_int_list(self.config.get('gwas_n_pcs_list', '')) or \
                              parse_int_list(self.config.get('gwas_n_pcs', '')) or \
                              parse_int_list(self.config.get('n_pcs_list', '0'))

    def get_vcf_header(self, vcf_file: str) -> List[str]:
        with _open_text_maybe_gz(vcf_file) as f:
            for line in f:
                if line.startswith('#CHROM'):
                    return line.strip().split('\t')
        raise ValueError("VCF header not found.")

    def vcf_to_dataframe(self, vcf_file: str) -> pd.DataFrame:
        header = self.get_vcf_header(vcf_file)
        # pandas can read gz if compression='infer' and extension is .gz, but we also allow magic detection above
        vcf_reader = pd.read_csv(vcf_file, sep='\t', comment='#', header=None, dtype={0: str}, compression='infer')
        vcf_reader.columns = header
        sample_columns = header[9:]
        genotype_df = vcf_reader[['ID'] + sample_columns]
        genotype_df = genotype_df.transpose()
        genotype_df.columns = genotype_df.iloc[0]
        genotype_df = genotype_df[1:]
        return genotype_df

    def convert_genotypes_to_numeric(self, genotype_df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert genotype strings to numeric values and fill NA using median imputation.
        Robust to VCF sample fields like "0/1:DP:..." by using only the GT part.
        """
        def _to_gt(x):
            if pd.isna(x):
                return np.nan
            s = str(x)
            if ":" in s:
                s = s.split(":", 1)[0]
            return s

        def _map_gt(gt):
            if pd.isna(gt):
                return np.nan
            if gt in ('0/0', '0|0'):
                return 1
            if gt in ('0/1', '1/0', '0|1', '1|0'):
                return 2
            if gt in ('1/1', '1|1'):
                return 3
            if gt in ('./.', '.|.'):
                return np.nan
            # Keep as-is (will fail cast if truly non-numeric); matches original behavior
            return gt

        genotype_numeric_df = genotype_df.apply(lambda col: col.map(lambda x: _map_gt(_to_gt(x))))
        genotype_numeric_df = genotype_numeric_df.astype(float)

        imputer = SimpleImputer(strategy='median')
        genotype_imputed = pd.DataFrame(
            imputer.fit_transform(genotype_numeric_df),
            columns=genotype_numeric_df.columns,
            index=genotype_numeric_df.index
        )
        return genotype_imputed

    def select_top_snps(self, gwas_results_file: str, top_n_snps: int = 10) -> List[str]:
        gwas_df = pd.read_csv(gwas_results_file, sep='\t')
        gwas_df = gwas_df[gwas_df['p_value'] > 0]
        gwas_df['-logp'] = -np.log10(gwas_df['p_value'])
        sorted_df = gwas_df.sort_values(by='-logp', ascending=False)
        top_n_unique_values = sorted_df['-logp'].unique()[:top_n_snps]
        top_df = sorted_df[sorted_df['-logp'].isin(top_n_unique_values)]
        top_snps = top_df['taxa'].tolist()
        return top_snps

    def select_markers_and_prepare_data_for_model_construction(
        self,
        file_fold_npc_index: str,
        file_fold_npc_nsnp_index: str,
        test_df: pd.DataFrame,
        top_n_snps: int = 10
    ) -> List[str]:
        """
        This is directly adapted from the reference script
        (select_markers_and_prepare_data_for_model_construction),
        except it is now located in Step 1 because it is "SNP -> numeric + formatting".
        """
        current_dir = os.path.join(self.gwas_output_dir, file_fold_npc_index)
        p_value_file = os.path.join(current_dir, f'train_{file_fold_npc_index}_{self.phenotype_column}_GWAS_result.txt')
        vcf_file = os.path.join(current_dir, f'train_{file_fold_npc_index}.vcf')

        top_snps = self.select_top_snps(p_value_file, top_n_snps)

        phenotype_df = pd.read_csv(self.phenotype_csv_path, sep=None, engine='python')
        phenotype_df = phenotype_df[['taxa', self.phenotype_column]]
        phenotype_df.set_index('taxa', inplace=True)

        train_genotype_df = self.vcf_to_dataframe(vcf_file)
        train_selected_snps = train_genotype_df[top_snps]
        train_selected_snps = self.convert_genotypes_to_numeric(train_selected_snps)

        test_selected_snps = test_df[top_snps]
        test_selected_snps = self.convert_genotypes_to_numeric(test_selected_snps)

        train_selected_snps = train_selected_snps.join(phenotype_df, how='inner')
        train_selected_snps.rename(columns={self.phenotype_column: 'Phenotype'}, inplace=True)

        test_selected_snps = test_selected_snps.join(phenotype_df, how='inner')
        test_selected_snps.rename(columns={self.phenotype_column: 'Phenotype'}, inplace=True)

        save_dir = os.path.join(self.model_output_dir, f'{file_fold_npc_nsnp_index}')
        os.makedirs(save_dir, exist_ok=True)

        train_csv_path = os.path.join(save_dir, f'train_selected_snps_{file_fold_npc_nsnp_index}.csv')
        train_selected_snps.to_csv(train_csv_path, index_label='taxa')

        test_csv_path = os.path.join(save_dir, f'test_selected_snps_{file_fold_npc_nsnp_index}.csv')
        test_selected_snps.to_csv(test_csv_path, index_label='taxa')

        return top_snps


def _save_numeric_npz(df_numeric: pd.DataFrame, out_npz: str) -> None:
    os.makedirs(os.path.dirname(out_npz) or ".", exist_ok=True)
    np.savez_compressed(
        out_npz,
        genotype=df_numeric.to_numpy(dtype=np.float32, copy=False),
        samples=np.array(df_numeric.index.astype(str)),
        snps=np.array(df_numeric.columns.astype(str)),
    )
