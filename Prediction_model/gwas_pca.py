#!/usr/bin/env python3
"""
Step 2: PCA + GWAS (BLINK)

Inputs:
  - VCF (config: vcf_file_path)
  - Numeric genotype NPZ (from Step 1 pca_numeric)

Outputs (same layout as reference script):
  output_dir/GWAS_dir/{file_fold_index}_{pc}PCs/
    train_{file_fold_index}_{pc}PCs.txt      (phenotype)
    train_{file_fold_index}_{pc}PCs.vcf      (subset VCF for BLINK)
    train_{file_fold_index}_{pc}PCs.cov      (covariates, if pc>0 and use_pcs)
    blink_run_train_{...}.log
    train_{...}_{phenotype_column}_GWAS_result.txt (BLINK output)

Notes:
  - This script keeps BLINK command style from the reference.
  - It fixes the VCF "format issue" by writing an UNCOMPRESSED .vcf file (bcftools -Ov),
    so downstream pandas-based VCF reading works reliably.
  - "PC number extracted for GWAS" is separated:
      * gwas_n_pcs_list (covariates used in GWAS; loop values)
      * pca_n_components (how many PCs to compute in PCA; must be >= max(gwas pcs), >= 15)
"""
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from Utils.utils import _resolve_outdir, parse_int_list, extract_vcf_samples


def _load_numeric_npz(npz_path: str) -> pd.DataFrame:
    data = np.load(npz_path, allow_pickle=True)
    genotype = data["genotype"]
    samples = data["samples"].astype(str)
    snps = data["snps"].astype(str)
    return pd.DataFrame(genotype, index=samples, columns=snps)


class GWAS_PCA:
    def __init__(self, config: dict):
        self.config = config
        self.species_name = self.config.get('species_name')
        self.phenotype_csv_path = self.config.get('phenotype_csv_path')
        self.vcf_file_path = self.config.get('vcf_file_path')
        self.chromosome_csv_path = self.config.get('chromosome_csv_path')
        self.phenotype_column = self.config.get('phenotype_column', 'Phenotype')

        self.use_pcs = self.config.get('use_pcs', 'True').lower() == 'true'

        # PC covariates for GWAS
        self.gwas_n_pcs_list = parse_int_list(self.config.get('gwas_n_pcs_list', '')) or \
                              parse_int_list(self.config.get('gwas_n_pcs', '')) or \
                              parse_int_list(self.config.get('n_pcs_list', '0'))

        # PCA extraction count (separate parameter)
        # Default: max(max(gwas pcs), 15)
        raw_pca_n = self.config.get('pca_n_components', '').strip()
        if raw_pca_n:
            try:
                self.pca_n_components = int(raw_pca_n)
            except Exception:
                self.pca_n_components = 15
        else:
            self.pca_n_components = 15

        self.output_dir = _resolve_outdir(
            self.config,
            key="prediction_output_dir",
            default="results",
            resolve=True,
            ensure_dir=True,
        )
        self.pca_scree_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="PCA_Scree_Plot",
            resolve=True,
            ensure_dir=True,
        )
        self.gwas_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="GWAS_dir",
            resolve=True,
            ensure_dir=True,
        )

        # BLINK binary path (keep reference default, allow override)
        self.blink_bin = self.config.get('blink_bin', '/home/abieskawa/tools/blink/blink1_5')

        # Chromosome mapping
        self.load_chromosome_recode()

    def load_chromosome_recode(self):
        self.chromosome_mapping = pd.read_csv(self.chromosome_csv_path)
        self.chr_name_to_recoded = dict(zip(self.chromosome_mapping['OriginalName'], self.chromosome_mapping['RecodedNumber']))
        self.recoded_to_label = dict(zip(self.chromosome_mapping['RecodedNumber'].astype(str), self.chromosome_mapping['Label']))

    def perform_pca(self, genotype_df: pd.DataFrame, n_components: int) -> Tuple[pd.DataFrame, np.ndarray]:
        scaler = StandardScaler()
        genotype_scaled = scaler.fit_transform(genotype_df)

        pca = PCA(n_components=n_components)
        principal_components = pca.fit_transform(genotype_scaled)
        explained_variance_ratio = pca.explained_variance_ratio_

        pc_columns = [f'PC{i+1}' for i in range(n_components)]
        principal_df = pd.DataFrame(data=principal_components, columns=pc_columns)
        principal_df.insert(0, 'taxa', genotype_df.index)
        return principal_df, explained_variance_ratio

    def blink_file_generator(self, train_samples: List[str], file_fold_npc_index: str, principal_df: pd.DataFrame, n_pcs: int):
        """
        Generate input files for BLINK GWAS analysis.
        Adapted from the reference blink_file_generator, with one key fix:
          - output_vcf is written as UNCOMPRESSED VCF (-Ov) to avoid downstream format issues.
        """
        train_samples_str = ','.join(train_samples)

        os.makedirs(os.path.join(self.gwas_output_dir, f'{file_fold_npc_index}'), exist_ok=True)
        current_dir = os.path.join(self.gwas_output_dir, f'{file_fold_npc_index}')

        # phenotype
        phenotype_df = pd.read_csv(self.phenotype_csv_path, sep=None, engine='python')
        phenotype_df = phenotype_df[['taxa', self.phenotype_column]]
        filtered_phenotype_df = phenotype_df[phenotype_df['taxa'].isin(train_samples)]
        phenotype_file = os.path.join(current_dir, f'train_{file_fold_npc_index}.txt')
        filtered_phenotype_df.to_csv(phenotype_file, sep='\t', index=False)

        # genotype vcf for BLINK
        vcf_file = self.vcf_file_path
        output_vcf = os.path.join(current_dir, f'train_{file_fold_npc_index}.vcf')

        mapping_file = os.path.join(current_dir, 'chr_mapping.txt')
        if not os.path.exists(mapping_file):
            self.chromosome_mapping[['OriginalName', 'RecodedNumber']].to_csv(
                mapping_file, sep='\t', header=False, index=False
            )

        original_chromosomes_to_keep = self.chromosome_mapping['OriginalName'].tolist()
        regions_str = ','.join(map(str, original_chromosomes_to_keep))

        # bgzip + index original vcf for fast bcftools access
        compressed_vcf = vcf_file + ".gz"
        if not os.path.exists(compressed_vcf):
            subprocess.run(f"bgzip -c {vcf_file} > {compressed_vcf}", shell=True, check=True)
            subprocess.run(f"bcftools index {compressed_vcf}", shell=True, check=True)
        else:
            # Ensure index exists
            if not (os.path.exists(compressed_vcf + ".csi") or os.path.exists(compressed_vcf + ".tbi")):
                subprocess.run(f"bcftools index {compressed_vcf}", shell=True, check=True)

        command = (
            f'bcftools view -s "{train_samples_str}" -r "{regions_str}" "{compressed_vcf}" | '
            f'bcftools annotate --rename-chrs "{mapping_file}" -Ov -o "{output_vcf}"'
        )
        subprocess.run(command, shell=True, check=True, executable="/bin/bash")

        # covariates
        if self.use_pcs and n_pcs > 0:
            principal_use = principal_df.iloc[:, :(n_pcs + 1)]
            filtered_principal_df = principal_use[principal_use['taxa'].isin(train_samples)]
            principal_file = os.path.join(current_dir, f'train_{file_fold_npc_index}.cov')
            filtered_principal_df.to_csv(principal_file, sep='\t', index=False)

    def run_gwas(self, file_fold_npc_index: str):
        current_dir = os.path.join(self.gwas_output_dir, f'{file_fold_npc_index}')
        os.makedirs(current_dir, exist_ok=True)
        tag = f"train_{file_fold_npc_index}"

        blink_command = (
            f'{self.blink_bin} '
            f'--gwas '
            f'--file {tag} '
            f'--vcf '
            f'--out {tag}'
        )

        log_path = os.path.join(current_dir, f'blink_run_{tag}.log')
        with open(log_path, 'w') as log_file:
            subprocess.run(
                blink_command,
                shell=True,
                check=True,
                cwd=current_dir,
                stdout=log_file,
                stderr=subprocess.STDOUT
            )

    def run_fold(self, file_fold_index: str, genotype_npz: str, gwas_pcs_override: Optional[str] = None, pca_components_override: Optional[int] = None):
        train_numeric = _load_numeric_npz(genotype_npz)
        train_samples = train_numeric.index.astype(str).tolist()

        gwas_pcs_list = parse_int_list(gwas_pcs_override) if (gwas_pcs_override and gwas_pcs_override.strip()) else self.gwas_n_pcs_list
        if not gwas_pcs_list:
            gwas_pcs_list = [0]

        # Determine PCA extraction count
        need_pca = self.use_pcs and max(gwas_pcs_list) > 0
        if need_pca:
            n_components = max(self.pca_n_components, 15, max(gwas_pcs_list))
            if pca_components_override is not None:
                n_components = max(n_components, int(pca_components_override))
            principal_df, explained_variance_ratio = self.perform_pca(train_numeric, n_components=n_components)
            # Save explained variance (optional)
            try:
                out_csv = os.path.join(self.pca_scree_output_dir, f'explained_variance_{file_fold_index}.csv')
                pd.DataFrame({"explained_variance_ratio": explained_variance_ratio}).to_csv(out_csv, index=False)
            except Exception:
                pass
        else:
            principal_df = pd.DataFrame({"taxa": train_numeric.index})

        for pc in gwas_pcs_list:
            file_fold_npc_index = f"{file_fold_index}_{pc}PCs"
            self.blink_file_generator(train_samples, file_fold_npc_index, principal_df, pc)
            self.run_gwas(file_fold_npc_index)


def run_gwas_pca_from_config(
    config: dict,
    *,
    vcf_path: Path,
    pheno_path: Path,
    pheno_col: str,
    gwas_pcs_override: str = "",
    pca_components_override: Optional[int] = None,
) -> None:
    from Prediction_model.snp_numeric_transformer import SNP_numerical, _save_numeric_npz

    vcf_path = Path(vcf_path).expanduser().resolve()
    pheno_path = Path(pheno_path).expanduser().resolve()

    if not vcf_path.exists():
        raise FileNotFoundError(f"vcf_file_path not found: {vcf_path}")
    if not pheno_path.exists():
        raise FileNotFoundError(f"phenotype_csv_path not found: {pheno_path}")

    pheno = pd.read_csv(pheno_path, sep=None, engine="python")
    if "taxa" not in pheno.columns:
        raise ValueError("phenotype_csv_path must contain a 'taxa' column.")
    if pheno_col not in pheno.columns:
        raise ValueError(f"phenotype column not found: {pheno_col}")

    pheno = pheno[["taxa", pheno_col]].dropna().set_index("taxa")
    vcf_samples = extract_vcf_samples(vcf_path)
    common = [s for s in vcf_samples if s in pheno.index]
    if not common:
        raise ValueError("No overlapping samples between VCF and phenotype CSV.")

    prediction_outdir = _resolve_outdir(config, key="prediction_output_dir", default="results", resolve=True)
    tmp_root = Path(prediction_outdir) / ".tmp" / "prediction_gwas"
    tmp_root.mkdir(parents=True, exist_ok=True)
    npz_path = tmp_root / "train_numeric.npz"

    config = dict(config)
    config.update(
        vcf_file_path=str(vcf_path),
        phenotype_csv_path=str(pheno_path),
        phenotype_column=pheno_col,
    )

    snp_tool = SNP_numerical(config)
    gwas_tool = GWAS_PCA(config)

    genotype_df = snp_tool.vcf_to_dataframe(str(vcf_path))
    present = [s for s in common if s in genotype_df.index]
    if not present:
        raise ValueError("No requested train samples were found in the VCF.")

    genotype_numeric = snp_tool.convert_genotypes_to_numeric(genotype_df.loc[present])
    _save_numeric_npz(genotype_numeric, str(npz_path))

    pca_components = pca_components_override
    if pca_components is None:
        raw_pca_components = (config.get("pca_n_components") or "").strip()
        if raw_pca_components:
            try:
                pca_components = int(raw_pca_components)
            except ValueError:
                pca_components = None

    gwas_tool.run_fold(
        "all",
        str(npz_path),
        gwas_pcs_override=gwas_pcs_override,
        pca_components_override=pca_components,
    )


class PredictionGWASCVStep:
    def __init__(self, fold_configs: List[Path], max_workers: int, worker: Callable[[Dict[str, object]], str]):
        self.fold_configs = fold_configs
        self.max_workers = max_workers
        self.worker = worker

    def run(self) -> List[str]:
        finished: List[str] = []
        with ProcessPoolExecutor(max_workers=int(self.max_workers)) as ex:
            futs = [
                ex.submit(self.worker, {"fold_config": str(cfg)})
                for cfg in self.fold_configs
            ]
            for fut in as_completed(futs):
                finished.append(fut.result())
        return finished
