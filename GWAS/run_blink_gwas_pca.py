#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
def main():
    parser = argparse.ArgumentParser(description="Run BLINK GWAS/PCA stage")
    parser.add_argument("config_file", help="Path to the configuration file")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from Prediction_model.gwas_pca import run_gwas_pca_from_config
    from Utils.utils import read_config

    config = read_config(args.config_file)

    vcf_path_raw = (config.get("vcf_file_path") or "").strip()
    pheno_path_raw = (config.get("phenotype_csv_path") or "").strip()
    pheno_col = (config.get("phenotype_column") or "Phenotype").strip()

    if not vcf_path_raw:
        raise ValueError("vcf_file_path is required for GWAS/PCA stage.")
    if not pheno_path_raw:
        raise ValueError("phenotype_csv_path is required for GWAS/PCA stage.")

    run_gwas_pca_from_config(
        config,
        vcf_path=vcf_path_raw,
        pheno_path=pheno_path_raw,
        pheno_col=pheno_col,
    )


if __name__ == "__main__":
    main()
