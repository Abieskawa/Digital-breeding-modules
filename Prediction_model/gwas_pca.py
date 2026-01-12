#!/usr/bin/env python3
"""Compatibility wrapper for GWAS/PCA analysis (moved to GWAS/run_blink_gwas_pca.py)."""
from GWAS.run_blink_gwas_pca import (
    GWAS_PCA,
    run_gwas_pca_from_config,
    run_global_pca_qc_from_config,
)

__all__ = ["GWAS_PCA", "run_gwas_pca_from_config", "run_global_pca_qc_from_config"]
