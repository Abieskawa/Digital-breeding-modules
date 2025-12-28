#!/usr/bin/env python3
"""
Backward-compatible exports for prediction step classes.
"""
from Prediction_model.gwas_pca import GWASPCA
from Prediction_model.model_construction import ModelConstruction

__all__ = ["GWASPCA", "ModelConstruction"]
