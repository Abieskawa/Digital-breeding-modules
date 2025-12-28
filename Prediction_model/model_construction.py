#!/usr/bin/env python3
from pathlib import Path

from Prediction_model.prediction_base import _PredictionBase


class ModelConstruction(_PredictionBase):
    def __init__(self, args, config_path: Path):
        super().__init__(args, config_path)

    def run(self):
        self._run_prediction_model(stage="model")

    def run_eval(self):
        self._run_eval_script("Evaluation/prediction_roc_plot.py")
