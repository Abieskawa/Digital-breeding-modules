#!/usr/bin/env python3
from Prediction_model.prediction_base import _PredictionBase


class GWASPCA(_PredictionBase):
    def run(self):
        self._run_prediction_model(stage="gwas")

    def run_eval(self):
        self._run_eval_script("Evaluation/gwas_manhattan_plot.py")
