#!/usr/bin/env python3
import os
from typing import Dict

from Prediction_model.prediction_model import SNPPredictionModelTool


class PredictionRocPlotRunner:
    def __init__(self, config: Dict[str, str]):
        self.config = config

    def run(self) -> None:
        tools = SNPPredictionModelTool(self.config)

        probabilities_files = [f for f in os.listdir(tools.model_output_dir) if f.startswith("Probabilities")]
        testy_files = [f for f in os.listdir(tools.model_output_dir) if f.startswith("Testy")]

        if not probabilities_files or not testy_files:
            print(f"No prediction probability outputs found under {tools.model_output_dir}")
            return

        tools.plot_performance_comparison()
