#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Plot ROC/AUC curves for prediction results")
    parser.add_argument("config_file", help="Path to the configuration file")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root))

    from Prediction_model.prediction_model import SNPPredictionModelTool
    from Utils.utils import read_config

    config = read_config(args.config_file)
    tools = SNPPredictionModelTool(config)

    probabilities_files = [f for f in os.listdir(tools.model_output_dir) if f.startswith("Probabilities")]
    testy_files = [f for f in os.listdir(tools.model_output_dir) if f.startswith("Testy")]

    if not probabilities_files or not testy_files:
        print(f"No prediction probability outputs found under {tools.model_output_dir}")
        return

    tools.plot_performance_comparison()


if __name__ == "__main__":
    main()
