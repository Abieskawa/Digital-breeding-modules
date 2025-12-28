#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path


def _collect_plot_targets(gwas_dir: str, phenotype_column: str):
    targets = []
    if not os.path.isdir(gwas_dir):
        return targets
    for name in sorted(os.listdir(gwas_dir)):
        current_dir = os.path.join(gwas_dir, name)
        if not os.path.isdir(current_dir):
            continue
        p_value_file = os.path.join(
            current_dir, f"train_{name}_{phenotype_column}_GWAS_result.txt"
        )
        if os.path.isfile(p_value_file):
            targets.append(name)
    return targets


def main():
    parser = argparse.ArgumentParser(description="Plot Manhattan/QQ for BLINK GWAS results")
    parser.add_argument("config_file", help="Path to the configuration file")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root))

    from Prediction_model.prediction_model import (
        SNPPredictionModelTool,
        parallel_process,
    )
    from Utils.utils import read_config

    config = read_config(args.config_file)
    tools = SNPPredictionModelTool(config)

    targets = _collect_plot_targets(tools.gwas_output_dir, tools.phenotype_column)
    if not targets:
        print(f"No GWAS results found under {tools.gwas_output_dir}")
        return

    parallel_process(
        targets,
        tools.plot_manhattan_and_qq,
        max_workers=tools.num_threads,
    )


if __name__ == "__main__":
    main()
