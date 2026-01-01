#!/usr/bin/env python3
import os
from typing import Dict, List

from Prediction_model.prediction_model import SNPPredictionModelTool
from Utils.utils import parallel_process


def _collect_plot_targets(gwas_dir: str, phenotype_column: str) -> List[str]:
    targets: List[str] = []
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


class GWASManhattanPlotRunner:
    def __init__(self, config: Dict[str, str]):
        self.config = config

    def run(self) -> None:
        tools = SNPPredictionModelTool(self.config)

        targets = _collect_plot_targets(tools.gwas_output_dir, tools.phenotype_column)
        if not targets:
            print(f"No GWAS results found under {tools.gwas_output_dir}")
            return

        parallel_process(
            targets,
            tools.plot_manhattan_and_qq,
            max_workers=tools.num_threads,
        )
