#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys


class _PredictionBase:
    def __init__(self, args, config_path: Path):
        if not config_path:
            raise ValueError("config_path is required for GWAS/prediction steps.")
        self.args = args
        self.config_path = Path(config_path).resolve()
        self.repo_root = Path(__file__).resolve().parents[1]
        self.prediction_python = getattr(args, "prediction_python", sys.executable)

    def _run_prediction_model(self, stage: str):
        script = self.repo_root / "Prediction_model" / "prediction_model.py"
        cmd = [self.prediction_python, str(script), str(self.config_path)]
        if stage:
            cmd.extend(["--stage", stage])
        subprocess.run(cmd, check=True, cwd=self.repo_root)

    def _run_eval_script(self, rel_script: str):
        script = self.repo_root / rel_script
        subprocess.run(
            [self.prediction_python, str(script), str(self.config_path)],
            check=True,
            cwd=self.repo_root,
        )
