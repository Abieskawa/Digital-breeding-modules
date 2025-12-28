#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Run BLINK GWAS/PCA stage")
    parser.add_argument("config_file", help="Path to the configuration file")
    args, passthrough = parser.parse_known_args()

    repo_root = Path(__file__).resolve().parents[1]
    model_script = repo_root / "Prediction_model" / "prediction_model.py"
    if not model_script.is_file():
        raise FileNotFoundError(f"Missing prediction script: {model_script}")

    cmd = [
        sys.executable,
        str(model_script),
        str(args.config_file),
        "--stage",
        "gwas",
        *passthrough,
    ]
    subprocess.run(cmd, check=True, cwd=repo_root)


if __name__ == "__main__":
    main()
