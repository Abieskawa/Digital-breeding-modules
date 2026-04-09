import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path

import pandas as pd

sys.modules.setdefault("matplotlib", types.ModuleType("matplotlib"))
sys.modules.setdefault("matplotlib.pyplot", types.ModuleType("matplotlib.pyplot"))
sys.modules.setdefault("seaborn", types.ModuleType("seaborn"))

REPO_ROOT = Path(__file__).resolve().parents[1]
EVAL_DIR = REPO_ROOT / "Evaluation"
if "Evaluation" not in sys.modules:
    eval_pkg = types.ModuleType("Evaluation")
    eval_pkg.__path__ = [str(EVAL_DIR)]
    sys.modules["Evaluation"] = eval_pkg

mapping_spec = importlib.util.spec_from_file_location("Evaluation.mapping_base", EVAL_DIR / "mapping_base.py")
mapping_module = importlib.util.module_from_spec(mapping_spec)
sys.modules["Evaluation.mapping_base"] = mapping_module
mapping_spec.loader.exec_module(mapping_module)

contam_spec = importlib.util.spec_from_file_location("Evaluation.contamination", EVAL_DIR / "contamination.py")
contam_module = importlib.util.module_from_spec(contam_spec)
sys.modules["Evaluation.contamination"] = contam_module
contam_spec.loader.exec_module(contam_module)

_build_non_target_summaries = contam_module._build_non_target_summaries
_resolve_kraken_db_path = contam_module._resolve_kraken_db_path


class ContaminationHelpersTest(unittest.TestCase):
    def test_build_non_target_summaries_filters_target_species(self):
        df = pd.DataFrame(
            [
                {
                    "name": "Oreochromis niloticus",
                    "sample": "S1",
                    "name_freq_pct": 92.0,
                    "n_records": 920,
                    "total_records": 1000,
                },
                {
                    "name": "Homo sapiens",
                    "sample": "S1",
                    "name_freq_pct": 4.6,
                    "n_records": 46,
                    "total_records": 1000,
                },
                {
                    "name": "Homo sapiens",
                    "sample": "S2",
                    "name_freq_pct": 3.2,
                    "n_records": 32,
                    "total_records": 1000,
                },
            ]
        )

        per_sample, cohort = _build_non_target_summaries(df, "Oreochromis_niloticus")

        self.assertEqual(list(per_sample["top_non_target_species"]), ["Homo sapiens", "Homo sapiens"])
        self.assertEqual(list(cohort["name"]), ["Homo sapiens"])

    def test_resolve_kraken_db_path_uses_existing_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            resolved = _resolve_kraken_db_path(tmpdir)
            self.assertEqual(resolved, Path(tmpdir).resolve())


if __name__ == "__main__":
    unittest.main()
