import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from DNA.variant_calling import VariantCalling


class VariantCallingDeepVariantTest(unittest.TestCase):
    def _make_args(self, tmpdir: str, *, mode: str) -> SimpleNamespace:
        tmp_path = Path(tmpdir)
        ref_path = tmp_path / "ref.fa"
        ref_path.write_text(">chr1\nACGT\n")
        align_dir = tmp_path / "align"
        variant_dir = tmp_path / "variant"
        align_dir.mkdir()
        variant_dir.mkdir()
        return SimpleNamespace(
            variant_outdir=variant_dir,
            align_outdir=align_dir,
            ref_fasta=ref_path,
            threads=4,
            deepvariant_mode=mode,
            deepvariant_use_gpu=(mode == "gpu"),
            deepvariant_num_shards=2,
            step3_max_concurrent_samples=1,
            step3_gpu_jobs=1,
            deepvariant_extra_args="",
            capture_bed=None,
            hwe_p="1e-5",
            geno_missing="0.1",
            maf=0.05,
            ld_prune=False,
            ld_method="indep",
            ld_window=50,
            ld_step=5,
            ld_param=2.0,
        )

    def test_gpu_mode_selects_gpu_backend_without_runtime_probe(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            vc = VariantCalling(self._make_args(tmpdir, mode="gpu"))
            vc._local_deepvariant_available = lambda: True
            self.assertEqual(vc._select_deepvariant_backend(), "local_gpu")

    def test_gpu_backend_command_does_not_add_use_gpu_flag(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            vc = VariantCalling(self._make_args(tmpdir, mode="gpu"))
            bam_path = vc.align_outdir / "sample.bam"
            bam_path.write_text("placeholder")
            cmd, _, _ = vc._build_deepvariant_cmd("sample", bam_path, "local_gpu")
            self.assertNotIn("--use_gpu", cmd)

    def test_missing_local_binary_raises_clear_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            vc = VariantCalling(self._make_args(tmpdir, mode="cpu"))
            vc._local_deepvariant_available = lambda: False
            with self.assertRaises(RuntimeError):
                vc._select_deepvariant_backend()


if __name__ == "__main__":
    unittest.main()
