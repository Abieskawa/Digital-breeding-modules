#!/usr/bin/env python3
from typing import Dict, Set, Optional, List
from Utils.utils import time_stamp, read_config
from Utils.bed_generator import build_capture_bed
from pathlib import Path
import sys
import shutil
import Preprocess.fastp_preprocessing as fastp_preprocessing
from DNA.ref_index import RefIndex
from DNA.dna_alignment import DNAAlignment
from DNA.variant_calling import VariantCalling
from Prediction_model.gwas_pca import GWASPCA
from Prediction_model.model_construction import ModelConstruction
import argparse
from Evaluation.fastq_qc import FastQCEvaluator
from Evaluation.contamination import ContaminationEvaluator
from Evaluation.evaluate_mapping import MappingYieldEvaluator
from Evaluation.variant_circos import VariantCircosEvaluator

def _coerce_int(name: str, raw: str, default: int) -> int:
    """
    Convert raw to int, quietly falling back to default when the field is blank/missing.
    Only warn when a non-empty value fails to parse.
    """
    if raw is None:
        return default
    s = str(raw).strip()
    if s == "":
        return default
    try:
        return int(s)
    except (TypeError, ValueError):
        print(f"Warning: Invalid integer for {name}, using default {default}.")
        return default

def _parse_steps(raw: str) -> Set[str]:
    """
    Accepts: '1', '2', '3', '1,2,3', '1-3'
    Returns a set like {'1','2','3'}
    """
    s = (raw or "").strip().lower()
    parts: Set[str] = set()
    for tok in s.replace(" ", "").split(","):
        if not tok:
            continue
        if "-" in tok:
            a, b = tok.split("-", 1)
            if a.isdigit() and b.isdigit():
                for i in range(int(a), int(b) + 1):
                    parts.add(str(i))
        elif tok.isdigit():
            parts.add(tok)
    return {p for p in parts}

def _coerce_bool(raw, default: bool) -> bool:
    """
    Convert common truthy/falsey strings to bool; fallback to default on blank.
    """
    if raw is None:
        return default
    s = str(raw).strip().lower()
    if s == "":
        return default
    if s in {"1", "true", "yes", "y", "on"}:
        return True
    if s in {"0", "false", "no", "n", "off"}:
        return False
    return default

def _parse_list(raw: str) -> Set[str]:
    """
    Parse comma/semicolon-separated list into a set of non-empty strings.
    """
    if raw is None:
        return set()
    cleaned = raw.replace(";", ",")
    vals = [tok.strip() for tok in cleaned.split(",")]
    return {v for v in vals if v}

def _parse_list_ordered(raw: str) -> List[str]:
    """
    Parse comma/semicolon-separated list into a de-duplicated list (preserve order).
    """
    if raw is None:
        return []
    cleaned = raw.replace(";", ",")
    vals = [tok.strip() for tok in cleaned.split(",")]
    ordered: List[str] = []
    seen: Set[str] = set()
    for v in vals:
        if not v or v in seen:
            continue
        seen.add(v)
        ordered.append(v)
    return ordered


class PredictionPipeline(object):
    def __init__(self, configure: Dict[str, str], config_path: Optional[Path] = None):
        self.config_path = config_path
        step_raw = configure.get("step", "0")
        self.steps = _parse_steps(step_raw)
        eva_step_raw = configure.get("eva_step", "0")
        self.eva_steps = _parse_steps(eva_step_raw)
        self._ev_qc = None
        self._ev_contam = None
        self._ev_align = None
        self._ev_var = None
        self._gwas_runner = None
        self._model_runner = None

        # Paths
        self.outdir = Path(configure.get("output_dir","DGBreeding")).resolve()
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.input_dir = Path(configure.get("input_dir",".")).resolve()
        if not self.input_dir.exists():
            raise FileNotFoundError(f"input_dir not found: {self.input_dir}")

        ref_raw = (configure.get("ref_fasta") or "").strip()
        self.ref_fasta = Path(ref_raw).expanduser().resolve() if ref_raw else None

        need_reference = bool({"1", "2", "3"} & self.steps or "3" in self.eva_steps)
        if need_reference:
            if not self.ref_fasta:
                raise ValueError(
                    "ref_fasta is required for steps 1/2/3 or evaluation step 3 but is missing in the config."
                )
            if not (self.ref_fasta.exists() and self.ref_fasta.is_file()):
                raise FileNotFoundError(f"ref_fasta not found or not a file: {self.ref_fasta}")

        gff_raw = (configure.get("gff") or "").strip()
        self.gff = Path(gff_raw).resolve() if gff_raw else None

        # Resources
        self.threads = _coerce_int("threads", configure.get("threads", ""), 96)
        self.mem = configure.get("mem","8G")

        # fastp
        self.fastp_threads = _coerce_int("fastp_threads", configure.get("fastp_threads", ""), self.threads)
        self.fastp_quality = _coerce_int("fastp_quality", configure.get("fastp_quality", ""), 20)
        self.fastp_minlength = _coerce_int("fastp_length", configure.get("fastp_length", ""), 30)
        self.fastp_average_qual = _coerce_int("fastp_average_qual", configure.get("fastp_average_qual", "20"), 20)
        self.fastp_trim_front = configure.get("fastp_trim_front", "10")
        self.fastp_outdir = self.outdir / "00_Preprocessed_DNA"
        self.fastp_report_dir = self.fastp_outdir / "fastp_reports"
        self.fastp_outdir.mkdir(parents=True, exist_ok=True)
        self.fastp_report_dir.mkdir(parents=True, exist_ok=True)

        # Alignment context
        self.population = configure.get("population", "diversity_panel")
        self.seq_platform = (configure.get("seq_platform", "HiSeq")).strip().lower()
        self.optical_distance = _coerce_int("optical_distance", configure.get("optical_distance", ""), 100)
        self.mark_duplicates = _coerce_bool(configure.get("mark_duplicates", "true"), True)
        self.filter_proper_pairs = _coerce_bool(configure.get("filter_proper_pairs", "true"), True)
        self.align_outdir = self.outdir / "01_DNAseq_alignment"
        self.align_outdir.mkdir(parents=True, exist_ok=True)
        # Optional: skip specific samples at alignment/variant stages
        self.skip_align_samples = _parse_list(configure.get("skip_alignment_samples", ""))
        # Default variant skip list to alignment skip list if not provided separately
        raw_variant_skip = configure.get("skip_variant_samples", "")
        self.skip_variant_samples = _parse_list(raw_variant_skip) or set(self.skip_align_samples)

        # Variant calling context (directory only; tool options use VariantCalling defaults
        # unless you add them to configure.txt)
        self.variant_outdir = self.outdir / "02_Variant_Calling"

        # Optional VC knobs (leave unset if not present; VariantCalling has sane defaults)
        self.deepvariant_image = configure.get("deepvariant_image", "google/deepvariant:1.9.0")
        self.model_type        = configure.get("model_type", "WGS")  # or WES
        self.capture_bed       = configure.get("capture_bed", "") or None
        self.generate_bed      = _coerce_bool(configure.get("generate_bed", ""), False)
        self.bed_top_n         = _coerce_int("bed_top_n", configure.get("bed_top_n", ""), 0)
        self.bed_chroms        = _parse_list_ordered(configure.get("bed_chroms", ""))
        self.bed_require_genes = _coerce_bool(configure.get("bed_require_genes", ""), False)
        self.bed_gene_biotype  = (configure.get("bed_gene_biotype", "") or "").strip() or None
        bed_output_raw = (configure.get("bed_output", "") or "").strip()
        self.bed_output = Path(bed_output_raw).expanduser().resolve() if bed_output_raw else None
        self.hwe_p             = configure.get("hwe_p", "1e-5")
        self.geno_missing      = configure.get("geno_missing", "0.10")
        ld_prune_raw = configure.get("ld_prune", "")
        if ld_prune_raw == "":
            ld_prune_raw = configure.get("ld_pruning", "")
        self.enable_ld_pruning = _coerce_bool(ld_prune_raw, True)
        ld_method_raw = (configure.get("ld_method", "") or "").strip()
        self.ld_method         = ld_method_raw or "indep"           # or indep-pairwise
        self.ld_window         = _coerce_int("ld_window", configure.get("ld_window",""), 50)
        self.ld_step           = _coerce_int("ld_step",   configure.get("ld_step",""),   5)
        try:
            self.ld_param = float(configure.get("ld_param", 2.0))
        except (TypeError, ValueError):
            self.ld_param = 2.0

        self.gff_biotype = configure.get("gff_biotype", "protein_coding")
        self.circos_window = _coerce_int("circos_window", configure.get("circos_window",""), 10000)  # in bp
        self.circos_tick_step   = _coerce_int("circos_tick_step",   configure.get("circos_tick_step",   ""), 0) or None
        chroms_raw = (configure.get("circos_chroms", "") or "").replace(";", ",")
        self.circos_chroms = [c.strip() for c in chroms_raw.split(",") if c.strip()]
        self.target_seq = _coerce_int("target_seq", configure.get("target_seq", ""), 0)
        # Default to Docker-mounted Kraken DB path; override via config if provided
        default_kraken_db = "/kraken_db"
        self.kraken_db = (configure.get("kraken_db", default_kraken_db) or default_kraken_db).strip()
        tag_outliers_raw = (configure.get("tag_outliers", "true") or "true").strip().lower()
        self.tag_outliers = tag_outliers_raw not in {"0", "false", "no", "off"}
        prediction_python_raw = (configure.get("prediction_python", "") or "").strip()
        if prediction_python_raw:
            self.prediction_python = prediction_python_raw
        else:
            ml_python = Path("/opt/conda/envs/DGbreeding-ml/bin/python")
            self.prediction_python = str(ml_python) if ml_python.exists() else sys.executable

    def _require_config_path(self) -> Path:
        if self.config_path is None:
            raise ValueError("config_file is required to run pipeline steps.")
        config_path = Path(self.config_path)
        if not config_path.exists():
            raise FileNotFoundError(f"config_file not found: {config_path}")
        if not config_path.is_file():
            raise FileNotFoundError(f"config_file is not a file: {config_path}")
        return config_path

    def _prepare_capture_bed(self) -> None:
        if self.capture_bed:
            cap_path = Path(self.capture_bed).expanduser().resolve()
            if cap_path.exists():
                self.capture_bed = cap_path
                return
            if not self.generate_bed:
                return

        if not self.generate_bed:
            return
        if not self.ref_fasta:
            raise ValueError("ref_fasta is required to generate capture BED.")
        if self.bed_require_genes and (not self.gff or not self.gff.exists()):
            raise ValueError("gff is required when bed_require_genes is true.")
        if self.bed_require_genes and (self.bed_top_n > 0 or self.bed_chroms):
            raise ValueError(
                "bed_require_genes=true cannot be combined with bed_top_n/bed_chroms. "
                "Choose one: gene-only (bed_require_genes=true, leave bed_top_n/bed_chroms blank) "
                "or top_n/chroms (bed_require_genes=false)."
            )

        if not self.bed_chroms and self.bed_top_n <= 0 and not self.bed_require_genes:
            time_stamp("Note: generate_bed enabled without filters; BED will include all contigs.")

        out_bed = self.bed_output or (self.variant_outdir / "glnexus.capture.bed")
        selected, missing, filtered = build_capture_bed(
            fasta=self.ref_fasta,
            out_bed=out_bed,
            top_n=self.bed_top_n,
            chroms=self.bed_chroms,
            require_genes=self.bed_require_genes,
            gff=self.gff if self.bed_require_genes else None,
            gene_biotype=self.bed_gene_biotype,
        )
        if not selected:
            time_stamp("Warning: Generated capture BED is empty; skipping --bed.")
            return
        if missing:
            preview = ", ".join(missing[:10])
            suffix = " ..." if len(missing) > 10 else ""
            time_stamp(f"Warning: {len(missing)} BED chroms not in FASTA: {preview}{suffix}")
        if filtered:
            preview = ", ".join(filtered[:10])
            suffix = " ..." if len(filtered) > 10 else ""
            time_stamp(f"{len(filtered)} BED chroms filtered (no genes in GFF): {preview}{suffix}")
        retained = len(selected)
        skipped_missing = len(missing)
        skipped_filtered = len(filtered)
        time_stamp(
            "Capture BED summary: retained {0} contigs; skipped {1} not in FASTA; skipped {2} no genes.".format(
                retained, skipped_missing, skipped_filtered
            )
        )
        time_stamp(f"Generated capture BED for GLnexus: {out_bed} ({retained} contigs)")
        self.capture_bed = out_bed

    def _cleanup_temp_cache(self) -> None:
        for name in (".tmp", ".cache"):
            path = self.outdir / name
            if not path.exists():
                continue
            if not path.is_dir():
                time_stamp(f"Warning: Expected directory for cleanup but found file: {path}")
                continue
            try:
                shutil.rmtree(path)
                time_stamp(f"Removed temp/cache directory: {path}")
            except Exception as exc:
                time_stamp(f"Warning: Failed to remove {path}: {exc}")

    def _run_core_step(self, step_id: str) -> None:
        if step_id == "1":
            time_stamp("Step 1: Preprocessing with fastp and index reference genome")
            RefIndex(self).run()
            fastp_preprocessing.run_fastp(
                wd=self.input_dir,
                outdir=self.fastp_outdir,
                threads=self.fastp_threads,
                min_length=self.fastp_minlength,
                qualified_phred=self.fastp_quality,
                average_qual=self.fastp_average_qual,
                lib_type="DNA",
                trim_front=self.fastp_trim_front,
                detect_adapter_pe=True,
                dedup=True,
                report_dir=self.fastp_report_dir,
            )
            time_stamp("Finished preprocessing with fastp")
            return
        if step_id == "2":
            time_stamp("Step 2: Aligning DNA reads to reference genome")
            DNAAlignment(self).run()
            time_stamp("Finished aligning DNA reads to reference genome")
            return
        if step_id == "3":
            time_stamp("Step 3: Running DeepVariant + GLnexus + plink")
            self._prepare_capture_bed()
            # Pass self directly; VariantCalling reads attributes it needs and uses defaults for the rest
            VariantCalling(self).run_all()
            time_stamp("Finished DeepVariant + GLnexus + plink")
            return
        if step_id == "4":
            if self._gwas_runner is None:
                self._gwas_runner = GWASPCA(self, self._require_config_path())
            time_stamp("Step 4: Running BLINK GWAS/PCA")
            self._gwas_runner.run()
            time_stamp("Finished BLINK GWAS/PCA")
            return
        if step_id == "5":
            if self._model_runner is None:
                self._model_runner = ModelConstruction(self, self._require_config_path())
            time_stamp("Step 5: Running SNP prediction model")
            self._model_runner.run()
            time_stamp("Finished SNP prediction model")
            return
        raise ValueError(f"Unknown core step: {step_id}")

    def _run_eval_step(self, step_id: str) -> None:
        if step_id not in self.eva_steps:
            return
        if step_id == "0":
            time_stamp("[eva0] FastQC + MultiQC")
            if self._ev_qc is None:
                self._ev_qc = FastQCEvaluator(self)
            self._ev_qc.evaluate_qc(fastq_stage="raw")
            return
        if step_id == "1":
            time_stamp("[eva1] FastQC + MultiQC + Kraken2 (all reads)")
            if self._ev_qc is None:
                self._ev_qc = FastQCEvaluator(self)
            if self._ev_contam is None:
                self._ev_contam = ContaminationEvaluator(self)
            self._ev_qc.evaluate_qc(fastq_stage="processed")
            self._ev_contam.evaluate_contamination()
            return
        if step_id == "2":
            time_stamp("[eva2] Mapping yield (unique/fastp)")
            _ev_align = MappingYieldEvaluator(self)
            _ev_align.evaluate_alignment()
            return
        if step_id == "3":
            time_stamp("[eva3] Variant stats + Circos")
            if self._ev_var is None:
                self._ev_var = VariantCircosEvaluator(self)
            self._ev_var.evaluate_variants()
            return
        if step_id == "4":
            if self._gwas_runner is None:
                self._gwas_runner = GWASPCA(self, self._require_config_path())
            time_stamp("[eva4] Manhattan/QQ plots")
            self._gwas_runner.run_eval()
            return
        if step_id == "5":
            if self._model_runner is None:
                self._model_runner = ModelConstruction(self, self._require_config_path())
            time_stamp("[eva5] Prediction ROC/AUC plots")
            self._model_runner.run_eval()
            return
        raise ValueError(f"Unknown evaluation step: {step_id}")

    def run(self):
        time_stamp("Starting Digital Breeding Prediction Pipeline")
        try:
            steps_to_run = self.steps - {"0"}
            if steps_to_run or self.eva_steps:
                self._require_config_path()
            if not steps_to_run and not self.eva_steps:
                print("No steps selected.")
                return
            if not steps_to_run and self.eva_steps:
                print("Evaluation-only mode.")

            # Evaluation step 0 can run before any analysis
            self._run_eval_step("0")

            for step_id in ("1", "2", "3", "4", "5"):
                if step_id in steps_to_run:
                    self._run_core_step(step_id)
                self._run_eval_step(step_id)
            time_stamp("Pipeline complete")
        finally:
            self._cleanup_temp_cache()


def main():
    parser = argparse.ArgumentParser(description="Digital Breeding Prediction Pipeline")
    parser.add_argument('--config_file', type=str, required=True, help='Path to the configuration file.')
    args = parser.parse_args()
    config_path = Path(args.config_file).resolve()
    configure = read_config(config_path)
    pipeline = PredictionPipeline(configure, config_path=config_path)
    pipeline.run()

if __name__ == "__main__":
    main()
