#!/usr/bin/env python3
from typing import List, Dict, Set
from Utils.utils import time_stamp
from pathlib import Path
import Preprocess.fastp_preprocessing as fastp_preprocessing
from DNA.ref_index import RefIndex
from DNA.dna_alignment import DNAAlignment
from DNA.variant_calling import VariantCalling
import argparse
from Evaluation.evaluate_mapping import MappingEvaluator
from Evaluation.variant_circos import VariantCircosEvaluator

# ----helpers-----
def load_config_file(config_file: Path) -> Dict[str, str]:
    """
    Simple key=value loader with support for blank lines, full-line comments,
    and inline comments after '#'.
    """
    cfg: Dict[str, str] = {}
    with open(config_file, "r") as cf:
        for raw in cf:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "#" in line:
                line = line.split("#", 1)[0].strip()
                if not line:
                    continue
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            cfg[k.strip()] = v.strip()
    return cfg


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


class PredictionPipeline(object):
    def __init__(self, configure: Dict[str,str]):
        self.step_raw = configure.get("step", "0")
        self.steps = _parse_steps(self.step_raw)
        self.eva_step_raw = configure.get("eva_step", "0")
        self.eva_steps = _parse_steps(self.eva_step_raw)

        # Paths
        self.outdir = Path(configure.get("output_dir","DGBreeding")).resolve()
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
        self.fastp_quality = _coerce_int("fastp_quality", configure.get("fastp_quality", ""), 20)
        self.fastp_minlength = _coerce_int("fastp_length", configure.get("fastp_length", ""), 15)
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
        self.hwe_p             = configure.get("hwe_p", "1e-5")
        self.geno_missing      = configure.get("geno_missing", "0.10")
        self.ld_method         = configure.get("ld_method", "indep")           # or indep-pairwise
        self.ld_window         = _coerce_int("ld_window", configure.get("ld_window",""), 50)
        self.ld_step           = _coerce_int("ld_step",   configure.get("ld_step",""),   5)
        try:
            self.ld_param = float(configure.get("ld_param", 2.0))
        except (TypeError, ValueError):
            self.ld_param = 2.0

        self.gff_biotype = configure.get("gff_biotype", "protein_coding")
        self.circos_window = _coerce_int("circos_window", configure.get("circos_window",""), 10000)  # in bp
        self.circos_gene_window = _coerce_int("circos_gene_window", configure.get("circos_gene_window", ""), 0) or None
        self.circos_var_window  = _coerce_int("circos_var_window",  configure.get("circos_var_window",  ""), 0) or None
        self.circos_tick_step   = _coerce_int("circos_tick_step",   configure.get("circos_tick_step",   ""), 0) or None
        chroms_raw = (configure.get("circos_chroms", "") or "").replace(";", ",")
        self.circos_chroms = [c.strip() for c in chroms_raw.split(",") if c.strip()]
        self.target_seq = _coerce_int("target_seq", configure.get("target_seq", ""), 0)
        # Default to Docker-mounted Kraken DB path; override via config if provided
        default_kraken_db = "/kraken_db"
        self.kraken_db = (configure.get("kraken_db", default_kraken_db) or default_kraken_db).strip()
        tag_outliers_raw = (configure.get("tag_outliers", "true") or "true").strip().lower()
        self.tag_outliers = tag_outliers_raw not in {"0", "false", "no", "off"}

    def run_step1_fastp(self):
        time_stamp("Step 1: Preprocessing with fastp and index reference genome")
        ri = RefIndex(self); ri.run()
        fastp_preprocessing.run_fastp(
            wd = self.input_dir,
            outdir = self.fastp_outdir,
            threads = self.threads,
            min_length = self.fastp_minlength,
            qualified_phred = self.fastp_quality,
            average_qual = self.fastp_average_qual,
            lib_type = "DNA",
            trim_front = self.fastp_trim_front,
            detect_adapter_pe = True,
            dedup = True,
            report_dir = self.fastp_report_dir
        )
        time_stamp("Finished preprocessing with fastp")

    def run_step2_alignement(self):
        time_stamp("Step 2: Aligning DNA reads to reference genome")
        dna_align = DNAAlignment(self); dna_align.run()
        time_stamp("Finished aligning DNA reads to reference genome")

    def run_step3_deepvariant(self):
        time_stamp("Step 3: Running DeepVariant + GLnexus + plink")
        # Pass self directly; VariantCalling reads attributes it needs and uses defaults for the rest
        vc = VariantCalling(self)
        vc.run_all()
        time_stamp("Finished DeepVariant + GLnexus + plink")

    def run(self):
        time_stamp("Starting Digital Breeding Prediction Pipeline")
        steps_to_run = self.steps - {"0"}
        if not steps_to_run and not self.eva_steps:
            print("No steps selected.")
            return
        if not steps_to_run and self.eva_steps:
            print("Evaluation-only mode.")

        ev_map = MappingEvaluator(self)
        ev_var = VariantCircosEvaluator(self)

        # Evaluation step 0 can run before any analysis
        if "0" in self.eva_steps:
            time_stamp("[eva0] FastQC + MultiQC")
            ev_map.evaluate_qc(fastq_stage="raw")

        if "1" in steps_to_run:
            self.run_step1_fastp()
        if "1" in self.eva_steps:
            time_stamp("[eva1] FastQC + MultiQC + Kraken2 (all reads)")
            ev_map.evaluate_qc(fastq_stage="processed")
            ev_map.evaluate_contamination()

        if "2" in steps_to_run:
            self.run_step2_alignement()
        if "2" in self.eva_steps:
            time_stamp("[eva2] Mapping yield (unique/fastp)")
            ev_map.evaluate_alignment()

        if "3" in steps_to_run:
            self.run_step3_deepvariant()
        if "3" in self.eva_steps:
            time_stamp("[eva3] Variant stats + Circos")
            ev_var.evaluate_variants()

        time_stamp("Pipeline complete")


def main():
    parser = argparse.ArgumentParser(description="Digital Breeding Prediction Pipeline")
    parser.add_argument('--config_file', type=str, required=True, help='Path to the configuration file.')
    args = parser.parse_args()
    configure = load_config_file(Path(args.config_file))
    pipeline = PredictionPipeline(configure)
    pipeline.run()

if __name__ == "__main__":
    main()
