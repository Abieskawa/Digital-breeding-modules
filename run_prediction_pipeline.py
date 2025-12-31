#!/usr/bin/env python3
import argparse, os, shutil, subprocess, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set
import pandas as pd

import Preprocess.fastp_preprocessing as fastp_preprocessing
from DNA.dna_alignment import DNAAlignment
from DNA.ref_index import RefIndex
from DNA.variant_calling import VariantCalling
from Evaluation.contamination import ContaminationEvaluator
from Evaluation.evaluate_mapping import MappingYieldEvaluator
from Evaluation.fastq_qc import FastQCEvaluator
from Evaluation.variant_circos import VariantCircosEvaluator
from Utils.bed_generator import build_capture_bed
from Utils.utils import _resolve_outdir, _ensure_bgzip_and_index, read_config, time_stamp, parse_int_list, coerce_bool
from Utils.cv_preparation import CV_preparation
from Prediction_model.gwas_pca import GWAS_PCA, PredictionGWASCVStep
from Prediction_model.model_construction import ModelConstruction, PredictionModelCVStep
from Prediction_model.snp_numeric_transformer import SNP_numerical, _save_numeric_npz

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


def _parse_list(raw: str) -> Set[str]:
    """
    Parse comma/semicolon-separated list into a set of non-empty strings.
    """
    return set(_parse_list_ordered(raw))


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


def _load_fold_job(job: Dict[str, object]):
    config_file = str(job["config_file"])
    repo_root = Path(str(job["repo_root"]))
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    config = read_config(config_file)
    return (
        config,
        str(job["file_fold_index"]),
        str(job.get("gwas_pcs", "") or ""),
        str(job.get("top_n_snps", "") or ""),
        str(job.get("models", "") or ""),
        int(job.get("inner_model_workers", 1)),
    )


def _cv_gwas_worker(job: Dict[str, object]) -> str:
    fold_cfg = CV_preparation.read_fold_config(Path(str(job["fold_config"])))
    config = read_config(fold_cfg["config_file"])
    file_fold_index = fold_cfg["file_fold_index"]
    gwas_pcs = str(fold_cfg.get("gwas_pcs", "") or "")
    pca_n_components = CV_preparation.parse_optional_int(fold_cfg.get("pca_n_components"))

    train_samples = CV_preparation.read_samples_file(Path(fold_cfg["train_samples_file"]))
    if not train_samples:
        raise ValueError(f"No train samples found in {fold_cfg['train_samples_file']}")

    snp_tool = SNP_numerical(config)
    gwas_tool = GWAS_PCA(config)

    prediction_outdir = _resolve_outdir(
        config,
        key="prediction_output_dir",
        default="results",
        resolve=True,
        ensure_dir=True,
    )
    tmp_dir = prediction_outdir / ".tmp" / "prediction_cv" / f"fold_{file_fold_index}"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    npz_path = tmp_dir / "train_numeric.npz"

    genotype_df = snp_tool.vcf_to_dataframe(snp_tool.vcf_file_path)
    present = [s for s in train_samples if s in genotype_df.index]
    if not present:
        raise ValueError("No requested train samples were found in the VCF.")
    genotype_numeric = snp_tool.convert_genotypes_to_numeric(genotype_df.loc[present])
    _save_numeric_npz(genotype_numeric, str(npz_path))

    gwas_tool.run_fold(
        file_fold_index,
        str(npz_path),
        gwas_pcs_override=gwas_pcs,
        pca_components_override=pca_n_components,
    )
    return file_fold_index


def _cv_model_prep_worker(job: Dict[str, object]) -> str:
    fold_cfg = CV_preparation.read_fold_config(Path(str(job["fold_config"])))
    config = read_config(fold_cfg["config_file"])
    file_fold_index = fold_cfg["file_fold_index"]
    gwas_pcs = str(fold_cfg.get("gwas_pcs", "") or "")
    top_n_snps = str(fold_cfg.get("top_n_snps", "") or "")
    models = str(fold_cfg.get("models", "") or "")
    inner_model_workers = int(job.get("inner_model_workers", 1))

    test_samples = CV_preparation.read_samples_file(Path(fold_cfg["test_samples_file"]))
    if not test_samples:
        raise ValueError(f"No test samples found in {fold_cfg['test_samples_file']}")

    snp_tool = SNP_numerical(config)
    model_tool = ModelConstruction(config)

    genotype_df = snp_tool.vcf_to_dataframe(snp_tool.vcf_file_path)
    present = [s for s in test_samples if s in genotype_df.index]
    if not present:
        raise ValueError("No requested test samples were found in the VCF.")
    test_df = genotype_df.loc[present]

    pcs = parse_int_list(gwas_pcs) if gwas_pcs.strip() else snp_tool.gwas_n_pcs_list
    nsnps = parse_int_list(top_n_snps) if top_n_snps.strip() else snp_tool.top_n_snps_list

    for pc in pcs:
        file_fold_npc_index = f"{file_fold_index}_{pc}PCs"
        for n_snp in nsnps:
            file_fold_npc_nsnp_index = f"{file_fold_npc_index}_{n_snp}SNPs"
            selected = snp_tool.select_markers_and_prepare_data_for_model_construction(
                file_fold_npc_index=file_fold_npc_index,
                file_fold_npc_nsnp_index=file_fold_npc_nsnp_index,
                test_df=test_df,
                top_n_snps=n_snp,
            )
            save_dir = os.path.join(snp_tool.model_output_dir, file_fold_npc_nsnp_index)
            with open(os.path.join(save_dir, "selected_snps.txt"), "w") as f:
                for snp in selected:
                    f.write(f"{snp}\n")

    models_list = [m.strip() for m in models.split(",") if m.strip()] if models.strip() else model_tool.models
    model_tool.run_fold(
        file_fold_index,
        pcs=pcs,
        nsnps=nsnps,
        models=models_list,
        max_workers=inner_model_workers,
    )
    return file_fold_index


def _cv_fold_model_worker(job: Dict[str, object]) -> str:
    """
    Worker for ONE fold job. Runs model training only (no GWAS/PCA/model_prep).
    """

    (
        config,
        file_fold_index,
        gwas_pcs,
        top_n_snps,
        models,
        inner_model_workers,
    ) = _load_fold_job(job)
    model_tool = ModelConstruction(config)

    pcs = parse_int_list(gwas_pcs) if gwas_pcs.strip() else model_tool.gwas_n_pcs_list
    nsnps = parse_int_list(top_n_snps) if top_n_snps.strip() else model_tool.top_n_snps_list
    models_list = [m.strip() for m in models.split(",") if m.strip()] if models.strip() else model_tool.models

    model_tool.run_fold(
        file_fold_index,
        pcs=pcs,
        nsnps=nsnps,
        models=models_list,
        max_workers=inner_model_workers,
    )
    return file_fold_index


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

        # Paths
        self.outdir = _resolve_outdir(
            configure,
            key="output_dir",
            default="DGBreeding",
            resolve=True,
            ensure_dir=True,
        )
        self.prediction_outdir = _resolve_outdir(
            configure,
            key="prediction_output_dir",
            base_outdir=self.outdir,
            default="results",
            resolve=True,
            ensure_dir=True,
        )
        self.input_dir = Path(configure.get("input_dir", ".")).resolve()
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
        self.mem = configure.get("mem", "8G")

        # fastp
        self.fastp_threads = _coerce_int("fastp_threads", configure.get("fastp_threads", ""), self.threads)
        self.fastp_quality = _coerce_int("fastp_quality", configure.get("fastp_quality", ""), 20)
        self.fastp_minlength = _coerce_int("fastp_length", configure.get("fastp_length", ""), 30)
        self.fastp_average_qual = _coerce_int("fastp_average_qual", configure.get("fastp_average_qual", "20"), 20)
        self.fastp_trim_front = configure.get("fastp_trim_front", "10")
        self.fastp_outdir = _resolve_outdir(
            base_outdir=self.outdir,
            subdir="00_Preprocessed_DNA",
            resolve=True,
            ensure_dir=True,
        )
        self.fastp_report_dir = _resolve_outdir(
            base_outdir=self.fastp_outdir,
            subdir="fastp_reports",
            resolve=True,
            ensure_dir=True,
        )

        # Alignment context
        self.population = configure.get("population", "diversity_panel")
        self.seq_platform = (configure.get("seq_platform", "HiSeq")).strip().lower()
        self.optical_distance = _coerce_int("optical_distance", configure.get("optical_distance", ""), 100)
        self.mark_duplicates = coerce_bool(configure.get("mark_duplicates", "true"), True)
        self.filter_proper_pairs = coerce_bool(configure.get("filter_proper_pairs", "true"), True)
        self.align_outdir = _resolve_outdir(
            base_outdir=self.outdir,
            subdir="01_DNAseq_alignment",
            resolve=True,
            ensure_dir=True,
        )
        # Optional: skip specific samples at alignment/variant stages
        self.skip_align_samples = _parse_list(configure.get("skip_alignment_samples", ""))
        # Default variant skip list to alignment skip list if not provided separately
        raw_variant_skip = configure.get("skip_variant_samples", "")
        self.skip_variant_samples = _parse_list(raw_variant_skip) or set(self.skip_align_samples)

        # Variant calling context (directory only; tool options use VariantCalling defaults
        # unless you add them to configure.txt)
        self.variant_outdir = _resolve_outdir(base_outdir=self.outdir, subdir="02_Variant_Calling", resolve=True)

        # Optional VC knobs (leave unset if not present; VariantCalling has sane defaults)
        self.deepvariant_image = configure.get("deepvariant_image", "google/deepvariant:1.9.0")
        self.deepvariant_use_gpu = coerce_bool(configure.get("deepvariant_use_gpu", ""), False)
        self.model_type = configure.get("model_type", "WGS")  # or WES
        self.capture_bed = configure.get("capture_bed", "") or None
        self.generate_bed = coerce_bool(configure.get("generate_bed", ""), False)
        self.bed_top_n = _coerce_int("bed_top_n", configure.get("bed_top_n", ""), 0)
        self.bed_chroms = _parse_list_ordered(configure.get("bed_chroms", ""))
        self.bed_require_genes = coerce_bool(configure.get("bed_require_genes", ""), False)
        self.bed_gene_biotype = (configure.get("bed_gene_biotype", "") or "").strip() or None
        bed_output_raw = (configure.get("bed_output", "") or "").strip()
        self.bed_output = Path(bed_output_raw).expanduser().resolve() if bed_output_raw else None
        self.hwe_p = configure.get("hwe_p", "1e-5")
        self.geno_missing = configure.get("geno_missing", "0.10")
        ld_prune_raw = configure.get("ld_prune", "")
        if ld_prune_raw == "":
            ld_prune_raw = configure.get("ld_pruning", "")
        self.enable_ld_pruning = coerce_bool(ld_prune_raw, True)
        ld_method_raw = (configure.get("ld_method", "") or "").strip()
        self.ld_method = ld_method_raw or "indep"  # or indep-pairwise
        self.ld_window = _coerce_int("ld_window", configure.get("ld_window", ""), 50)
        self.ld_step = _coerce_int("ld_step", configure.get("ld_step", ""), 5)
        try:
            self.ld_param = float(configure.get("ld_param", 2.0))
        except (TypeError, ValueError):
            self.ld_param = 2.0

        self.gff_biotype = configure.get("gff_biotype", "protein_coding")
        self.circos_window = _coerce_int("circos_window", configure.get("circos_window", ""), 10000)  # in bp
        self.circos_tick_step = _coerce_int("circos_tick_step", configure.get("circos_tick_step", ""), 0) or None
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

        # --- Prediction / CV knobs (used in step 4/5 CV orchestration) ---
        self.pred_vcf_file = Path((configure.get("vcf_file_path") or "").strip()).expanduser().resolve() if (configure.get("vcf_file_path") or "").strip() else None
        self.pred_pheno_csv = Path((configure.get("phenotype_csv_path") or "").strip()).expanduser().resolve() if (configure.get("phenotype_csv_path") or "").strip() else None
        self.pred_pheno_col = (configure.get("phenotype_column", "Phenotype") or "Phenotype").strip()
        self.pred_data_type = (configure.get("data_type", "binary") or "binary").strip().lower()

        self.pred_num_splits = _coerce_int("num_splits", configure.get("num_splits", ""), 5)
        self.pred_random_state = _coerce_int("random_state", configure.get("random_state", ""), 2024)
        folds_path_raw = (configure.get("cv_folds_file") or "").strip()
        if folds_path_raw:
            self.pred_cv_folds_file = Path(folds_path_raw).expanduser().resolve()
        else:
            self.pred_cv_folds_file = self.prediction_outdir / "cv_folds.tsv"
        # PCs used as GWAS covariates (separate from PCA extraction count)
        self.pred_gwas_pcs = (configure.get("gwas_n_pcs_list") or configure.get("gwas_n_pcs") or configure.get("n_pcs_list") or "0").strip()
        self.pred_pca_n_components = (configure.get("pca_n_components", "") or "").strip()

        self.pred_top_n_snps = (configure.get("top_n_snps_list", "10") or "10").strip()
        self.pred_models = (configure.get("models", "SVM,RandomForest,XGBoost") or "SVM,RandomForest,XGBoost").strip()
        self.pred_num_threads = _coerce_int("num_threads", configure.get("num_threads", ""), max(1, min(8, self.threads)))
        self.pred_cv_workers = _coerce_int("cv_workers", configure.get("cv_workers", ""), max(1, min(self.pred_num_splits, self.pred_num_threads)))

        # Inside each fold, keep model training sequential by default (your principle: per split parallel)
        self.pred_inner_model_workers = _coerce_int("inner_model_workers", configure.get("inner_model_workers", ""), 1)
        step5_mode_raw = (configure.get("step5_mode") or "full").strip().lower()
        if step5_mode_raw not in {"full", "model_only"}:
            time_stamp(f"Warning: Unknown step5_mode '{step5_mode_raw}', defaulting to 'full'.")
            step5_mode_raw = "full"
        self.pred_step5_mode = step5_mode_raw

    def _require_config_path(self) -> Path:
        if self.config_path is None:
            raise ValueError("config_file is required to run pipeline steps.")
        config_path = Path(self.config_path)
        if not config_path.exists():
            raise FileNotFoundError(f"config_file not found: {config_path}")
        if not config_path.is_file():
            raise FileNotFoundError(f"config_file is not a file: {config_path}")
        return config_path

    def _require_prediction_inputs(self) -> tuple:
        if self.pred_vcf_file is None or not self.pred_vcf_file.exists():
            raise FileNotFoundError("vcf_file_path is required and must exist.")
        if self.pred_pheno_csv is None or not self.pred_pheno_csv.exists():
            raise FileNotFoundError("phenotype_csv_path is required and must exist.")
        return self.pred_vcf_file, self.pred_pheno_csv

    def _prepare_prediction_cv_configs(
        self,
        *,
        create_if_missing: bool,
        ensure_bgzip: bool = False,
    ) -> List[Path]:
        cfg_path = self._require_config_path()
        vcf_path, pheno_path = self._require_prediction_inputs()
        if ensure_bgzip:
            _ensure_bgzip_and_index(vcf_path)
        cv_config_dir = self.prediction_outdir / "cv_configs"
        return CV_preparation.prepare_cv_fold_configs(
            cfg_path=cfg_path,
            vcf_path=vcf_path,
            pheno_path=pheno_path,
            pheno_col=self.pred_pheno_col,
            data_type=self.pred_data_type,
            num_splits=self.pred_num_splits,
            random_state=self.pred_random_state,
            cv_folds_file=self.pred_cv_folds_file,
            cv_config_dir=cv_config_dir,
            gwas_pcs=self.pred_gwas_pcs,
            pca_n_components=self.pred_pca_n_components if self.pred_pca_n_components else None,
            top_n_snps=self.pred_top_n_snps,
            models=self.pred_models,
            create_if_missing=create_if_missing,
        )

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
        for base in {self.outdir, self.prediction_outdir}:
            for name in (".tmp", ".cache"):
                path = base / name
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

    def _run_prediction_eval_script(self, rel_script: str) -> None:
        
        repo_root = Path(__file__).resolve().parent
        cmd = [self.prediction_python, str(repo_root / rel_script), str(self._require_config_path())]
        subprocess.run(cmd, check=True, cwd=repo_root)

    def _run_prediction_gwas_step4(self) -> None:
        """
        Step 4: Run GWAS/PCA per CV fold (train split only).
        """
        
        fold_configs = self._prepare_prediction_cv_configs(
            create_if_missing=True,
            ensure_bgzip=True,
        )

        time_stamp(f"[step4] CV fold configs ready: {len(fold_configs)} folds")
        runner = PredictionGWASCVStep(fold_configs, self.pred_cv_workers, _cv_gwas_worker)
        finished = runner.run()
        time_stamp(f"[step4] Completed folds: {len(finished)}")

    def _run_prediction_cv_step5(self) -> None:
        """
        Step 5: CV orchestration (fold-level parallel) for model prep + training.
        """
        
        if self.pred_step5_mode == "model_only":
            cfg_path = self._require_config_path()
            model_dir = self.prediction_outdir / "Model_result"
            if not model_dir.exists():
                raise FileNotFoundError(f"Model_result not found: {model_dir}")

            fold_map: Dict[str, Dict[str, Set[int]]] = {}
            for entry in model_dir.iterdir():
                if not entry.is_dir():
                    continue
                parts = entry.name.split("_")
                if len(parts) < 3:
                    continue
                fold = "_".join(parts[:-2])
                pc_part = parts[-2]
                snp_part = parts[-1]
                if not (pc_part.endswith("PCs") and snp_part.endswith("SNPs")):
                    continue
                try:
                    pc = int(pc_part[:-3])
                    n_snp = int(snp_part[:-4])
                except ValueError:
                    continue
                if fold not in fold_map:
                    fold_map[fold] = {"pcs": set(), "snps": set()}
                fold_map[fold]["pcs"].add(pc)
                fold_map[fold]["snps"].add(n_snp)

            if not fold_map:
                raise ValueError(f"No model input directories found under {model_dir}")

            jobs = []
            for fold in sorted(fold_map):
                pcs = sorted(fold_map[fold]["pcs"])
                nsnps = sorted(fold_map[fold]["snps"])
                jobs.append(
                    {
                        "config_file": str(cfg_path),
                        "repo_root": str(Path(__file__).resolve().parent),
                        "file_fold_index": fold,
                        "gwas_pcs": ",".join(str(p) for p in pcs),
                        "top_n_snps": ",".join(str(n) for n in nsnps),
                        "models": self.pred_models,
                        "inner_model_workers": int(self.pred_inner_model_workers),
                    }
                )

            time_stamp(f"[step5] mode=model_only: {len(jobs)} folds discovered under {model_dir}")

            finished: List[str] = []
            with ProcessPoolExecutor(max_workers=int(self.pred_cv_workers)) as ex:
                futs = [ex.submit(_cv_fold_model_worker, j) for j in jobs]
                for fut in as_completed(futs):
                    finished.append(fut.result())

            time_stamp(f"[step5] Completed folds: {len(finished)}")
            return

        fold_configs = self._prepare_prediction_cv_configs(create_if_missing=False)

        time_stamp(f"[step5] CV fold configs ready: {len(fold_configs)} folds")
        runner = PredictionModelCVStep(
            fold_configs,
            self.pred_cv_workers,
            self.pred_inner_model_workers,
            _cv_model_prep_worker,
        )
        finished = runner.run()
        time_stamp(f"[step5] Completed folds: {len(finished)}")

        # Summarize metrics (replicates cat_result behavior; avoids re-running training)
        model_dir = self.prediction_outdir / "Model_result"
        if model_dir.exists():
            result_files = [p for p in model_dir.iterdir() if p.is_file() and p.name.startswith("model_evaluation_results_") and p.suffix == ".csv"]
            if result_files:
                dfs = [pd.read_csv(p) for p in result_files]
                combined = pd.concat(dfs, ignore_index=True)
                combined = combined.sort_values(by=["n_pcs", "n_snps", "model_name"])
                mean_df = combined.groupby(["n_pcs", "n_snps", "model_name"])[["accuracy", "specificity", "mcc", "f1_score"]].mean().reset_index()
                combined.to_csv(self.prediction_outdir / "model_evaluation_results_all.csv", index=False)
                mean_df.to_csv(self.prediction_outdir / "model_evaluation_results_mean.csv", index=False)
                time_stamp("[step5] Wrote model_evaluation_results_all.csv and model_evaluation_results_mean.csv")

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
            time_stamp("Step 4: Running BLINK GWAS/PCA (per CV fold)")
            self._run_prediction_gwas_step4()
            time_stamp("Finished BLINK GWAS/PCA")
            return
        if step_id == "5":
            time_stamp("Step 5: Running SNP prediction model (CV orchestration in pipeline)")
            self._run_prediction_cv_step5()
            time_stamp("Finished SNP prediction model (CV orchestration)")
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
            time_stamp("[eva4] Manhattan/QQ plots")
            self._run_prediction_eval_script("Evaluation/gwas_manhattan_plot.py")
            return
        if step_id == "5":
            time_stamp("[eva5] Prediction ROC/AUC plots")
            self._run_prediction_eval_script("Evaluation/prediction_roc_plot.py")
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
