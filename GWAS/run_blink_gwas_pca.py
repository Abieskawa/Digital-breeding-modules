#!/usr/bin/env python3
"""
Step 4: PCA + GWAS (BLINK)

Inputs:
  - VCF (config: vcf_file_path)
  - Numeric genotype NPZ (from Step 1 pca_numeric)

Outputs:
  output_dir/GWAS_dir/{file_fold_index}_{pc}PCs/
    train_{file_fold_index}_{pc}PCs.txt      (phenotype)
    train_{file_fold_index}_{pc}PCs.vcf      (subset VCF for BLINK)
    train_{file_fold_index}_{pc}PCs.cov      (covariates, if pc>0 and use_pcs)
    blink_run_train_{...}.log
    train_{...}_{phenotype_column}_GWAS_result.txt (BLINK output)
"""
import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from Utils.utils import (
    _resolve_outdir,
    _ensure_bgzip_and_index,
    parse_int_list,
    extract_vcf_samples,
    encode_phenotype_series,
)


def _load_numeric_npz(npz_path: str) -> pd.DataFrame:
    data = np.load(npz_path, allow_pickle=True)
    genotype = data["genotype"]
    samples = data["samples"].astype(str)
    snps = data["snps"].astype(str)
    return pd.DataFrame(genotype, index=samples, columns=snps)


class GWAS_PCA:
    def __init__(self, config: dict):
        self.config = config
        self.species_name = self.config.get("species_name")
        self.phenotype_csv_path = self.config.get("phenotype_csv_path")
        self.vcf_file_path = self.config.get("vcf_file_path")
        self.chromosome_csv_path = self.config.get("chromosome_csv_path")
        self.phenotype_column = self.config.get("phenotype_column", "Phenotype")
        self.data_type = self.config.get("data_type", "binary")

        self.use_pcs = self.config.get("use_pcs", "True").lower() == "true"

        self.gwas_n_pcs_list = (
            parse_int_list(self.config.get("gwas_n_pcs_list", ""))
            or parse_int_list(self.config.get("gwas_n_pcs", ""))
            or parse_int_list(self.config.get("n_pcs_list", ""))
        )

        raw_pca_n = self.config.get("pca_n_components", "").strip()
        if raw_pca_n:
            try:
                self.pca_n_components = int(raw_pca_n)
            except Exception:
                max_pcs = max(self.gwas_n_pcs_list) if self.gwas_n_pcs_list else 0
                self.pca_n_components = max(15, max_pcs)
        else:
            max_pcs = max(self.gwas_n_pcs_list) if self.gwas_n_pcs_list else 0
            self.pca_n_components = max(15, max_pcs)

        self.output_dir = _resolve_outdir(
            self.config,
            key="prediction_output_dir",
            default="results",
            resolve=True,
            ensure_dir=True,
        )
        pca_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="evaluation/PCA_Scree_Plot",
            resolve=True,
            ensure_dir=True,
        )
        self.pca_scree_output_dir = pca_output_dir
        self.pca_eval_output_dir = pca_output_dir
        self.gwas_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="GWAS_dir",
            resolve=True,
            ensure_dir=True,
        )

        self.blink_bin = self._resolve_blink_bin(
            self.config.get("blink_bin", "/home/abieskawa/tools/blink/blink1_5")
        )

        self.load_chromosome_recode()

    @staticmethod
    def _ensure_executable(path: Path) -> Path:
        if os.access(path, os.X_OK):
            return path
        try:
            path.chmod(path.stat().st_mode | 0o111)
        except OSError as exc:
            raise PermissionError(
                f"BLINK binary exists but is not executable: {path}. chmod +x failed: {exc}"
            ) from exc
        if not os.access(path, os.X_OK):
            raise PermissionError(f"BLINK binary exists but is not executable: {path}")
        return path

    @classmethod
    def _resolve_blink_bin(cls, raw_blink_bin: str) -> str:
        blink_bin = (raw_blink_bin or "").strip()
        if not blink_bin:
            blink_bin = "/home/abieskawa/tools/blink/blink1_5"

        path_candidate = Path(blink_bin).expanduser()
        if path_candidate.exists():
            return str(cls._ensure_executable(path_candidate))

        resolved = shutil.which(blink_bin)
        if resolved:
            return resolved

        raise FileNotFoundError(
            f"BLINK binary not found: {blink_bin}. Set blink_bin=... in the config or install BLINK on PATH."
        )

    def load_chromosome_recode(self):
        self.chromosome_mapping = pd.read_csv(self.chromosome_csv_path)
        self.chr_name_to_recoded = dict(
            zip(
                self.chromosome_mapping["OriginalName"],
                self.chromosome_mapping["RecodedNumber"],
            )
        )
        self.recoded_to_label = dict(
            zip(
                self.chromosome_mapping["RecodedNumber"].astype(str),
                self.chromosome_mapping["Label"],
            )
        )

    def perform_pca(self, genotype_df: pd.DataFrame, n_components: int) -> Tuple[pd.DataFrame, np.ndarray]:
        max_components = min(genotype_df.shape[0], genotype_df.shape[1])
        if max_components < 1:
            raise ValueError("Not enough data to perform PCA.")
        n_components = min(int(n_components), max_components)

        scaler = StandardScaler()
        genotype_scaled = scaler.fit_transform(genotype_df)

        pca = PCA(n_components=n_components)
        principal_components = pca.fit_transform(genotype_scaled)
        explained_variance_ratio = pca.explained_variance_ratio_

        pc_columns = [f"PC{i + 1}" for i in range(n_components)]
        principal_df = pd.DataFrame(data=principal_components, columns=pc_columns)
        principal_df.insert(0, "taxa", genotype_df.index)
        return principal_df, explained_variance_ratio

    def _load_phenotypes(self) -> pd.DataFrame:
        phenotype_df = pd.read_csv(self.phenotype_csv_path, sep=None, engine="python")
        if "taxa" not in phenotype_df.columns:
            raise ValueError("phenotype_csv_path must contain a 'taxa' column.")
        if self.phenotype_column not in phenotype_df.columns:
            raise ValueError(f"phenotype column not found: {self.phenotype_column}")
        phenotype_df = phenotype_df[["taxa", self.phenotype_column]].dropna()
        phenotype_df["taxa"] = phenotype_df["taxa"].astype(str)
        return phenotype_df.set_index("taxa")

    def _resolve_pca_components(
        self,
        gwas_pcs_list: List[int],
        pca_components_override: Optional[int] = None,
    ) -> int:
        max_pcs = max(gwas_pcs_list) if gwas_pcs_list else 0
        n_components = max(self.pca_n_components, 15, max_pcs)
        if pca_components_override is not None:
            n_components = max(n_components, int(pca_components_override))
        return n_components

    def _write_pca_outputs(
        self,
        tag: str,
        principal_df: pd.DataFrame,
        explained_variance_ratio: np.ndarray,
    ) -> None:
        try:
            out_csv = os.path.join(self.pca_eval_output_dir, f"explained_variance_{tag}.csv")
            pd.DataFrame({"explained_variance_ratio": explained_variance_ratio}).to_csv(
                out_csv,
                index=False,
            )
            comp_csv = os.path.join(self.pca_eval_output_dir, f"pca_components_{tag}.csv")
            principal_df.to_csv(comp_csv, index=False)
        except Exception:
            pass

    def run_global_qc(self, *, pca_components_override: Optional[int] = None) -> None:
        from Prediction_model.snp_numeric_transformer import SNP_numerical
        from Evaluation.gwas_pca_eval import PCAEvaluator

        phenotype_df = self._load_phenotypes()
        genotype_tool = SNP_numerical(self.config)
        genotype_df = genotype_tool.vcf_to_dataframe(self.vcf_file_path)

        present = [s for s in phenotype_df.index if s in genotype_df.index]
        if not present:
            raise ValueError("No overlapping samples between VCF and phenotype CSV.")

        genotype_numeric = genotype_tool.convert_genotypes_to_numeric(genotype_df.loc[present])

        n_components = self._resolve_pca_components(
            self.gwas_n_pcs_list,
            pca_components_override=pca_components_override,
        )

        principal_df, explained_variance_ratio = self.perform_pca(
            genotype_numeric,
            n_components=n_components,
        )
        self._write_pca_outputs("global_qc", principal_df, explained_variance_ratio)

        pca_eval = PCAEvaluator(self.config)
        pca_eval.plot_global_qc(principal_df, explained_variance_ratio)

    def blink_file_generator(
        self,
        train_samples: List[str],
        file_fold_npc_index: str,
        principal_df: pd.DataFrame,
        n_pcs: int,
    ):
        """
        Generate input files for BLINK GWAS analysis.
        Adapted from the reference blink_file_generator, with one key fix:
          - output_vcf is written as UNCOMPRESSED VCF (-Ov) to avoid downstream format issues.
        """
        train_samples_str = ",".join(train_samples)

        os.makedirs(os.path.join(self.gwas_output_dir, f"{file_fold_npc_index}"), exist_ok=True)
        current_dir = os.path.join(self.gwas_output_dir, f"{file_fold_npc_index}")

        phenotype_df = self._load_phenotypes().reset_index()
        phenotype_df[self.phenotype_column], _ = encode_phenotype_series(
            phenotype_df[self.phenotype_column],
            log_prefix=f"GWAS phenotype ({file_fold_npc_index})",
        )
        filtered_phenotype_df = phenotype_df[phenotype_df["taxa"].isin(train_samples)]
        phenotype_file = os.path.join(current_dir, f"train_{file_fold_npc_index}.txt")
        filtered_phenotype_df.to_csv(phenotype_file, sep="\t", index=False)

        vcf_file = self.vcf_file_path
        output_vcf = os.path.join(current_dir, f"train_{file_fold_npc_index}.vcf")

        mapping_file = os.path.join(current_dir, "chr_mapping.txt")
        if not os.path.exists(mapping_file):
            self.chromosome_mapping[["OriginalName", "RecodedNumber"]].to_csv(
                mapping_file,
                sep="\t",
                header=False,
                index=False,
            )

        original_chromosomes_to_keep = self.chromosome_mapping["OriginalName"].tolist()
        regions_str = ",".join(map(str, original_chromosomes_to_keep))

        compressed_vcf = str(_ensure_bgzip_and_index(Path(vcf_file)))

        if shutil.which("bcftools") is None:
            raise FileNotFoundError("bcftools is required for step 4 but was not found on PATH.")

        command = (
            f'bcftools view -s "{train_samples_str}" -r "{regions_str}" "{compressed_vcf}" | '
            f'bcftools annotate --rename-chrs "{mapping_file}" -Ov -o "{output_vcf}"'
        )
        subprocess.run(command, shell=True, check=True, executable="/bin/bash")

        if self.use_pcs and n_pcs > 0:
            principal_use = principal_df.iloc[:, :(n_pcs + 1)]
            filtered_principal_df = principal_use[principal_use["taxa"].isin(train_samples)]
            principal_file = os.path.join(current_dir, f"train_{file_fold_npc_index}.cov")
            filtered_principal_df.to_csv(principal_file, sep="\t", index=False)

    def run_gwas(self, file_fold_npc_index: str):
        current_dir = os.path.join(self.gwas_output_dir, f"{file_fold_npc_index}")
        os.makedirs(current_dir, exist_ok=True)
        tag = f"train_{file_fold_npc_index}"

        blink_command = (
            f"{self.blink_bin} "
            f"--gwas "
            f"--file {tag} "
            f"--vcf "
            f"--out {tag}"
        )

        log_path = os.path.join(current_dir, f"blink_run_{tag}.log")
        try:
            with open(log_path, "w") as log_file:
                subprocess.run(
                    blink_command,
                    shell=True,
                    check=True,
                    cwd=current_dir,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                )
        except subprocess.CalledProcessError as exc:
            log_tail = ""
            try:
                with open(log_path, "r") as log_file:
                    lines = log_file.readlines()[-20:]
                log_tail = "".join(lines).strip()
            except OSError:
                log_tail = ""
            if log_tail:
                raise RuntimeError(
                    f"BLINK failed for {file_fold_npc_index}. See {log_path}. Log tail:\n{log_tail}"
                ) from exc
            raise RuntimeError(
                f"BLINK failed for {file_fold_npc_index}. See {log_path}."
            ) from exc

    def run_fold(
        self,
        file_fold_index: str,
        genotype_npz: str,
        gwas_pcs_override: Optional[str] = None,
        pca_components_override: Optional[int] = None,
    ):
        train_numeric = _load_numeric_npz(genotype_npz)
        train_samples = train_numeric.index.astype(str).tolist()

        gwas_pcs_list = (
            parse_int_list(gwas_pcs_override)
            if (gwas_pcs_override and gwas_pcs_override.strip())
            else self.gwas_n_pcs_list
        )
        if not gwas_pcs_list:
            gwas_pcs_list = [0]

        n_components = self._resolve_pca_components(
            gwas_pcs_list,
            pca_components_override=pca_components_override,
        )

        principal_df, explained_variance_ratio = self.perform_pca(
            train_numeric,
            n_components=n_components,
        )
        self._write_pca_outputs(file_fold_index, principal_df, explained_variance_ratio)

        for pc in gwas_pcs_list:
            file_fold_npc_index = f"{file_fold_index}_{pc}PCs"
            self.blink_file_generator(train_samples, file_fold_npc_index, principal_df, pc)
            self.run_gwas(file_fold_npc_index)


def run_gwas_pca_from_config(
    config: dict,
    *,
    vcf_path: Path,
    pheno_path: Path,
    pheno_col: str,
    gwas_pcs_override: str = "",
    pca_components_override: Optional[int] = None,
) -> None:
    from Prediction_model.snp_numeric_transformer import SNP_numerical
    from Utils.utils import _save_numeric_npz

    vcf_path = Path(vcf_path).expanduser().resolve()
    pheno_path = Path(pheno_path).expanduser().resolve()

    if not vcf_path.exists():
        raise FileNotFoundError(f"vcf_file_path not found: {vcf_path}")
    if not pheno_path.exists():
        raise FileNotFoundError(f"phenotype_csv_path not found: {pheno_path}")

    pheno = pd.read_csv(pheno_path, sep=None, engine="python")
    if "taxa" not in pheno.columns:
        raise ValueError("phenotype_csv_path must contain a 'taxa' column.")
    if pheno_col not in pheno.columns:
        raise ValueError(f"phenotype column not found: {pheno_col}")

    pheno = pheno[["taxa", pheno_col]].dropna().set_index("taxa")
    vcf_samples = extract_vcf_samples(vcf_path)
    common = [s for s in vcf_samples if s in pheno.index]
    if not common:
        raise ValueError("No overlapping samples between VCF and phenotype CSV.")

    prediction_outdir = _resolve_outdir(config, key="prediction_output_dir", default="results", resolve=True)
    tmp_root = Path(prediction_outdir) / ".tmp" / "prediction_gwas"
    tmp_root.mkdir(parents=True, exist_ok=True)
    npz_path = tmp_root / "train_numeric.npz"

    config = dict(config)
    config.update(
        vcf_file_path=str(vcf_path),
        phenotype_csv_path=str(pheno_path),
        phenotype_column=pheno_col,
    )

    snp_tool = SNP_numerical(config)
    gwas_tool = GWAS_PCA(config)

    genotype_df = snp_tool.vcf_to_dataframe(str(vcf_path))
    present = [s for s in common if s in genotype_df.index]
    if not present:
        raise ValueError("No requested train samples were found in the VCF.")

    genotype_numeric = snp_tool.convert_genotypes_to_numeric(genotype_df.loc[present])
    _save_numeric_npz(genotype_numeric, str(npz_path))

    pca_components = pca_components_override
    if pca_components is None:
        raw_pca_components = (config.get("pca_n_components") or "").strip()
        if raw_pca_components:
            try:
                pca_components = int(raw_pca_components)
            except ValueError:
                pca_components = None

    gwas_tool.run_fold(
        "all",
        str(npz_path),
        gwas_pcs_override=gwas_pcs_override,
        pca_components_override=pca_components,
    )


def run_global_pca_qc_from_config(
    config: dict,
    *,
    vcf_path: Path,
    pheno_path: Path,
    pheno_col: str,
    pca_components_override: Optional[int] = None,
) -> None:
    vcf_path = Path(vcf_path).expanduser().resolve()
    pheno_path = Path(pheno_path).expanduser().resolve()

    if not vcf_path.exists():
        raise FileNotFoundError(f"vcf_file_path not found: {vcf_path}")
    if not pheno_path.exists():
        raise FileNotFoundError(f"phenotype_csv_path not found: {pheno_path}")

    config = dict(config)
    config.update(
        vcf_file_path=str(vcf_path),
        phenotype_csv_path=str(pheno_path),
        phenotype_column=pheno_col,
    )

    gwas_tool = GWAS_PCA(config)
    gwas_tool.run_global_qc(pca_components_override=pca_components_override)


def run_blink_gwas_pca(
    config_file: str,
    *,
    global_qc: bool = False,
    pca_components: Optional[int] = None,
    gwas_pcs: str = "",
) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from Utils.utils import read_config

    config = read_config(config_file)

    vcf_path_raw = (config.get("vcf_file_path") or "").strip()
    pheno_path_raw = (config.get("phenotype_csv_path") or "").strip()
    pheno_col = (config.get("phenotype_column") or "Phenotype").strip()

    if not vcf_path_raw:
        raise ValueError("vcf_file_path is required for GWAS/PCA stage.")
    if not pheno_path_raw:
        raise ValueError("phenotype_csv_path is required for GWAS/PCA stage.")

    if global_qc:
        run_global_pca_qc_from_config(
            config,
            vcf_path=vcf_path_raw,
            pheno_path=pheno_path_raw,
            pheno_col=pheno_col,
            pca_components_override=pca_components,
        )
    else:
        run_gwas_pca_from_config(
            config,
            vcf_path=vcf_path_raw,
            pheno_path=pheno_path_raw,
            pheno_col=pheno_col,
            gwas_pcs_override=gwas_pcs,
            pca_components_override=pca_components,
        )
