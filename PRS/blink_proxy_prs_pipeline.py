"""
BLINK → clump → LD proxy check → impact-variant listing.

Designed to be imported and called from an existing pipeline.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


@dataclass
class PipelineConfig:
    blink_tsv: Path
    gwas_vcf: Path
    full_vcf: Path
    highimpact_list: Path
    outdir: Path
    plink2_path: str = "plink2"
    top_n: int = 100
    clump_kb: int = 250
    clump_r2: float = 0.2
    tag_kb: int = 250
    tag_r2: float = 0.2
    sample_check_n: int = 10
    threads: Optional[int] = None


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _run_cmd(cmd: List[str], log_path: Path) -> None:
    _ensure_dir(log_path.parent)
    with open(log_path, "w", encoding="utf-8") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
        subprocess.run(cmd, check=True, stdout=log, stderr=log, text=True)


def _plink2_base(cfg: PipelineConfig) -> List[str]:
    cmd = [cfg.plink2_path]
    if cfg.threads is not None:
        cmd += ["--threads", str(int(cfg.threads))]
    return cmd


def _convert_vcf_to_pfile(cfg: PipelineConfig, vcf: Path, out_prefix: Path, log_path: Path) -> None:
    cmd = _plink2_base(cfg) + [
        "--vcf",
        str(vcf),
        "--max-alleles",
        "2",
        "--allow-extra-chr",
        "--double-id",
        "--make-pgen",
        "--out",
        str(out_prefix),
    ]
    _run_cmd(cmd, log_path)



def _read_pvar_ids(cfg: PipelineConfig, pfile_prefix: Path, tmp_dir: Path, log_path: Path) -> set:
    pvar_path = pfile_prefix.with_suffix(".pvar")
    if pvar_path.exists():
        df = pd.read_csv(pvar_path, sep="\t", usecols=["ID"], dtype=str)
        return set(df["ID"].dropna().tolist())

    pvar_zst = pfile_prefix.with_suffix(".pvar.zst")
    if not pvar_zst.exists():
        raise FileNotFoundError(f"pvar not found for prefix: {pfile_prefix}")

    decompressed = tmp_dir / f"{pfile_prefix.name}.pvar"
    _decompress_zst(cfg, pvar_zst, decompressed, log_path)
    df = pd.read_csv(decompressed, sep="\t", usecols=["ID"], dtype=str)
    return set(df["ID"].dropna().tolist())


def _decompress_zst(cfg: PipelineConfig, zst_path: Path, out_path: Path, log_path: Path) -> Path:
    out_prefix = out_path
    cmd = _plink2_base(cfg) + ["--zst-decompress", str(zst_path), "--out", str(out_prefix)]
    _run_cmd(cmd, log_path)

    if out_path.exists():
        return out_path

    fallback = zst_path.with_suffix("")
    if fallback.exists():
        shutil.copy2(fallback, out_path)
        return out_path

    raise FileNotFoundError(f"Failed to decompress {zst_path}")


def _read_clumps(cfg: PipelineConfig, clump_prefix: Path, tmp_dir: Path, log_path: Path) -> pd.DataFrame:
    clumps_path = clump_prefix.with_suffix(".clumps")
    if clumps_path.exists():
        return pd.read_csv(clumps_path, delim_whitespace=True, comment="#", dtype=str)

    clumps_zst = clump_prefix.with_suffix(".clumps.zst")
    if clumps_zst.exists():
        decompressed = tmp_dir / f"{clump_prefix.name}.clumps"
        _decompress_zst(cfg, clumps_zst, decompressed, log_path)
        shutil.copy2(decompressed, clumps_path)
        return pd.read_csv(decompressed, delim_whitespace=True, comment="#", dtype=str)

    raise FileNotFoundError(f"No clumps output found for prefix: {clump_prefix}")


def _read_vcor(cfg: PipelineConfig, vcor_prefix: Path, tmp_dir: Path, log_path: Path) -> pd.DataFrame:
    vcor_path = vcor_prefix.with_suffix(".vcor")
    if vcor_path.exists():
        return pd.read_csv(vcor_path, delim_whitespace=True, comment="#", dtype=str)

    vcor_zst = vcor_prefix.with_suffix(".vcor.zst")
    if vcor_zst.exists():
        decompressed = tmp_dir / f"{vcor_prefix.name}.vcor"
        _decompress_zst(cfg, vcor_zst, decompressed, log_path)
        shutil.copy2(decompressed, vcor_path)
        return pd.read_csv(decompressed, delim_whitespace=True, comment="#", dtype=str)

    raise FileNotFoundError(f"No vcor output found for prefix: {vcor_prefix}")


def _parse_vcor_columns(df: pd.DataFrame) -> Tuple[str, str, str]:
    col_a = None
    col_b = None
    if "ID_A" in df.columns and "ID_B" in df.columns:
        col_a, col_b = "ID_A", "ID_B"
    elif "ID1" in df.columns and "ID2" in df.columns:
        col_a, col_b = "ID1", "ID2"
    elif "SNP_A" in df.columns and "SNP_B" in df.columns:
        col_a, col_b = "SNP_A", "SNP_B"

    if col_a is None:
        raise ValueError("Unable to locate ID_A/ID_B columns in vcor output")

    r2_col = None
    for cand in df.columns:
        if cand.upper().startswith("R2"):
            r2_col = cand
            break
    if r2_col is None:
        raise ValueError("Unable to locate R2 column in vcor output")

    return col_a, col_b, r2_col


def _load_impact_variants(path: Path) -> Dict[str, str]:
    """Load impacts from TSV/CSV (variant_id + impact) or a 1-column list (assumed HIGH)."""
    if not path.exists():
        raise FileNotFoundError(f"Impact list not found: {path}")

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        first = handle.readline().strip().lower()
    if not first:
        return {}
    tokens = first.replace(",", " ").replace("\t", " ").split()
    has_header = any(
        token in {"impact", "variant_id", "variant", "snp", "id"}
        for token in tokens
    )

    df = pd.read_csv(path, sep=None, engine="python", dtype=str, header=0 if has_header else None)
    if df.empty:
        return {}

    if has_header:
        cols_lower = {str(c).lower(): c for c in df.columns}
        id_col = None
        for key in ("variant_id", "id", "variant", "snp"):
            if key in cols_lower:
                id_col = cols_lower[key]
                break

        impact_col = None
        for key in ("impact", "ann_impact"):
            if key in cols_lower:
                impact_col = cols_lower[key]
                break

        if id_col is None:
            return {}
    else:
        id_col = 0
        impact_col = 1 if df.shape[1] > 1 else None

    if impact_col is None:
        ids = df[id_col].dropna().astype(str).tolist()
        return {vid.strip(): "HIGH" for vid in ids if str(vid).strip()}

    mapping: Dict[str, str] = {}
    for vid, impact in zip(df[id_col], df[impact_col]):
        if vid is None or impact is None:
            continue
        vid = str(vid).strip()
        if not vid:
            continue
        imp = str(impact).strip().upper()
        if imp not in {"HIGH", "MODERATE"}:
            continue
        if imp == "HIGH" or vid not in mapping:
            mapping[vid] = imp
    return mapping


def run_pipeline(cfg: PipelineConfig) -> Dict[str, str]:
    outdir = Path(cfg.outdir).resolve()
    logs_dir = outdir / "logs"
    tmp_dir = outdir / "tmp"
    plink_dir = outdir / "plink"
    _ensure_dir(logs_dir)
    _ensure_dir(tmp_dir)
    _ensure_dir(plink_dir)

    gwas_prefix = plink_dir / "gwas"
    full_prefix = plink_dir / "full"

    # Step 1: convert VCFs to pfiles
    _convert_vcf_to_pfile(cfg, Path(cfg.gwas_vcf), gwas_prefix, logs_dir / "01_make_pgen_gwas.log")
    _convert_vcf_to_pfile(cfg, Path(cfg.full_vcf), full_prefix, logs_dir / "02_make_pgen_full.log")

    # Step 2: load BLINK results and create sumstats
    blink_df = pd.read_csv(cfg.blink_tsv, sep="\t", dtype=str)
    if not {"taxa", "p_value"}.issubset(blink_df.columns):
        raise ValueError("BLINK TSV must have columns: taxa, p_value")
    blink_df["p_value"] = pd.to_numeric(blink_df["p_value"], errors="coerce")
    blink_df = blink_df.groupby("taxa", as_index=False)["p_value"].min()
    if blink_df.empty:
        raise ValueError("BLINK TSV contains no usable records after parsing")
    sumstats_path = outdir / "blink.sumstats.tsv"
    sumstats_df = blink_df[["taxa", "p_value"]].rename(columns={"taxa": "SNP", "p_value": "P"})
    sumstats_df = sumstats_df.dropna(subset=["P"])
    sumstats_df.to_csv(sumstats_path, sep="\t", index=False)

    # Step 3: fail-fast ID check
    gwas_ids = _read_pvar_ids(cfg, gwas_prefix, tmp_dir, logs_dir / "03_read_pvar_gwas.log")
    full_ids = _read_pvar_ids(cfg, full_prefix, tmp_dir, logs_dir / "04_read_pvar_full.log")
    check_count = min(int(cfg.sample_check_n), len(blink_df))
    sample_ids = blink_df["taxa"].head(check_count).tolist()
    missing_gwas = [s for s in sample_ids if s not in gwas_ids]
    missing_full = [s for s in sample_ids if s not in full_ids]
    if missing_gwas or missing_full:
        raise ValueError(
            "BLINK SNP IDs do not match pvar IDs. "
            f"Missing in GWAS pvar: {missing_gwas}. "
            f"Missing in FULL pvar: {missing_full}."
        )

    # Step 4: clumping on GWAS genotype set
    clump_prefix = outdir / "blink_clump"
    clump_cmd = _plink2_base(cfg) + [
        "--pfile",
        str(gwas_prefix),
        "--clump",
        str(sumstats_path),
        "--clump-id-field",
        "SNP",
        "--clump-p-field",
        "P",
        "--clump-kb",
        str(int(cfg.clump_kb)),
        "--clump-r2",
        str(cfg.clump_r2),
        "--clump-p1",
        "1",
        "--clump-p2",
        "1",
        "--clump-unphased",
        "--out",
        str(clump_prefix),
    ]
    _run_cmd(clump_cmd, logs_dir / "05_clump.log")

    clumps_df = _read_clumps(cfg, clump_prefix, tmp_dir, logs_dir / "06_read_clumps.log")
    if "SNP" in clumps_df.columns:
        lead_ids = clumps_df["SNP"].dropna().unique().tolist()
    elif "ID" in clumps_df.columns:
        lead_ids = clumps_df["ID"].dropna().unique().tolist()
    else:
        raise ValueError("Unable to find SNP column in clumps output")
    if not lead_ids:
        raise ValueError("No lead SNPs found in clumps output")

    # Step 5: rank leads and select top_n
    lead_df = pd.DataFrame({"lead_snp": lead_ids})
    lead_df = lead_df.merge(
        blink_df[["taxa", "p_value"]].rename(columns={"taxa": "lead_snp"}),
        on="lead_snp",
        how="left",
    )
    if lead_df["p_value"].isna().any():
        missing = lead_df[lead_df["p_value"].isna()]["lead_snp"].head(5).tolist()
        raise ValueError(f"Lead SNPs missing p_value in BLINK results: {missing}")
    lead_df = lead_df.sort_values(by="p_value", ascending=True)
    lead_df["rank"] = range(1, len(lead_df) + 1)
    lead_ranked_path = outdir / "lead_ranked.tsv"
    lead_df.to_csv(lead_ranked_path, sep="\t", index=False)

    top_n = min(int(cfg.top_n), len(lead_df))
    top_leads = lead_df.head(top_n)
    lead_snplist_path = outdir / "lead.snplist"
    top_leads["lead_snp"].to_csv(lead_snplist_path, index=False, header=False)

    # Step 6: LD proxy check vs impact variants (HIGH/MODERATE)
    ld_prefix = outdir / "lead_ld"
    ld_cmd = _plink2_base(cfg) + [
        "--pfile",
        str(full_prefix),
        "--r2-unphased",
        "--ld-snp-list",
        str(lead_snplist_path),
        "--ld-window-kb",
        str(int(cfg.tag_kb)),
        "--ld-window-r2",
        str(cfg.tag_r2),
        "--out",
        str(ld_prefix),
    ]
    _run_cmd(ld_cmd, logs_dir / "07_ld_proxy.log")

    vcor_df = _read_vcor(cfg, ld_prefix, tmp_dir, logs_dir / "08_read_vcor.log")
    id_a_col, id_b_col, r2_col = _parse_vcor_columns(vcor_df)

    impact_map = _load_impact_variants(Path(cfg.highimpact_list))
    if not impact_map:
        print("Warning: impact list is empty; no HIGH/MODERATE variants will be reported.")
    lead_set = set(top_leads["lead_snp"].tolist())
    proxy_map: Dict[str, List[Tuple[str, float, str]]] = {lead: [] for lead in lead_set}

    for _, row in vcor_df.iterrows():
        id_a = row[id_a_col]
        id_b = row[id_b_col]
        if id_a in lead_set:
            lead_id = id_a
            partner = id_b
        elif id_b in lead_set:
            lead_id = id_b
            partner = id_a
        else:
            continue
        impact = impact_map.get(partner)
        if impact:
            r2_val = float(row[r2_col])
            proxy_map[lead_id].append((partner, r2_val, impact))

    summary_rows = []
    impact_rows = []
    for _, rec in top_leads.iterrows():
        lead_id = rec["lead_snp"]
        rank = int(rec["rank"])
        blink_p = rec["p_value"]
        proxies = proxy_map.get(lead_id, [])
        lead_impact = impact_map.get(lead_id)
        if lead_impact:
            proxies = proxies + [(lead_id, 1.0, lead_impact)]

        proxy_impacts: Dict[str, str] = {}
        proxy_r2: Dict[str, float] = {}
        for pid, r2_val, impact in proxies:
            if pid not in proxy_r2 or r2_val > proxy_r2[pid]:
                proxy_r2[pid] = r2_val
            if impact == "HIGH" or pid not in proxy_impacts:
                proxy_impacts[pid] = impact

        high_ids = sorted([pid for pid, imp in proxy_impacts.items() if imp == "HIGH"])
        moderate_ids = sorted([pid for pid, imp in proxy_impacts.items() if imp == "MODERATE"])
        has_impact = 1 if proxy_impacts else 0
        has_high = 1 if high_ids else 0
        max_r2_high = max([proxy_r2[pid] for pid in high_ids], default=0.0)
        max_r2_any = max(proxy_r2.values(), default=0.0)
        summary_rows.append(
            {
                "lead_snp": lead_id,
                "rank": rank,
                "blink_p": blink_p,
                "has_highimpact_proxy": has_high,
                "has_impact_proxy": has_impact,
                "num_highimpact_proxies": len(high_ids),
                "num_moderate_proxies": len(moderate_ids),
                "max_r2_among_highimpact_proxies": max_r2_high,
                "max_r2_among_impact_proxies": max_r2_any,
                "highimpact_proxy_ids": ",".join(high_ids),
                "moderateimpact_proxy_ids": ",".join(moderate_ids),
            }
        )

        for pid, impact in proxy_impacts.items():
            impact_rows.append(
                {
                    "lead_snp": lead_id,
                    "rank": rank,
                    "blink_p": blink_p,
                    "proxy_variant": pid,
                    "impact": impact,
                    "r2": proxy_r2.get(pid, 0.0),
                    "is_lead_variant": 1 if pid == lead_id else 0,
                }
            )

    summary_df = pd.DataFrame(summary_rows)
    summary_path = outdir / "lead_impact_proxy_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    impact_df = pd.DataFrame(
        impact_rows,
        columns=[
            "lead_snp",
            "rank",
            "blink_p",
            "proxy_variant",
            "impact",
            "r2",
            "is_lead_variant",
        ],
    )
    impact_list_path = outdir / "impact_variants_near_leads.tsv"
    impact_df.to_csv(impact_list_path, sep="\t", index=False)
    if impact_df.empty:
        print("No HIGH/MODERATE impact variants found near top GWAS leads.")

    return {
        "sumstats": str(sumstats_path),
        "clumps": str(clump_prefix.with_suffix(".clumps")),
        "lead_snplist": str(lead_snplist_path),
        "lead_ranked": str(lead_ranked_path),
        "lead_vcor": str(ld_prefix.with_suffix(".vcor")),
        "lead_impact_summary": str(summary_path),
        "impact_variants_near_leads": str(impact_list_path),
        "logs_dir": str(logs_dir),
    }


"""
Example usage (no CLI parsing):

from pathlib import Path
from blink_proxy_prs_pipeline import PipelineConfig, run_pipeline

cfg = PipelineConfig(
    blink_tsv=Path("/data/GWAS/train_1_0PCs_Phenotype_GWAS_result.txt"),
    gwas_vcf=Path("/data/GWAS/train_1_0PCs.vcf"),
    full_vcf=Path("/data/02_Variant_Calling/cohort.filtered.vcf.gz"),
    highimpact_list=Path("/data/impact_variants.tsv"),
    outdir=Path("/data/PRS"),
    plink2_path="/opt/conda/envs/DGbreeding-biotools/bin/plink2",
)

outputs = run_pipeline(cfg)
"""
