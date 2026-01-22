"""
Run SnpEff and extract HIGH/MODERATE impact variant IDs.

Designed to be imported and called from an existing pipeline.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict


@dataclass
class SnpEffConfig:
    vcf: Path
    snpeff_jar_dir: Path
    species: str
    outdir: Path
    prefix: str = "snpeff"
    memory_gb: int = 8


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _run_cmd(cmd: str, log_path: Path) -> None:
    _ensure_dir(log_path.parent)
    with open(log_path, "w", encoding="utf-8") as log:
        log.write(cmd + "\n")
        log.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=log, stderr=log, text=True)


def _extract_impact_ids(ann_vcf: Path) -> Dict[str, object]:
    impact_rank = {"HIGH": 0, "MODERATE": 1}
    impact_ids = {impact: set() for impact in impact_rank}
    total_variants = 0

    with open(ann_vcf, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            total_variants += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            var_id = parts[2]
            info = parts[7]
            ann_field = None
            for entry in info.split(";"):
                if entry.startswith("ANN="):
                    ann_field = entry[4:]
                    break
            if not ann_field:
                continue
            ann_entries = ann_field.split(",")
            best_impact = None
            for ann in ann_entries:
                fields = ann.split("|")
                if len(fields) <= 2:
                    continue
                impact = fields[2].strip().upper()
                if impact not in impact_rank:
                    continue
                if best_impact is None or impact_rank[impact] < impact_rank[best_impact]:
                    best_impact = impact
                    if best_impact == "HIGH":
                        break
            if best_impact:
                if var_id == ".":
                    raise ValueError(f"{best_impact}-impact variant missing ID in VCF")
                impact_ids[best_impact].add(var_id)

    high_ids = impact_ids["HIGH"]
    moderate_ids = impact_ids["MODERATE"]
    return {
        "impact_ids": impact_ids,
        "total_variants": total_variants,
        "highimpact_variants": len(high_ids),
        "moderateimpact_variants": len(moderate_ids),
    }


def run_snpeff_and_extract_highimpact(cfg: SnpEffConfig) -> Dict[str, str]:
    outdir = Path(cfg.outdir).resolve()
    _ensure_dir(outdir)
    logs_dir = outdir / "logs"
    _ensure_dir(logs_dir)

    ann_vcf = outdir / "snpeff.ann.vcf"
    stats_html = outdir / "snpEff_summary.html"
    stats_csv = outdir / f"{cfg.prefix}.varann_summary.csv"

    snpeff_jar = Path(cfg.snpeff_jar_dir) / "snpEff" / "snpEff.jar"
    if not snpeff_jar.exists():
        raise FileNotFoundError(f"snpEff.jar not found: {snpeff_jar}")

    cmd = (
        f"java -Xmx{int(cfg.memory_gb)}g -jar {snpeff_jar} -v "
        f"-stats {stats_html} "
        f"-csvStats {stats_csv} "
        f"{cfg.species} {cfg.vcf} > {ann_vcf}"
    )
    _run_cmd(cmd, logs_dir / "snpeff.log")

    info = _extract_impact_ids(ann_vcf)
    impact_ids = info["impact_ids"]
    highimpact_ids = sorted(impact_ids.get("HIGH", []))
    moderateimpact_ids = sorted(impact_ids.get("MODERATE", []))

    highimpact_snplist = outdir / "highimpact.snplist"
    with open(highimpact_snplist, "w", encoding="utf-8") as handle:
        for vid in highimpact_ids:
            handle.write(f"{vid}\n")

    moderateimpact_snplist = outdir / "moderateimpact.snplist"
    with open(moderateimpact_snplist, "w", encoding="utf-8") as handle:
        for vid in moderateimpact_ids:
            handle.write(f"{vid}\n")

    impact_tsv = outdir / "impact_variants.tsv"
    with open(impact_tsv, "w", encoding="utf-8") as handle:
        handle.write("variant_id\timpact\n")
        for vid in highimpact_ids:
            handle.write(f"{vid}\tHIGH\n")
        for vid in moderateimpact_ids:
            handle.write(f"{vid}\tMODERATE\n")

    counts_tsv = outdir / "highimpact_counts.tsv"
    with open(counts_tsv, "w", encoding="utf-8") as handle:
        handle.write("metric\tcount\n")
        handle.write(f"total_variants\t{info['total_variants']}\n")
        handle.write(f"highimpact_variants\t{info['highimpact_variants']}\n")
        handle.write(f"highimpact_unique_ids\t{len(highimpact_ids)}\n")
        handle.write(f"moderateimpact_variants\t{info['moderateimpact_variants']}\n")
        handle.write(f"moderateimpact_unique_ids\t{len(moderateimpact_ids)}\n")
        handle.write(f"high_moderate_unique_ids\t{len(set(highimpact_ids) | set(moderateimpact_ids))}\n")

    return {
        "snpeff_ann_vcf": str(ann_vcf),
        "highimpact_snplist": str(highimpact_snplist),
        "moderateimpact_snplist": str(moderateimpact_snplist),
        "impact_variants_tsv": str(impact_tsv),
        "counts_tsv": str(counts_tsv),
        "logs_dir": str(logs_dir),
    }


"""
Example usage (no CLI parsing):

from pathlib import Path
from snpeff_highimpact import SnpEffConfig, run_snpeff_and_extract_highimpact

cfg = SnpEffConfig(
    vcf=Path("/data/02_Variant_Calling/cohort.filtered.vcf.gz"),
    snpeff_jar_dir=Path("/opt/snpeff"),
    species="Oreochromis_niloticus",
    outdir=Path("/data/variant_annotation"),
    prefix="tilapia",
    memory_gb=16,
)

outputs = run_snpeff_and_extract_highimpact(cfg)
"""
