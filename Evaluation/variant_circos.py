#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from pycirclize import Circos

from Utils.utils import (
    _append_log_line,
    _load_chromosome_recode,
    _optional_int,
    _resolve_outdir,
    fasta_lengths,
    setup_eval_env,
    time_stamp,
)


def _format_bp(x: int) -> str:
    if x >= 1_000_000_000:
        return f"{x / 1e9:.0f} Gb"
    if x >= 1_000_000:
        return f"{x / 1e6:.0f} Mb"
    if x >= 1_000:
        return f"{x / 1e3:.0f} kb"
    return f"{x} bp"


def _fixed_ticks(seq_len: int, step_bp: Optional[int] = None):
    """Tick positions/labels (0..L) using a user-specified fixed step. Always include 0."""
    if seq_len <= 0:
        return [], []
    if step_bp and step_bp > 0:
        step = max(1, int(step_bp))
        pos = list(range(0, seq_len, step))
    else:
        pos = [0, seq_len]
    labels = [_format_bp(p) for p in pos]
    if labels:
        labels[0] = "0 Mb"
    return pos, labels


def _add_track_yticks(track, tick_max: float) -> None:
    if tick_max <= 0:
        return
    ticks = [0.0, tick_max / 2.0, tick_max]
    labels = [str(int(round(t))) for t in ticks]
    yticks = getattr(track, "yticks", None)
    if yticks is None:
        return
    try:
        yticks(ticks, labels=labels, label_size=7, tick_length=0.8)
    except TypeError:
        try:
            yticks(ticks, labels=labels)
        except TypeError:
            try:
                yticks(ticks)
            except Exception:
                return


def make_circos_for_vcf(
    fasta_lens: Dict[str, int],
    gene_starts: Dict[str, List[int]],
    vcf_pos: Dict[str, List[int]],
    window: int,
    out_png: Path,
    tick_step: Optional[int] = None,
    dpi: int = 300,
) -> None:
    sectors = {k: v for k, v in fasta_lens.items() if v > 0}
    if not sectors:
        return

    binned: Dict[str, Dict] = {}
    gwin = window
    vwin = window
    for chrom, L in sectors.items():
        nbin_g = max(1, int(np.ceil(L / gwin)))
        edges_g = np.linspace(0, L, nbin_g + 1)
        centers_g = (edges_g[:-1] + edges_g[1:]) / 2.0

        gene_vals = np.zeros(nbin_g, dtype=int)
        if chrom in gene_starts and gene_starts[chrom]:
            gene_vals, _ = np.histogram(
                np.asarray(gene_starts[chrom], dtype=int),
                bins=nbin_g,
                range=(1, L),
            )

        nbin_v = max(1, int(np.ceil(L / vwin)))
        edges_v = np.linspace(0, L, nbin_v + 1)
        centers_v = (edges_v[:-1] + edges_v[1:]) / 2.0

        var_vals = np.zeros(nbin_v, dtype=int)
        if chrom in vcf_pos and vcf_pos[chrom]:
            var_vals, _ = np.histogram(
                np.asarray(vcf_pos[chrom], dtype=int),
                bins=nbin_v,
                range=(1, L),
            )

        binned[chrom] = dict(
            gene_centers=centers_g,
            gene_vals=gene_vals,
            gene_width=gwin,
            var_centers=centers_v,
            var_vals=var_vals,
            var_width=vwin,
            L=L,
        )

    gene_maxes = [int(v["gene_vals"].max()) for v in binned.values() if v["gene_vals"].size]
    var_maxes = [int(v["var_vals"].max()) for v in binned.values() if v["var_vals"].size]
    global_gene_max = max(gene_maxes) if gene_maxes else 0
    global_var_max = max(var_maxes) if var_maxes else 0

    n_sectors = len(sectors)
    space = 5 if n_sectors <= 10 else 3 if n_sectors <= 22 else 1
    circos = Circos(sectors, space=space)

    chr_names = list(sectors.keys())
    palette = ["#4c72b0", "#55a868", "#c44e52", "#8172b3", "#ccb974", "#64b5cd"]
    chr_name2color = {ch: palette[i % len(palette)] for i, ch in enumerate(chr_names)}

    first_sector = circos.sectors[0] if circos.sectors else None
    for sector in circos.sectors:
        chrom = sector.name
        data = binned.get(chrom)
        if data is None:
            continue

        L = data["L"]
        gene_centers = data["gene_centers"]
        gene_vals = data["gene_vals"]
        gene_width = data["gene_width"]
        var_centers = data["var_centers"]
        var_vals = data["var_vals"]
        var_width = data["var_width"]

        ideogram = sector.add_track((92, 100), r_pad_ratio=0.02)
        color = chr_name2color.get(chrom, "lightgrey")
        ideogram.axis(fc=color, ec="none", alpha=0.9)

        mid = 0.5 * (sector.deg_lim[0] + sector.deg_lim[1])
        ideogram.text(
            chrom,
            r=96,
            size=10,
            fontweight="bold",
            fontstyle="italic",
            rotation=mid + 90,
            rotation_mode="anchor",
            va="center",
            ha="center",
            zorder=50,
        )

        tick_outer = sector.add_track((100, 100), r_pad_ratio=0.0)
        tick_pos, tick_labs = _fixed_ticks(L, step_bp=tick_step)
        if tick_labs and tick_pos:
            tick_outer.xticks(
                tick_pos,
                labels=tick_labs,
                label_size=6,
                label_orientation="vertical",
                tick_length=0.9,
            )

        var_cap = float(global_var_max or 1)
        var_track = sector.add_track((80, 92), r_pad_ratio=0.05)
        var_track.axis(fc="#f5f5f5", ec="none", alpha=1.0)
        if var_vals.any():
            var_track.bar(
                var_centers,
                var_vals.astype(float),
                vmin=0.0,
                vmax=float(var_cap),
                width=float(var_width),
                color="#de425b",
                ec="none",
            )
        if sector == first_sector:
            _add_track_yticks(var_track, var_cap)

        gene_cap = float(global_gene_max or 1)
        gene_track = sector.add_track((70, 80), r_pad_ratio=0.05)
        gene_track.axis(fc="#f5f5f5", ec="none", alpha=1.0)
        if gene_vals.any():
            gene_track.bar(
                gene_centers,
                gene_vals.astype(float),
                vmin=0.0,
                vmax=float(gene_cap),
                width=float(gene_width),
                color="#4c9f50",
                ec="none",
            )
        if sector == first_sector:
            _add_track_yticks(gene_track, gene_cap)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    legend_handles = [
        Patch(facecolor="#de425b", label="Variant density"),
        Patch(facecolor="#4c9f50", label="Gene density"),
    ]
    if hasattr(circos, "plotfig"):
        fig = circos.plotfig()
        fig.legend(handles=legend_handles, loc="upper right", frameon=False)
        fig.savefig(str(out_png), dpi=int(dpi))
        plt.close(fig)
    else:
        circos.savefig(str(out_png), dpi=int(dpi))


def gff_gene_starts(gff: Path, include_biotype: Optional[str]) -> Dict[str, List[int]]:
    starts: Dict[str, List[int]] = defaultdict(list)
    opn = gzip.open if str(gff).endswith(".gz") else open
    with opn(gff, "rt", errors="replace") as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            seqid, feat, start_s = c[0], c[2], c[3]
            if feat not in {"gene", "mRNA", "transcript"}:
                continue
            try:
                start = int(start_s)
            except Exception:
                continue
            attrs = {kv.split("=", 1)[0]: kv.split("=", 1)[1] for kv in c[8].split(";") if "=" in kv}
            biotype = (
                attrs.get("gene_biotype")
            )
            if include_biotype and biotype is not None and biotype != include_biotype:
                continue
            starts[seqid].append(start)
    return starts


def vcf_positions(vcf: Path) -> Dict[str, List[int]]:
    pos: Dict[str, List[int]] = defaultdict(list)
    opn = gzip.open if str(vcf).endswith(".gz") else open
    with opn(vcf, "rt", errors="replace") as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            col = ln.rstrip("\n").split("\t")
            if len(col) < 2:
                continue
            try:
                pos[col[0]].append(int(col[1]))
            except Exception:
                pass
    return pos


def _write_chrom_stats(
    vcf: Path,
    out_stats: Path,
    *,
    contig_lengths: Dict[str, int],
    label_map: Optional[Dict[str, str]] = None,
    allowed_contigs: Optional[Set[str]] = None,
    contig_order: Optional[List[str]] = None,
    filtered_vcf: Optional[Path] = None,
) -> Dict[str, object]:
    counts: Dict[str, int] = defaultdict(int)
    total = 0
    out_stats.parent.mkdir(parents=True, exist_ok=True)

    opn = gzip.open if str(vcf).endswith(".gz") else open
    vcf_out = None
    if filtered_vcf:
        filtered_vcf.parent.mkdir(parents=True, exist_ok=True)
        vcf_out = open(filtered_vcf, "w")
    try:
        with opn(vcf, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    if vcf_out:
                        vcf_out.write(line)
                    continue
                chrom = line.split("\t", 1)[0]
                if allowed_contigs and chrom not in allowed_contigs:
                    continue
                counts[chrom] += 1
                total += 1
                if vcf_out:
                    vcf_out.write(line)
    finally:
        if vcf_out:
            vcf_out.close()

    if contig_order:
        order = [c for c in contig_order if not allowed_contigs or c in allowed_contigs]
    elif allowed_contigs:
        order = sorted(list(allowed_contigs))
    else:
        order = sorted(counts.keys())

    with open(out_stats, "w") as out:
        out.write(f"#total_snps={total}\n")
        out.write("chrom\tlabel\tlength_bp\tsnp_count\tsnp_density_per_mb\n")
        for chrom in order:
            if allowed_contigs and chrom not in allowed_contigs:
                continue
            snps = counts.get(chrom, 0)
            length = contig_lengths.get(chrom, 0)
            label = label_map.get(chrom, chrom) if label_map else chrom
            if length:
                density = snps * 1_000_000.0 / length
                out.write(f"{chrom}\t{label}\t{length}\t{snps}\t{density:.6f}\n")
            else:
                out.write(f"{chrom}\t{label}\tNA\t{snps}\tNA\n")

    return {"total_snps": total, "contigs": len(order)}


def prepare_circos_reference(
    ref_fasta: Optional[Path],
    gff: Optional[Path],
    gff_biotype: Optional[str],
    allowed_chroms: Optional[Set[str]],
    recode_map: Optional[Dict[str, str]] = None,
) -> Tuple[bool, Dict[str, int], Dict[str, List[int]], Dict[str, str], List[Tuple[str, int]]]:
    if not (ref_fasta and Path(ref_fasta).exists()):
        return False, {}, {}, {}, []

    fasta_lens_raw = fasta_lengths(ref_fasta)
    if allowed_chroms:
        ordered = [(k, v) for k, v in fasta_lens_raw.items() if k in allowed_chroms]
    else:
        ordered = sorted(fasta_lens_raw.items(), key=lambda kv: kv[1], reverse=True)

    if not ordered:
        return False, {}, {}, {}, []

    used_labels: Set[str] = set()
    rename_map: Dict[str, str] = {}
    for i, (name, _) in enumerate(ordered):
        if recode_map:
            candidate = recode_map.get(name) or name
        else:
            candidate = f"chr{i+1}"
        base_label = str(candidate)
        suffix = 1
        while candidate in used_labels:
            suffix += 1
            candidate = f"{base_label}_{suffix}"
        rename_map[name] = candidate
        used_labels.add(candidate)
    fasta_lens = {rename_map[name]: length for name, length in ordered}

    gene_raw = gff_gene_starts(gff, gff_biotype) if gff and Path(gff).exists() else defaultdict(list)
    gene_starts = {rename_map[name]: gene_raw.get(name, []) for name, _ in ordered}

    return True, fasta_lens, gene_starts, rename_map, ordered


class VariantCircosEvaluator:
    """
    Step 3: Variant stats + Circos (gene + variant density).
    """

    def __init__(self, args):
        self.outdir    = _resolve_outdir(base_outdir=getattr(args, "outdir", "."), resolve=True)
        self.variant_dir = _resolve_outdir(base_outdir=getattr(args, "variant_outdir", self.outdir), resolve=True)
        self.ref_fasta = Path(getattr(args, "ref_fasta", "")).expanduser().resolve() if getattr(args, "ref_fasta", None) else None
        self.gff       = Path(getattr(args, "gff", "")).expanduser().resolve() if getattr(args, "gff", None) else None
        chrom_raw = getattr(args, "chromosome_csv_path", "")
        if isinstance(chrom_raw, Path):
            chrom_raw = str(chrom_raw)
        chrom_raw = (chrom_raw or "").strip()
        self.chromosome_csv_path = Path(chrom_raw).expanduser().resolve() if chrom_raw else None

        _, _, self.report_dir = setup_eval_env(self.outdir)

        self.gff_biotype   = getattr(args, "gff_biotype", None)
        self.circos_window = _optional_int(getattr(args, "circos_window", None), default=10000, none_if_blank=False) or 10000
        self.circos_tick_step = _optional_int(getattr(args, "circos_tick_step", None), default=None, none_if_blank=True)
        # report_dir already created by setup_eval_env

    def evaluate_variants(self) -> pd.DataFrame:
        if not self.variant_dir.exists():
            time_stamp(f"[variants] variant_dir not found: {self.variant_dir}")
            return pd.DataFrame()

        recode_map = _load_chromosome_recode(self.chromosome_csv_path)
        ref_lengths = fasta_lengths(self.ref_fasta) if (self.ref_fasta and self.ref_fasta.exists()) else {}
        allowed_chroms = set(recode_map.keys()) if recode_map else None
        do_circos, fasta_lens, gene_starts, rename_map, contig_order = prepare_circos_reference(
            ref_fasta=self.ref_fasta,
            gff=self.gff,
            gff_biotype=self.gff_biotype,
            allowed_chroms=allowed_chroms,
            recode_map=recode_map,
        )
        if do_circos and not contig_order:
            time_stamp("[variants] circos skipped: no contigs after filtering")
            do_circos = False
            fasta_lens, gene_starts, rename_map = {}, {}, {}

        vcfs = sorted(list(self.variant_dir.glob("*.vcf"))) + sorted(list(self.variant_dir.glob("*.vcf.gz")))
        vcfs = [p for p in vcfs if ".g.vcf" not in p.name and not p.name.endswith(".g.vcf.gz")]
        vcfs = [p for p in vcfs if p.name.startswith("cohort")]
        if not vcfs:
            time_stamp(f"[variants] no cohort* VCF/VCF.GZ in {self.variant_dir}")
            return pd.DataFrame()

        rows: List[Dict] = []
        stats_dir = self.report_dir / "variant_stats"
        circos_dir = self.report_dir / "variant_circos"
        circos_dir.mkdir(parents=True, exist_ok=True)
        contigs_filter = set(recode_map.keys()) if recode_map else set()
        contig_order_names = [name for name, _ in contig_order] if contig_order else None
        if do_circos and contig_order:
            map_path = circos_dir / "contig_name_map.tsv"
            with open(map_path, "w") as mh:
                mh.write("original\tproxy\tlength\n")
                for name, length in contig_order:
                    mh.write(f"{name}\t{rename_map.get(name, name)}\t{length}\n")
            time_stamp(f"[variants] circos contig map: {map_path}")

        for vcf in vcfs:
            stem = vcf.name[:-7] if str(vcf).endswith(".vcf.gz") else vcf.stem
            summary = {"vcf": str(vcf)}
            summary.update(
                _write_chrom_stats(
                    vcf,
                    stats_dir / f"{stem}.stats",
                    contig_lengths=ref_lengths,
                    label_map=recode_map,
                    allowed_contigs=None,
                    contig_order=None,
                )
            )
            if contigs_filter:
                target_summary = _write_chrom_stats(
                    vcf,
                    stats_dir / f"{stem}.filtered_target_seq.stats",
                    contig_lengths=ref_lengths,
                    label_map=recode_map,
                    allowed_contigs=contigs_filter,
                    contig_order=contig_order_names,
                    filtered_vcf=stats_dir / f"{stem}.filtered_target_seq.vcf",
                )
                summary["total_snps_target_seq"] = target_summary.get("total_snps", 0)
            rows.append(summary)

            if do_circos:
                vpos_raw = vcf_positions(vcf)
                vpos = {}
                for k, v in vpos_raw.items():
                    if k in rename_map:
                        vpos[rename_map[k]] = v
                if not vpos:
                    time_stamp(f"[variants] circos skipped for {vcf}: no variants on selected chromosomes")
                    continue

                for cn in fasta_lens:
                    _append_log_line(
                        self.report_dir / "logs",
                        "pycirclize_circos",
                        f"{stem} contig {cn}: genes={len(gene_starts.get(cn, []))}, vars={len(vpos.get(cn, []))}",
                    )

                out_png = circos_dir / f"{stem}.circos.png"
                time_stamp(f"[variants] circos params: tick_step={self.circos_tick_step}, window={self.circos_window}")
                make_circos_for_vcf(
                    fasta_lens=fasta_lens,
                    gene_starts=gene_starts,
                    vcf_pos=vpos,
                    window=self.circos_window,
                    out_png=out_png,
                    tick_step=self.circos_tick_step,
                )
                time_stamp(f"[variants] circos: {out_png}")

        df = pd.DataFrame(rows)
        out_csv = self.report_dir / "variant_summary.csv"
        if not df.empty:
            df.to_csv(out_csv, index=False)
            time_stamp(f"[variants] wrote {out_csv}")
        return df
