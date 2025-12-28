#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import os
import subprocess as sbp
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from pycirclize import Circos

from Utils.utils import call_log, time_stamp


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


def make_circos_for_vcf(
    fasta_lens: Dict[str, int],
    gene_starts: Dict[str, List[int]],
    vcf_pos: Dict[str, List[int]],
    window: int,
    out_png: Path,
    gene_window: Optional[int] = None,
    var_window: Optional[int] = None,
    tick_step: Optional[int] = None,
    gene_cap_override: Optional[float] = None,
    var_cap_override: Optional[float] = None,
    dpi: int = 300,
) -> None:
    sectors = {k: v for k, v in fasta_lens.items() if v > 0}
    if not sectors:
        return

    binned: Dict[str, Dict] = {}
    global_gene_max = 0
    global_var_max = 0
    gwin = gene_window if gene_window is not None else window
    vwin = var_window if var_window is not None else window
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
        if gene_vals.size:
            global_gene_max = max(global_gene_max, int(gene_vals.max()))
        if var_vals.size:
            global_var_max = max(global_var_max, int(var_vals.max()))

    n_sectors = len(sectors)
    space = 5 if n_sectors <= 10 else 3 if n_sectors <= 22 else 1
    circos = Circos(sectors, space=space)

    chr_names = list(sectors.keys())
    palette = ["#4c72b0", "#55a868", "#c44e52", "#8172b3", "#ccb974", "#64b5cd"]
    chr_name2color = {ch: palette[i % len(palette)] for i, ch in enumerate(chr_names)}

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

        var_cap = var_cap_override if (var_cap_override and var_cap_override > 0) else float(global_var_max or 1)
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

        gene_cap = gene_cap_override if (gene_cap_override and gene_cap_override > 0) else float(global_gene_max or 1)
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

    out_png.parent.mkdir(parents=True, exist_ok=True)
    circos.savefig(str(out_png), dpi=int(dpi))


def fasta_lengths(fa: Path) -> Dict[str, int]:
    lens: Dict[str, int] = {}
    opn = gzip.open if str(fa).endswith(".gz") else open
    with opn(fa, "rt", errors="replace") as f:
        name = None
        size = 0
        for ln in f:
            if ln.startswith(">"):
                if name:
                    lens[name] = size
                name = ln[1:].split()[0]
                size = 0
            else:
                size += len(ln.strip())
        if name:
            lens[name] = size
    return lens


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
                or attrs.get("biotype")
                or attrs.get("gene_type")
                or attrs.get("transcript_biotype")
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


def prepare_circos_reference(
    ref_fasta: Optional[Path],
    gff: Optional[Path],
    gff_biotype: Optional[str],
    target_seq: Optional[int],
    allowed_chroms: Optional[Set[str]],
) -> Tuple[bool, Dict[str, int], Dict[str, List[int]], Dict[str, str], List[Tuple[str, int]]]:
    if not (ref_fasta and Path(ref_fasta).exists()):
        return False, {}, {}, {}, []

    fasta_lens_raw = fasta_lengths(ref_fasta)
    if target_seq:
        ordered = sorted(fasta_lens_raw.items(), key=lambda kv: kv[1], reverse=True)[: target_seq or 0]
    elif allowed_chroms:
        ordered = [(k, v) for k, v in fasta_lens_raw.items() if k in allowed_chroms]
    else:
        ordered = sorted(fasta_lens_raw.items(), key=lambda kv: kv[1], reverse=True)

    if not ordered:
        return False, {}, {}, {}, []

    rename_map = {name: f"chr{i+1}" for i, (name, _) in enumerate(ordered)}
    fasta_lens = {rename_map[name]: length for name, length in ordered}

    gene_raw = gff_gene_starts(gff, gff_biotype) if gff and Path(gff).exists() else defaultdict(list)
    gene_starts = {rename_map[name]: gene_raw.get(name, []) for name, _ in ordered}

    return True, fasta_lens, gene_starts, rename_map, ordered


def _append_log_line(log_dir: Path, step: str, msg: str) -> Path:
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{step}.log"
    with open(log_path, "a", encoding="utf-8") as fh:
        fh.write(msg.rstrip() + "\n")
    return log_path


def _optional_int(val, default=None, none_if_blank: bool = True):
    if val in (None, "", False):
        return None if none_if_blank else default
    try:
        iv = int(val)
    except (TypeError, ValueError):
        return default
    if none_if_blank and iv == 0:
        return None
    return iv


class VariantCircosEvaluator:
    """
    Step 3: Variant stats + Circos (gene + variant density).
    """

    def _setup_local_env(self) -> None:
        tmp_dir = self.outdir / ".tmp"
        mpl_dir = self.outdir / ".cache" / "matplotlib"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        mpl_dir.mkdir(parents=True, exist_ok=True)
        for key in ("TMPDIR", "TEMP", "TMP"):
            os.environ[key] = str(tmp_dir)
        os.environ["MPLCONFIGDIR"] = str(mpl_dir)
        os.environ.setdefault("MPLBACKEND", "Agg")
        self.tmp_dir = tmp_dir
        self.mpl_config_dir = mpl_dir

    def __init__(self, args):
        self.threads   = int(getattr(args, "threads", 1))
        self.outdir    = Path(getattr(args, "outdir", ".")).resolve()
        self.variant_dir = Path(getattr(args, "variant_outdir", self.outdir)).resolve()
        self.ref_fasta = Path(getattr(args, "ref_fasta", "")).expanduser().resolve() if getattr(args, "ref_fasta", None) else None
        self.gff       = Path(getattr(args, "gff", "")).expanduser().resolve() if getattr(args, "gff", None) else None

        self._setup_local_env()

        self.gff_biotype   = getattr(args, "gff_biotype", None)
        self.circos_window = _optional_int(getattr(args, "circos_window", None), default=10000, none_if_blank=False) or 10000
        self.circos_gene_window = _optional_int(getattr(args, "circos_gene_window", None), default=None, none_if_blank=True)
        self.circos_var_window  = _optional_int(getattr(args, "circos_var_window",  None), default=None, none_if_blank=True)
        self.circos_tick_step = _optional_int(getattr(args, "circos_tick_step", None), default=None, none_if_blank=True)
        self.circos_chroms = list(getattr(args, "circos_chroms", []) or [])
        try:
            gcap = getattr(args, "circos_gene_max", None)
            self.circos_gene_max = float(gcap) if gcap not in (None, "", False) else None
        except Exception:
            self.circos_gene_max = None
        try:
            vcap = getattr(args, "circos_var_max", None)
            self.circos_var_max = float(vcap) if vcap not in (None, "", False) else None
        except Exception:
            self.circos_var_max = None
        try:
            self.target_seq  = int(getattr(args, "target_seq", 0))
        except Exception:
            self.target_seq  = 0

        self.report_dir = self.outdir / "evaluation"
        self.report_dir.mkdir(parents=True, exist_ok=True)
        (self.report_dir / "logs").mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _bcftools_stats_one(vcf: Path, stats_out: Path) -> Optional[Dict]:
        stats_out.parent.mkdir(parents=True, exist_ok=True)
        cmd = f"bcftools stats -s - {vcf} > {stats_out}"
        try:
            sbp.run(cmd, shell=True, check=True, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)
        except sbp.CalledProcessError:
            call_log(stats_out.parent, "bcftools_stats", cmd)
            return None

    def evaluate_variants(self) -> pd.DataFrame:
        if not self.variant_dir.exists():
            time_stamp(f"[variants] variant_dir not found: {self.variant_dir}")
            return pd.DataFrame()

        allowed_chroms = set(self.circos_chroms) if self.circos_chroms else None
        max_chroms = self.target_seq if self.target_seq and self.target_seq > 0 else None

        do_circos, fasta_lens, gene_starts, rename_map, contig_order = prepare_circos_reference(
            ref_fasta=self.ref_fasta,
            gff=self.gff,
            gff_biotype=self.gff_biotype,
            target_seq=max_chroms,
            allowed_chroms=allowed_chroms,
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
        if do_circos and contig_order:
            map_path = circos_dir / "contig_name_map.tsv"
            with open(map_path, "w") as mh:
                mh.write("original\tproxy\tlength\n")
                for name, length in contig_order:
                    mh.write(f"{name}\t{rename_map.get(name, name)}\t{length}\n")
            time_stamp(f"[variants] circos contig map: {map_path}")

        for vcf in vcfs:
            stem = vcf.name[:-7] if str(vcf).endswith(".vcf.gz") else vcf.stem
            rec = self._bcftools_stats_one(vcf, stats_dir / f"{stem}.stats")
            if rec:
                rows.append(rec)

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
                gw = self.circos_gene_window or self.circos_window
                vw = self.circos_var_window or self.circos_window
                time_stamp(f"[variants] circos params: tick_step={self.circos_tick_step}, gene_window={gw}, var_window={vw}")
                make_circos_for_vcf(
                    fasta_lens=fasta_lens,
                    gene_starts=gene_starts,
                    vcf_pos=vpos,
                    window=self.circos_window,
                    gene_window=self.circos_gene_window or self.circos_window,
                    var_window=self.circos_var_window or self.circos_window,
                    out_png=out_png,
                    tick_step=self.circos_tick_step,
                    gene_cap_override=self.circos_gene_max,
                    var_cap_override=self.circos_var_max,
                )
                time_stamp(f"[variants] circos: {out_png}")

        df = pd.DataFrame(rows)
        out_csv = self.report_dir / "variant_summary.csv"
        if not df.empty:
            df.to_csv(out_csv, index=False)
            time_stamp(f"[variants] wrote {out_csv}")
        return df
