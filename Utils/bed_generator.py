#!/usr/bin/env python3
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple
import gzip


def _open_text(path: Path):
    opn = gzip.open if str(path).endswith(".gz") else open
    return opn(path, "rt", errors="replace")


def fasta_lengths(fa: Path) -> Dict[str, int]:
    lens: Dict[str, int] = {}
    with _open_text(fa) as f:
        name = None
        size = 0
        for line in f:
            if line.startswith(">"):
                if name:
                    lens[name] = size
                name = line[1:].split()[0]
                size = 0
            else:
                size += len(line.strip())
        if name:
            lens[name] = size
    return lens


def gff_contigs_with_genes(gff: Path, include_biotype: Optional[str] = None) -> Set[str]:
    contigs: Set[str] = set()
    with _open_text(gff) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, feat = cols[0], cols[2]
            if feat not in {"gene"}:
                continue
            if include_biotype:
                attrs = {kv.split("=", 1)[0]: kv.split("=", 1)[1] for kv in cols[8].split(";") if "=" in kv}
                biotype = (
                    attrs.get("gene_biotype")
                )
                if biotype is not None and biotype != include_biotype:
                    continue
            contigs.add(seqid)
    return contigs


def build_capture_bed(
    fasta: Path,
    out_bed: Path,
    top_n: int = 0,
    chroms: Optional[Iterable[str]] = None,
    require_genes: bool = False,
    gff: Optional[Path] = None,
    gene_biotype: Optional[str] = None,
) -> Tuple[List[Tuple[str, int]], List[str], List[str]]:
    lens = fasta_lengths(fasta)
    if not lens:
        raise ValueError(f"No contigs found in FASTA: {fasta}")

    gene_contigs = None
    if require_genes:
        if not gff or not Path(gff).exists():
            raise ValueError("gff is required to filter contigs with genes.")
        gene_contigs = gff_contigs_with_genes(gff, gene_biotype)

    explicit: List[str] = []
    for raw in chroms or []:
        val = str(raw).strip()
        if val and val not in explicit:
            explicit.append(val)

    missing = [c for c in explicit if c not in lens]

    selected: Dict[str, int] = {}

    for name in explicit:
        if name in lens:
            selected.setdefault(name, lens[name])

    sorted_lens = None
    if top_n and top_n > 0:
        sorted_lens = sorted(lens.items(), key=lambda kv: (-kv[1], kv[0]))
        for name, length in sorted_lens[:top_n]:
            selected.setdefault(name, length)

    if not selected:
        if sorted_lens is None:
            sorted_lens = sorted(lens.items(), key=lambda kv: (-kv[1], kv[0]))
        selected = {name: length for name, length in sorted_lens}

    selected_order = list(selected.items())
    selected_set = set(selected)

    filtered: List[str] = []
    if gene_contigs is not None:
        filtered = sorted([name for name in selected_set if name not in gene_contigs])
        if filtered:
            selected_order = [(name, length) for name, length in selected_order if name in gene_contigs]
            selected_set = {name for name, _ in selected_order}

    out_bed.parent.mkdir(parents=True, exist_ok=True)
    with open(out_bed, "w") as out:
        for name, length in selected_order:
            if length <= 0:
                continue
            out.write(f"{name}\t0\t{length}\n")

    return selected_order, missing, filtered
