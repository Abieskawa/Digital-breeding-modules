#!/usr/bin/env python3
import sys
import gzip
import csv
import os
from bisect import bisect_left
from collections import defaultdict

# Increase CSV field size limit if needed
csv.field_size_limit(sys.maxsize)

def parse_model_filename(filename):
    """
    Parse (model, fold, pc, n_snps) from a filename like:
      "ShapTop10_XGBoost_Fold_3-6_PC0_20SNPs.csv"
    or a path like:
      "some/path/ShapTop10_XGBoost_Fold_3-6_PC0_20SNPs.csv".
    Adjust the logic if your filenames differ.
    """
    base_name = os.path.basename(filename)
    root, _ = os.path.splitext(base_name)
    parts = root.split('_')

    meta = {
        'model':  None,
        'fold':   None,
        'pc':     None,
        'n_snps': None
    }

    if len(parts) > 1:
        meta['model'] = parts[1]

    try:
        fold_idx = parts.index("Fold")
        meta['fold'] = parts[fold_idx + 1]
    except ValueError:
        pass

    for p in parts:
        if p.startswith('PC'):
            meta['pc'] = p.replace('PC', '')
    for p in parts:
        if p.endswith("SNPs"):
            meta['n_snps'] = p.replace('SNPs', '')
    return meta

def open_file(fname):
    """Open a file in text mode, supporting gz if needed."""
    if fname.endswith('.gz'):
        return gzip.open(fname, 'rt')
    return open(fname, 'r')

def load_nonzero_snps(model_csv, topx_option='all'):
    """
    Reads CSV with columns: feature, importance.
    Only keeps rows where importance != 0.
    If topx_option is an integer, only the top X unique importance scores (including ties) are retained.
    Returns a list of (feature, importance) tuples sorted in descending order by importance.
    """
    topx = None
    if topx_option.lower() != 'all':
        try:
            topx = int(topx_option)
        except ValueError:
            print(f"Warning: --topx must be an integer or 'all'. Got '{topx_option}'. Using 'all'.", file=sys.stderr)
            topx = None

    rows = []
    with open(model_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                imp = float(row['importance'])
            except (ValueError, KeyError):
                continue
            if imp == 0.0:
                continue
            feature = row['feature']
            rows.append((feature, imp))

    rows.sort(key=lambda x: x[1], reverse=True)

    if topx is None:
        return rows

    selected = []
    unique_importances = []
    for (feature, imp) in rows:
        if imp not in unique_importances:
            unique_importances.append(imp)
            if len(unique_importances) > topx:
                break
        if len(unique_importances) <= topx:
            selected.append((feature, imp))
    return selected

def load_vcf_by_ids(vcf_file, id_list):
    """
    For each variant ID in id_list, find (chrom, pos) from the VCF (3rd column = ID).
    Returns a dictionary mapping variant_id to (chrom, pos).
    """
    results = {}
    with open_file(vcf_file) as f:
        header_found = False
        chrom_idx = pos_idx = id_idx = None
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header_found = True
                headers = line.strip().split('\t')
                headers[0] = headers[0].lstrip('#')
                chrom_idx = headers.index('CHROM')
                pos_idx   = headers.index('POS')
                id_idx    = headers.index('ID')
                continue

            if not header_found:
                continue

            fields = line.strip().split('\t')
            var_id = fields[id_idx]
            if var_id in id_list:
                chrom = fields[chrom_idx]
                pos   = int(fields[pos_idx])
                results[var_id] = (chrom, pos)
    return results

def load_gff(gff_file, molas=False):
    """
    Load features from a GFF file.

    In standard mode (molas=False) only gene features are loaded.
    In MOLAS mode (molas=True) only mRNA features are loaded, and the mRNA label is cleaned by:
      - Removing the "rna-" prefix (if present)
      - Removing any version suffix (e.g. turning "rna-XM_019367679.2" into "XM_019367679")
    
    In standard mode the gene-level labels are cleaned by stripping a leading "gene-" from both
    the gene_id and gene_name.

    Returns a list of dictionaries with keys: chrom, start, end, strand, gene_id, gene_name.
    """
    features = []
    with open_file(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2].lower()
            if molas:
                if feature_type != 'mrna':
                    continue
            else:
                if feature_type != 'gene':
                    continue

            start = int(fields[3])
            end   = int(fields[4])
            attr_dict = {}
            for attr in fields[8].split(';'):
                if '=' in attr:
                    k, v = attr.split('=', 1)
                    attr_dict[k.strip()] = v.strip()
            
            if molas:
                # Get the mRNA label from "Name" or fallback to "ID"
                name_val = attr_dict.get('Name', attr_dict.get('ID', 'NA'))
                if name_val.startswith("rna-"):
                    name_val = name_val[4:]
                clean_label = name_val.split('.')[0]
                features.append({
                    'chrom':     fields[0],
                    'start':     start,
                    'end':       end,
                    'strand':    fields[6],
                    'gene_id':   clean_label,
                    'gene_name': clean_label
                })
            else:
                gene_id = attr_dict.get('ID', 'NA')
                gene_name = attr_dict.get('Name', 'NA')
                # Remove "gene-" prefix if present
                if gene_id.startswith("gene-"):
                    gene_id = gene_id[5:]
                if gene_name.startswith("gene-"):
                    gene_name = gene_name[5:]
                features.append({
                    'chrom':     fields[0],
                    'start':     start,
                    'end':       end,
                    'strand':    fields[6],
                    'gene_id':   gene_id,
                    'gene_name': gene_name
                })
    return features

def group_genes_for_updown_search(genes):
    """
    Create two sorted lists for each chromosome:
      - by_end: sorted by gene.end
      - by_start: sorted by gene.start
    Returns a dictionary mapping chromosome to a dictionary with keys 'by_end' and 'by_start'.
    Each element in these lists is a tuple (coordinate, gene_dict).
    """
    chrom_map = defaultdict(lambda: {'by_end': [], 'by_start': []})
    for g in genes:
        c = g['chrom']
        chrom_map[c]['by_end'].append((g['end'], g))
        chrom_map[c]['by_start'].append((g['start'], g))

    for c in chrom_map:
        chrom_map[c]['by_end'].sort(key=lambda x: x[0])
        chrom_map[c]['by_start'].sort(key=lambda x: x[0])

    return chrom_map

def find_nearest_genes(chrom_map, chrom, pos):
    """
    1) If pos is inside one or more genes, return those gene dictionaries.
    2) Otherwise, return up to two dictionaries for the upstream and downstream genes:
       - Upstream: gene with the largest end < pos
       - Downstream: gene with the smallest start > pos
    """
    if chrom not in chrom_map:
        return []

    inrange = []
    for (_, gene_dict) in chrom_map[chrom]['by_start']:
        if gene_dict['start'] <= pos <= gene_dict['end']:
            inrange.append(gene_dict)
    if inrange:
        return inrange

    end_list = chrom_map[chrom]['by_end']
    ends_only = [x[0] for x in end_list]
    idx = bisect_left(ends_only, pos)
    upstream = None
    if idx > 0:
        cand_end, cand_gene = end_list[idx - 1]
        if cand_end < pos:
            upstream = cand_gene

    start_list = chrom_map[chrom]['by_start']
    starts_only = [x[0] for x in start_list]
    jdx = bisect_left(starts_only, pos)
    downstream = None
    if jdx < len(start_list):
        cand_start, cand_gene = start_list[jdx]
        if cand_start > pos:
            downstream = cand_gene

    out = []
    if upstream is not None:
        out.append(upstream)
    if downstream is not None:
        out.append(downstream)
    return out

def annotate_variants(
    model_csv: str,
    vcf: str,
    gff: str,
    *,
    output: str = "",
    topx: str = "all",
    molas: bool = False,
):
    meta = parse_model_filename(model_csv)
    model   = meta.get('model',  '-')
    fold    = meta.get('fold',   '-')
    pc_val  = meta.get('pc',     '-')
    n_snps  = meta.get('n_snps', '-')

    print(f"Reading model CSV: {model_csv}", file=sys.stderr)
    rows = load_nonzero_snps(model_csv, topx_option=topx)
    if not rows:
        print("No non-zero importance SNPs found after filtering.", file=sys.stderr)
        return []

    features = [r[0] for r in rows]

    print("Loading VCF...", file=sys.stderr)
    id_to_pos = load_vcf_by_ids(vcf, features)
    if not id_to_pos:
        print("No matching IDs found in the VCF.", file=sys.stderr)
        return []

    print("Loading GFF...", file=sys.stderr)
    genes = load_gff(gff, molas=molas)
    print("Indexing genes/mRNA entries...", file=sys.stderr)
    chrom_map = group_genes_for_updown_search(genes)

    out_handle = open(output, 'w') if output else sys.stdout

    print("Variant_ID\tmodel\tfold\tpc\tn_snps\tGene_ID\tGene_Name\tChrom\tGene_Start\tGene_End\tStrand\timportance",
          file=out_handle)

    all_genes_for_printing = []

    for (feature, importance) in rows:
        if feature not in id_to_pos:
            continue

        chrom, pos = id_to_pos[feature]
        neighbors = find_nearest_genes(chrom_map, chrom, pos)

        if not neighbors:
            print(f"{feature}\t{model}\t{fold}\t{pc_val}\t{n_snps}\t-\t-\t-\t-\t-\t-\t{importance}",
                  file=out_handle)
            continue

        for gene_dict in neighbors:
            print(
                f"{feature}\t{model}\t{fold}\t{pc_val}\t{n_snps}\t"
                f"{gene_dict['gene_id']}\t{gene_dict['gene_name']}\t"
                f"{gene_dict['chrom']}\t{gene_dict['start']}\t{gene_dict['end']}\t{gene_dict['strand']}\t"
                f"{importance}",
                file=out_handle
            )
            all_genes_for_printing.append(gene_dict['gene_id'])

    if output:
        out_handle.close()

    print("\nDone annotating. Gene IDs found:", file=sys.stderr)
    unique_genes = sorted(set(all_genes_for_printing))
    for gid in unique_genes:
        print(gid, file=sys.stderr)
    return unique_genes
