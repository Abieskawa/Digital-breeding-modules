#!/usr/bin/env python3
"""
STAR ⟩ StringTie versatile RNA seq pipeline

Outputs (you may request any combination):
  --transcript-fpkm    per-sample transcript FPKM values
  --gene-fpkm          per-sample gene FPKM values
  --gene-tpm           per-sample gene TPM values
  --transcript-counts  transcript count matrix via prepDE.py per sample
  --gene-counts        gene count matrix via prepDE.py per sample
  --map-only           perform mapping and filtering only; skip StringTie and prepDE
  --merge-expression   merge and process per-mode tables into a single expression table (FPKM/TPM/counts)
Each output mode writes into its own directory under --out-dir; per-sample count matrices are placed at --out-dir.
"""
import argparse, subprocess, glob, os, sys, math

# ---------- helpers ----------

def compute_genomeSAindex_nbases(fasta):
    total = 0
    with open(fasta) as fh:
        for ln in fh:
            if ln.startswith('>'): continue
            total += len(ln.strip())
    val = int(math.log2(total) / 2 - 1)
    return val if val < 14 else 14


def run(cmd, *, cwd=None, use_shell=False, **kw):
    if use_shell:
        print(cmd, flush=True)
        subprocess.run(cmd, shell=True, check=True, cwd=cwd,
                       executable='/bin/bash', **kw)
    else:
        print(' '.join(cmd), flush=True)
        subprocess.run(cmd, check=True, cwd=cwd, **kw)


def run_version(cmd):
    print('#', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=False)

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description='STAR → StringTie multi-mode pipeline')
    ap.add_argument('--reads-dir', required=True)
    ap.add_argument('--genome-fasta', required=True)
    ap.add_argument('--annotation', required=True,
                    help='GTF annotation file (GFF input is deprecated)')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--map-only', action='store_true',
                    help='Only perform mapping and filtering; skip StringTie and prepDE')
    ap.add_argument('--transcript-fpkm',   action='store_true')
    ap.add_argument('--gene-fpkm',         action='store_true')
    ap.add_argument('--gene-tpm',          action='store_true')
    ap.add_argument('--transcript-counts', action='store_true')
    ap.add_argument('--gene-counts',       action='store_true')
    ap.add_argument('--merge-expression',  action='store_true')
    args = ap.parse_args()

    print('Command invoked:', ' '.join(sys.argv), flush=True)

    if args.annotation.lower().endswith(('.gff', '.gff3')):
        sys.exit('Error: GFF input is not suggested to run STAR; please provide a GTF (*.gtf) file.')

    if not args.map_only and not any([args.transcript_fpkm, args.gene_fpkm,
                                      args.gene_tpm, args.transcript_counts,
                                      args.gene_counts]):
        sys.exit('Error: select at least one output mode or use --map-only')

    star_index    = os.path.join(args.out_dir, 'STAR_index')
    align_dir     = os.path.join(args.out_dir, 'alignments'); os.makedirs(align_dir, exist_ok=True)
    stringtie_dir = os.path.join(args.out_dir, 'stringtie'); os.makedirs(stringtie_dir, exist_ok=True)

    mode_dirs = {}
    if args.transcript_fpkm:   mode_dirs['transcript_fpkm']   = os.path.join(args.out_dir, 'transcript_fpkm')
    if args.transcript_counts: mode_dirs['transcript_counts'] = os.path.join(args.out_dir, 'transcript_counts')
    if args.gene_fpkm:         mode_dirs['gene_fpkm']         = os.path.join(args.out_dir, 'gene_fpkm')
    if args.gene_tpm:          mode_dirs['gene_tpm']          = os.path.join(args.out_dir, 'gene_tpm')
    if args.gene_counts:       mode_dirs['gene_counts']       = os.path.join(args.out_dir, 'gene_counts')
    for d in mode_dirs.values(): os.makedirs(d, exist_ok=True)

    for tool in [['STAR','--version'], ['stringtie','--version'], ['samtools','--version'], ['sort','--version']]:
        run_version(tool)
    print(flush=True)

    idx_file = os.path.join(star_index, 'SAindex')
    if os.path.exists(idx_file):
        print(f"STAR index exists.", flush=True)
    if not os.path.exists(idx_file):
        cmd = [
            'STAR','--runMode','genomeGenerate',
            '--genomeFastaFiles', args.genome_fasta,
            '--genomeSAindexNbases', str(compute_genomeSAindex_nbases(args.genome_fasta)),
            '--runThreadN', str(args.threads),
            '--genomeDir', star_index,
            '--sjdbGTFfile', args.annotation
        ]
        print('Building STAR index with command:', ' '.join(cmd), flush=True)
        os.makedirs(star_index, exist_ok=True)
        run(cmd)
    else:
        print("STAR index exists, skipping build", flush=True)

    # process each sample (support both .fastq* and .fq* extensions)
    fastq1 = glob.glob(os.path.join(args.reads_dir, '*1.cleaned.fastq*'))
    fq1    = glob.glob(os.path.join(args.reads_dir, '*1.cleaned.fq*'))
    r1_files = sorted(set(fastq1 + fq1))
    samples = []
    for r1 in r1_files:
        r2 = r1.replace('1.cleaned', '2.cleaned')
        paired = os.path.exists(r2)
        sample = os.path.basename(r1).split('1.cleaned')[0]    
        if not paired:
            print(f"Single-end sample for {sample}, no R2 found: {r2}", flush=True)
        samples.append(sample)

        pref = os.path.join(align_dir, sample + '.')
        bam  = pref + 'Aligned.sortedByCoord.out.bam'
        if not os.path.exists(bam):
            cmd = [
                'STAR','--runMode','alignReads','--twopassMode','Basic',
                '--runThreadN', str(args.threads),
                '--genomeDir', star_index,
                '--outSAMtype','BAM','SortedByCoordinate',
                '--quantMode','TranscriptomeSAM','GeneCounts',
                '--outFileNamePrefix', pref,
                '--outSAMstrandField','intronMotif'
            ]
            if paired:
                cmd += ['--readFilesIn', r1, r2]
            else:
                cmd += ['--readFilesIn', r1]

            if r1.endswith('.gz'):
                cmd += ['--readFilesCommand','zcat']
            print(f"Running STAR alignment for {sample} with command:", ' '.join(cmd), flush=True)
            run(cmd)
        else:
            print(f"Skipping STAR for {sample}, BAM exists", flush=True)

        bam_flag2 = pref + 'flag2.bam'
        samtools_cmd = ['samtools','view','-b','-f','2','-@', str(args.threads), bam, '-o', bam_flag2]
        print('Filtering BAM with command:', ' '.join(samtools_cmd), flush=True)
        run(samtools_cmd)
        
        if args.map_only:
            continue

        # run StringTie
        raw_dir = os.path.join(stringtie_dir, sample); os.makedirs(raw_dir, exist_ok=True)
        cmd = [
            'stringtie', bam_flag2, '-G', args.annotation,
            '-o', os.path.join(raw_dir,'stringtie.gtf'),
            '-p', str(args.threads), '-e', '-l', sample,
            '-b', raw_dir, '-f','0.15','-m','200','-a','10',
            '-j','1','-c','2','-g','50','-M','0.95',
            '-A', os.path.join(raw_dir,'abundance.tab')
        ]
        run(cmd)

        if args.transcript_fpkm:
            src = os.path.join(raw_dir,'t_data.ctab')
            out = os.path.join(mode_dirs['transcript_fpkm'], sample + '.transcript_fpkm.tab')
            run(f"sort -k6,6 {src} | cut -f6,12 > {out}", use_shell=True)
        if args.gene_fpkm:
            src = os.path.join(raw_dir,'abundance.tab')
            out = os.path.join(mode_dirs['gene_fpkm'], sample + '.gene_fpkm.tab')
            run(f"sort -k1,1 {src} | cut -f1,8 > {out}", use_shell=True)
        if args.gene_tpm:
            src = os.path.join(raw_dir,'abundance.tab')
            out = os.path.join(mode_dirs['gene_tpm'], sample + '.gene_tpm.tab')
            run(f"sort -k1,1 {src} | cut -f1,9 > {out}", use_shell=True)

    # generate count matrices
    if args.transcript_counts or args.gene_counts:
        mapping_file = os.path.join(args.out_dir, 'prepDE_input.txt')
        with open(mapping_file, 'w') as fh:
            for sample in samples:
                gtf_file = os.path.abspath(os.path.join(stringtie_dir, sample, 'stringtie.gtf'))
                fh.write(f"{sample}\t{gtf_file}\n")
                print('Mapping entry:', sample, '->', gtf_file, flush=True)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        prepde_script = os.path.join(script_dir, 'prepDE.py')
        prepde_cmd = ['python3', '-W', 'ignore::SyntaxWarning', prepde_script, '-i', mapping_file]
        if args.gene_counts:
            gene_matrix = os.path.join(mode_dirs['gene_counts'], 'gene_count_matrix.csv')
            prepde_cmd += ['-g', gene_matrix]
        if args.transcript_counts:
            trans_matrix = os.path.join(mode_dirs['transcript_counts'], 'transcript_count_matrix.csv')
            prepde_cmd += ['-t', trans_matrix]
        print('Running prepDE command:', ' '.join(prepde_cmd), flush=True)
        try:
            run(prepde_cmd)
            print("Successfully generated count matrices", flush=True)
            # create sorted copies
            if args.gene_counts:
                sorted_gene_counts = os.path.join(mode_dirs['gene_counts'], 'gene_count_matrix.sorted.csv')
                run(f"head -n1 {gene_matrix} > {sorted_gene_counts}; tail -n+2 {gene_matrix} | sort -t, -k1,1 >> {sorted_gene_counts}", use_shell=True)
            if args.transcript_counts:
                sorted_trans_counts = os.path.join(mode_dirs['transcript_counts'], 'transcript_count_matrix.sorted.csv')
                run(f"head -n1 {trans_matrix} > {sorted_trans_counts}; tail -n+2 {trans_matrix} | sort -t, -k1,1 >> {sorted_trans_counts}", use_shell=True)
        except subprocess.CalledProcessError as e:
            sys.exit(e.returncode)

    # merge per-option outputs into tables
    if args.merge_expression:
        # transcript FPKM
        if args.transcript_fpkm:
            suffix = '.t.fpkm'
            ids, expr = [], {}
            for sample in samples:
                fn = os.path.join(mode_dirs['transcript_fpkm'], f"{sample}.transcript_fpkm.tab")
                if os.path.exists(fn):
                    with open(fn) as fh:
                        first = fh.readline()
                        fields = first.strip().split('\t')
                        is_data = False
                        for v in fields[1:]:
                            try:
                                float(v)
                                is_data = True
                                break
                            except ValueError:
                                pass
                        lines = [first] + fh.readlines() if is_data else fh.readlines()
                    for ln in lines:
                        id_, v = ln.strip().split('\t')
                        expr.setdefault(id_, {})[sample+suffix] = v
                        if id_ not in ids: ids.append(id_)
            out = os.path.join(args.out_dir, 'transcript_fpkm.table')
            with open(out, 'w') as fh:
                fh.write('ID\t' + '\t'.join([s+suffix for s in samples]) + '\n')
                for id_ in sorted(ids):
                    fh.write(id_ + '\t' + '\t'.join(expr[id_].get(s+suffix, '.') for s in samples) + '\n')
        # gene FPKM
        if args.gene_fpkm:
            suffix = '.g.fpkm'
            ids, expr = [], {}
            for sample in samples:
                fn = os.path.join(mode_dirs['gene_fpkm'], f"{sample}.gene_fpkm.tab")
                if os.path.exists(fn):
                    with open(fn) as fh:
                        first = fh.readline()
                        fields = first.strip().split('\t')
                        is_data = False
                        for v in fields[1:]:
                            try:
                                float(v)
                                is_data = True
                                break
                            except ValueError:
                                pass
                        lines = [first] + fh.readlines() if is_data else fh.readlines()
                    for ln in lines:
                        id_, v = ln.strip().split('\t')
                        expr.setdefault(id_, {})[sample+suffix] = v
                        if id_ not in ids: ids.append(id_)
            out = os.path.join(args.out_dir, 'gene_fpkm.table')
            with open(out, 'w') as fh:
                fh.write('ID\t' + '\t'.join([s+suffix for s in samples]) + '\n')
                for id_ in sorted(ids):
                    fh.write(id_ + '\t' + '\t'.join(expr[id_].get(s+suffix, '.') for s in samples) + '\n')
        # gene TPM
        if args.gene_tpm:
            suffix = '.g.tpm'
            ids, expr = [], {}
            for sample in samples:
                fn = os.path.join(mode_dirs['gene_tpm'], f"{sample}.gene_tpm.tab")
                if os.path.exists(fn):
                    with open(fn) as fh:
                        first = fh.readline()
                        fields = first.strip().split('\t')
                        is_data = False
                        for v in fields[1:]:
                            try:
                                float(v)
                                is_data = True
                                break
                            except ValueError:
                                pass
                        lines = [first] + fh.readlines() if is_data else fh.readlines()
                    for ln in lines:
                        id_, v = ln.strip().split('\t')
                        expr.setdefault(id_, {})[sample+suffix] = v
                        if id_ not in ids: ids.append(id_)
            out = os.path.join(args.out_dir, 'gene_tpm.table')
            with open(out, 'w') as fh:
                fh.write('ID\t' + '\t'.join([s+suffix for s in samples]) + '\n')
                for id_ in sorted(ids):
                    fh.write(id_ + '\t' + '\t'.join(expr[id_].get(s+suffix, '.') for s in samples) + '\n')  # ensure sort by first column
        # gene counts
        if args.gene_counts and 'sorted_gene_counts' in locals():
            gene_cols, gene_data = [], {}
            with open(sorted_gene_counts) as fh:
                hdr = next(fh).strip().split(',')[1:]
                for sample in hdr:
                    gene_cols.append(sample+'.g.count')
                for ln in fh:
                    parts = ln.strip().split(',')
                    id_ = parts[0]
                    if id_.startswith('<class'): continue
                    for sample, val in zip(hdr, parts[1:]):
                        gene_data.setdefault(id_, {})[sample+'.g.count'] = val
            out = os.path.join(args.out_dir, 'gene_counts.table')
            with open(out, 'w') as fh:
                fh.write('ID\t' + '\t'.join(gene_cols) + '\n')
                for id_ in sorted(gene_data):
                    fh.write(id_ + '\t' + '\t'.join(gene_data[id_].get(c) or '.' for c in gene_cols) + '\n')
        # transcript counts
        if args.transcript_counts and 'sorted_trans_counts' in locals():
            t_cols, t_data = [], {}
            with open(sorted_trans_counts) as fh:
                hdr = next(fh).strip().split(',')[1:]
                for sample in hdr:
                    t_cols.append(sample+'.t.count')
                for ln in fh:
                    parts = ln.strip().split(',')
                    id_ = parts[0]
                    for sample, val in zip(hdr, parts[1:]):
                        t_data.setdefault(id_, {})[sample+'.t.count'] = val
            out = os.path.join(args.out_dir, 'transcript_counts.table')
            with open(out, 'w') as fh:
                fh.write('ID\t' + '\t'.join(t_cols) + '\n')
                for id_ in sorted(t_data):
                    fh.write(id_ + '\t' + '\t'.join(t_data[id_].get(c) or '.' for c in t_cols) + '\n')

    print('Everything Finished.')

if __name__=='__main__':
    main()
