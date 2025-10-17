#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import math
import re

def compute_genomeSAindex_nbases(genome_path):
    """
    Compute the STAR parameter genomeSAindexNbases based on the genome length.
    
    The value is calculated as:
        floor(log2(total_genome_length) / 2) - 1
    and then taking the minimum between 14 and that value.
    """
    total_length = 0
    with open(genome_path, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                continue
            total_length += len(line.strip())
    computed = int(math.log2(total_length) / 2) - 1
    return min(14, computed)

def find_fastq_file(directory, pattern):
    """
    Search for a FASTQ file in 'directory' that starts with 'pattern'
    and ends with one of the allowed extensions.
    Allowed extensions: .fastq, .fq, .fastq.gz, .fq.gz
    """
    allowed_ext = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    for ext in allowed_ext:
        candidate = os.path.join(directory, pattern + ext)
        if os.path.isfile(candidate):
            return candidate
    return None

def find_isoseq_file(directory, pattern):
    """
    Search for an iso-seq file in 'directory' that starts with 'pattern'
    and ends with '_polished.hq.fasta'.
    """
    candidate = os.path.join(directory, pattern + "_polished.hq.fasta")
    if os.path.isfile(candidate):
        return candidate
    return None

def extract_sample_basename(filename):
    """
    Extract the sample basename and an optional ".cleaned" suffix from a filename.
    
    This function ignores files with "fastqc" (case-insensitive) and supports:
      - Paired-end: e.g. a_muscle_02_1.cleaned.fq.gz or a_muscle_02_1.fq.gz
      - Single-end: e.g. a_muscle_02.cleaned.fq.gz or a_muscle_02.fq.gz
    
    Returns a tuple (sample, suffix) where suffix is either ".cleaned" or an empty string.
    """
    if "fastqc" in filename.lower():
        return None, None
    # Paired-end pattern: capture sample and optional ".cleaned"
    paired_pattern = re.compile(
        r"^(?P<sample>.+)_[12](?P<cleaned>\.cleaned)?\.(?:fq|fastq)(?:\.gz)?$", re.IGNORECASE)
    m = paired_pattern.match(filename)
    if m:
        sample = m.group("sample")
        suffix = m.group("cleaned") if m.group("cleaned") is not None else ""
        return sample, suffix
    # Single-end pattern: capture sample and optional ".cleaned"
    single_pattern = re.compile(
        r"^(?P<sample>.+)(?P<cleaned>\.cleaned)?\.(?:fq|fastq)(?:\.gz)?$", re.IGNORECASE)
    m = single_pattern.match(filename)
    if m:
        sample = m.group("sample")
        suffix = m.group("cleaned") if m.group("cleaned") is not None else ""
        return sample, suffix
    return None, None

def extract_isoseq_basename(filename):
    """
    Extract the sample basename from an iso-seq filename.
    
    Supports files ending with '_polished.hq.fasta'.
    Returns the sample name or None if pattern doesn't match.
    """
    isoseq_pattern = re.compile(r"^(?P<sample>.+)_polished\.hq\.fasta$", re.IGNORECASE)
    m = isoseq_pattern.match(filename)
    if m:
        return m.group("sample")
    return None

def run_with_conda_env(cmd, env_name):
    """
    Run a command within a specific conda environment.
    """
    conda_cmd = ["conda", "run", "-n", env_name] + cmd
    print("Running command:", " ".join(conda_cmd))
    subprocess.run(conda_cmd, check=True)

def process_short_read_mode(args):
    """
    Process samples in short-read mode using STAR alignment.
    """
    print("Running in SHORT-READ mode")
    
    # Generate STAR genome index if not present.
    star_index_path = os.path.join(args.genomedir, "SAindex")
    if os.path.isdir(args.genomedir) and os.path.exists(star_index_path):
        print("STAR genome index already generated. Skipping.")
    else:
        print("STAR genome index not found. Generating...")
        genome_saindex_nbases = compute_genomeSAindex_nbases(args.genomepath)
        star_cmd = [
            "STAR",
            "--runMode", "genomeGenerate",
            "--genomeFastaFiles", args.genomepath,
            "--genomeSAindexNbases", str(genome_saindex_nbases),
            "--runThreadN", str(args.threads),
            "--genomeDir", args.genomedir
        ]
        run_with_conda_env(star_cmd, "assembly_annotation")

    # List files in the working directory.
    all_files = os.listdir(args.wd)
    print("Files found in working directory:", all_files)
    
    # Build a dictionary mapping each unique sample to its detected suffix (".cleaned" or "").
    unique_samples = {}  # key: sample, value: suffix
    allowed_ext = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    for filename in all_files:
        if any(filename.lower().endswith(ext) for ext in allowed_ext):
            sample, suffix = extract_sample_basename(filename)
            if sample:
                # If multiple files for a sample exist with inconsistent naming, warn.
                if sample in unique_samples and unique_samples[sample] != suffix:
                    print(f"Warning: Inconsistent naming for sample {sample}.")
                unique_samples[sample] = suffix
    unique_sample_names = sorted(unique_samples.keys())
    print("Unique sample basenames found:", unique_sample_names)

    def get_read_files_command(file_path):
        # If file is gzipped, add zcat to uncompress on the fly.
        return ["--readFilesCommand", "zcat"] if file_path and file_path.lower().endswith(".gz") else []

    # Process each unique sample.
    for sample in unique_sample_names:
        suffix = unique_samples[sample]
        # Look for paired-end files using the detected suffix.
        paired1 = find_fastq_file(args.wd, f"{sample}_1{suffix}")
        paired2 = find_fastq_file(args.wd, f"{sample}_2{suffix}")
        # Alternatively, for single-end files:
        single  = find_fastq_file(args.wd, f"{sample}{suffix}")

        out_prefix   = os.path.join(args.out_dir, f"{sample}_2pass")
        log_progress = f"{out_prefix}Log.progress.out"
        sj_tab       = f"{out_prefix}SJ.out.tab"
        bam_sorted   = os.path.join(args.out_dir, f"{sample}_2passAligned.sortedByCoord.out.bam")
        bam_flag2    = os.path.join(args.out_dir, f"{sample}_2pass_flag2.bam")

        if paired1 and paired2:
            read_files_command = get_read_files_command(paired1)
            if os.path.isfile(log_progress) and os.path.isfile(sj_tab):
                print(f"{sample}: Paired FASTQ files have been mapped already.")
                samtools_cmd = [
                    "samtools", "view", "-b", "-f", "2",
                    "-@", str(args.threads),
                    bam_sorted,
                    "-o", bam_flag2
                ]
                run_with_conda_env(samtools_cmd, "assembly_annotation")
            else:
                star_cmd = [
                    "STAR",
                    "--runMode", "alignReads",
                    "--twopassMode", "Basic",
                    "--genomeDir", args.genomedir,
                    "--runThreadN", str(args.threads)
                ] + read_files_command + [
                    "--readFilesIn", paired1, paired2,
                    "--outSAMstrandField", "intronMotif",
                    "--outFileNamePrefix", out_prefix,
                    "--outSAMtype", "BAM", "SortedByCoordinate"
                ]
                run_with_conda_env(star_cmd, "assembly_annotation")
                samtools_cmd = [
                    "samtools", "view", "-b", "-f", "2",
                    "-@", str(args.threads),
                    bam_sorted,
                    "-o", bam_flag2
                ]
                run_with_conda_env(samtools_cmd, "assembly_annotation")
        elif single:
            read_files_command = get_read_files_command(single)
            if os.path.isfile(log_progress) and os.path.isfile(sj_tab):
                print(f"{sample}: Single FASTQ file has been mapped already.")
                samtools_cmd = [
                    "samtools", "view", "-b", "-f", "2",
                    "-@", str(args.threads),
                    bam_sorted,
                    "-o", bam_flag2
                ]
                run_with_conda_env(samtools_cmd, "annotation")
            else:
                star_cmd = [
                    "STAR",
                    "--runMode", "alignReads",
                    "--twopassMode", "Basic",
                    "--genomeDir", args.genomedir,
                    "--runThreadN", str(args.threads)
                ] + read_files_command + [
                    "--readFilesIn", single,
                    "--outSAMstrandField", "intronMotif",
                    "--outFileNamePrefix", out_prefix,
                    "--outSAMtype", "BAM", "SortedByCoordinate"
                ]
                run_with_conda_env(star_cmd, "assembly_annotation")
                samtools_cmd = [
                    "samtools", "view", "-b", "-f", "2",
                    "-@", str(args.threads),
                    bam_sorted,
                    "-o", bam_flag2
                ]
                run_with_conda_env(samtools_cmd, "assembly_annotation")
        else:
            print(f"Error: Could not find paired or single-end FASTQ files for sample {sample}.")
            sys.exit(1)

def process_isoseq_mode(args):
    """
    Process samples in iso-seq mode using minimap2 alignment.
    """
    print("Running in ISO-SEQ mode")
    
    # List files in the working directory.
    all_files = os.listdir(args.wd)
    print("Files found in working directory:", all_files)
    
    # Find iso-seq files ending with '_polished.hq.fasta'
    isoseq_samples = []
    for filename in all_files:
        if filename.lower().endswith("_polished.hq.fasta"):
            sample = extract_isoseq_basename(filename)
            if sample:
                isoseq_samples.append(sample)
    
    isoseq_samples = sorted(set(isoseq_samples))
    print("Iso-seq sample basenames found:", isoseq_samples)
    
    if not isoseq_samples:
        print("Error: No iso-seq files (*_polished.hq.fasta) found in the working directory.")
        sys.exit(1)
    
    # Process each iso-seq sample
    for sample in isoseq_samples:
        isoseq_file = find_isoseq_file(args.wd, sample)
        if not isoseq_file:
            print(f"Error: Could not find iso-seq file for sample {sample}.")
            continue
        
        sam_output = os.path.join(args.out_dir, f"{sample}.sam")
        bam_output = os.path.join(args.out_dir, f"{sample}.bam")
        
        # Check if already processed
        if os.path.isfile(bam_output):
            print(f"{sample}: Iso-seq file has been mapped already.")
            continue
        
        # Step 1: minimap2 alignment (in minimap2 conda env)
        minimap2_cmd = [
            "minimap2", "-t", str(args.threads), "-ax", "splice:hq",
            "-o", sam_output, args.genomepath, isoseq_file
        ]
        
        print(f"Running minimap2 for sample {sample}...")
        conda_cmd = ["conda", "run", "-n", "minimap2"] + minimap2_cmd
        print("Running command:", " ".join(conda_cmd))
        subprocess.run(conda_cmd, check=True)
        
        # Step 2: samtools view to convert SAM to BAM (in annotation conda env)
        samtools_cmd = [
            "samtools", "view", "--threads", str(args.threads), "-Sb",
            "-o", bam_output, sam_output
        ]
        
        print(f"Converting SAM to BAM for sample {sample}...")
        conda_cmd = ["conda", "run", "-n", "assembly_annotation"] + samtools_cmd
        print("Running command:", " ".join(conda_cmd))
        subprocess.run(conda_cmd, check=True)
        
        # Clean up SAM file to save space
        if os.path.isfile(bam_output) and os.path.getsize(bam_output) > 0:
            os.remove(sam_output)
            print(f"Removed intermediate SAM file: {sam_output}")

def main():
    parser = argparse.ArgumentParser(
        description="RNA mapping script supporting both short-read (STAR) and iso-seq (minimap2) modes. All parameters are required."
    )
    parser.add_argument("--mode", required=True, choices=["short-read", "iso-seq"],
                        help="Mapping mode: 'short-read' for STAR alignment or 'iso-seq' for minimap2 alignment")
    parser.add_argument("--genomepath", required=True,
                        help="Path to the genome FASTA file")
    parser.add_argument("--genomedir", required=False,
                        help="Directory to hold the STAR genome index (required only in short-read mode)")
    parser.add_argument("--wd", required=True,
                        help="Working directory containing RNA-seq files")
    parser.add_argument("--out_dir", required=True,
                        help="Output directory for RNA mapping")
    parser.add_argument("--threads", required=True, type=int,
                        help="Number of threads to use")
    args = parser.parse_args()

    print("Command line:", " ".join(sys.argv))
    os.makedirs(args.out_dir, exist_ok=True)
    
    if args.mode == "short-read":
        if not args.genomedir:
            print("Error: --genomedir is required for short-read mode")
            sys.exit(1)
        os.makedirs(args.genomedir, exist_ok=True)
        process_short_read_mode(args)
    elif args.mode == "iso-seq":
        process_isoseq_mode(args)

if __name__ == '__main__':
    main()