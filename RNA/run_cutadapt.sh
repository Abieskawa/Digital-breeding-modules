#!/bin/bash
# run_cutadapt.sh
#
# This script uses cutadapt to process paired‐end FASTQ files.
#
# It supports front‐trimming for each sample via the -u (read1) and -U (read2)
# options. The front‐trim values can be supplied either as a single number (applied
# globally) or via a file (specified with -f) that contains 2 columns:
#
#   2 columns: <sample_filename_key> <trimValue>
#       => For example, if your FASTQ files are named SRR26245751_1.fastq and
#          SRR26245751_2.fastq, the file should contain one line for each read:
#
#          SRR26245751_1 10
#          SRR26245751_2 10
#
#       => Alternatively, if using _R1/_R2 naming convention:
#          4202P-HPGonad_R1 10
#          4202P-HPGonad_R2 10
#
# Adapter sequences can be specified via -a1 and -a2. If not provided, the script will
# look for a file named "adapter_list.txt" in the same directory as this script; if that
# file isn't found, default adapters are used:
#
#   For read1 (-a):
#     AGATCGGAAGAG
#     AAAAAAAAAAAA
#     GGGGGGGGGGGG
#
#   For read2 (-A):
#     AGATCGGAAGAG
#     AAAAAAAAAAAA
#     GGGGGGGGGGGG
#
# Supported file extensions for paired‐end inputs are: .fastq, .fq, .fastq.gz, and .fq.gz
#
# Example (global front trimming):
#   nohup bash run_cutadapt.sh -d /path/to/05_raw_RNA -o /path/to/05_cleaned_RNA \
#         -t 128 -l 5 -q 20 -Q 20 -f 10 -a1 YOUR_ADAPTER1 -a2 YOUR_ADAPTER2 > run_cutadapt.log 2>&1 &
#
# Example (per-sample trimming using a trim file):
#   nohup bash run_cutadapt.sh -d /path/to/05_raw_RNA -o /path/to/05_cleaned_RNA \
#         -t 128 -l 5 -q 20 -Q 20 -f sample_trim.txt > run_cutadapt.log 2>&1 &

############################################
# Defaults for optional parameters
############################################
threads=16
min_length=5         # Minimum read length (-m option)
quality1=20          # Quality cutoff for read1 (-q option)
quality2=20          # Quality cutoff for read2 (-Q option)
trim_front=""        # Either a single number (global) or a filename with per-read trim values

# Adapter parameters. If provided, these override the defaults.
adapter1=""
adapter2=""

# Required parameters
wd=""
output_dir=""

############################################
# Function to display usage information
############################################
usage() {
  cat <<EOF
Usage: $0 -d <working_directory> -o <output_directory> -t <threads> -l <min_length> \\
          -q <quality_R1> -Q <quality_R2> [-f <trim_front_value_or_file>] \\
          [-a1 <adapter_read1>] [-a2 <adapter_read2>]

This script expects paired‐end files named as: <basename>_1.<ext> and <basename>_2.<ext>,
where <ext> can be one of: fastq, fq, fastq.gz, or fq.gz.

Example (global front trimming):
  $0 -d /path/to/05_raw_RNA -o /path/to/05_cleaned_RNA -t 16 -l 5 -q 20 -Q 20 -f 10 -a1 YOUR_ADAPTER1 -a2 YOUR_ADAPTER2

Example (per-sample trimming):
  $0 -d /path/to/05_raw_RNA -o /path/to/05_cleaned_RNA -t 16 -l 5 -q 20 -Q 20 -f sample_trim.txt

If you supply a trim file to -f, it must have 2 columns per line:
  <sample_filename_key> <trimValue>
For paired‐end data, the keys must match the actual filenames:
  e.g.,
    SRR26245751_1 10
    SRR26245751_2 10
  or using _R1/_R2 convention:
    4202P-HPGonad_R1 10
    4202P-HPGonad_R2 10

Options:
  -d      Working directory (required)
  -o      Output directory (required; must already exist)
  -t      Number of threads (default: $threads)
  -l      Minimum read length (default: $min_length)
  -q      Quality cutoff for read1 (default: $quality1)
  -Q      Quality cutoff for read2 (default: $quality2)
  -f      (Optional) Either a single number or a filename with per-read trimming values
  -a1     (Optional) Adapter sequence for read1
  -a2     (Optional) Adapter sequence for read2
  -? or --help   Display this help message
EOF
  exit 1
}

############################################
# If no arguments provided, show usage
############################################
if [ $# -eq 0 ]; then
  usage
fi

############################################
# Save the original command line for logging
############################################
orig_cmd="$0 $@"
echo "Script invoked with: $orig_cmd"

############################################
# Parse command-line options
############################################
while [ $# -gt 0 ]; do
  case "$1" in
    -d)
      wd="$2"
      shift 2
      ;;
    -o)
      output_dir="$2"
      shift 2
      ;;
    -t)
      threads="$2"
      shift 2
      ;;
    -l)
      min_length="$2"
      shift 2
      ;;
    -q)
      quality1="$2"
      shift 2
      ;;
    -Q)
      quality2="$2"
      shift 2
      ;;
    -f)
      trim_front="$2"
      shift 2
      ;;
    -a1)
      adapter1="$2"
      shift 2
      ;;
    -a2)
      adapter2="$2"
      shift 2
      ;;
    -\? | --help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

############################################
# Validate required parameters
############################################
if [[ -z "$wd" ]]; then
  echo "Error: Working directory (-d) not specified."
  usage
fi
if [[ -z "$output_dir" ]]; then
  echo "Error: Output directory (-o) not specified."
  usage
fi

############################################
# Determine script directory before changing directories
############################################
# Handle both absolute and relative paths robustly
if command -v realpath >/dev/null 2>&1 && realpath "$0" >/dev/null 2>&1; then
  script_dir=$(cd "$(dirname "$(realpath "$0")")" && pwd)
else
  # Fallback method if realpath fails or is not available
  script_dir=$(cd "$(dirname "$0")" && pwd)
fi

############################################
# Convert directories to absolute paths
############################################
initial_dir=$(pwd)
if [[ ! -d "$wd" ]]; then
  echo "Error: Working directory '$wd' does not exist!"
  exit 1
fi
wd=$(cd "$wd" && pwd)

# Make output_dir absolute if it isn't already
if [[ "${output_dir:0:1}" != "/" ]]; then
  output_dir="$initial_dir/$output_dir"
fi

############################################
# Process the -f parameter (front trim)
############################################
# If -f is provided, it can be a file (with 2 columns) or a global number.
declare global_trim_R1=""
declare global_trim_R2=""
declare -A per_sample_trim_R1
declare -A per_sample_trim_R2

if [[ -n "$trim_front" ]]; then
  if [[ -f "$trim_front" ]]; then
    # Parse the file line by line.
    while IFS= read -r line || [[ -n "$line" ]]; do
      # Remove carriage return if present (handles DOS/Windows line endings)
      line=$(echo "$line" | tr -d '\r')
      [[ -z "$line" || "$line" =~ ^# ]] && continue
      fields=( $line )
      if [[ ${#fields[@]} -eq 2 ]]; then
        sample="${fields[0]}"
        trim_value="${fields[1]}"
        # Support both _1/_2 and _R1/_R2 naming conventions
        if [[ "$sample" =~ _1$ ]] || [[ "$sample" =~ _R1$ ]]; then
          per_sample_trim_R1["$sample"]="$trim_value"
        elif [[ "$sample" =~ _2$ ]] || [[ "$sample" =~ _R2$ ]]; then
          per_sample_trim_R2["$sample"]="$trim_value"
        else
          echo "Warning: Trim key '$sample' does not match expected pattern (_1/_2 or _R1/_R2); skipping."
        fi
      else
        echo "Warning: Unexpected number of columns in trim file line: $line"
      fi
    done < "$trim_front"
  else
    # Otherwise, treat -f as a global trim value.
    global_trim_R1="$trim_front"
    global_trim_R2="$trim_front"
  fi
fi

############################################
# Change to the working directory
############################################
cd "$wd" || { echo "Error: Failed to change directory to '$wd'"; exit 1; }

# Enable nullglob so that unmatched globs expand to null
shopt -s nullglob

############################################
# Process each sample
############################################
# This script expects paired‐end files named as: <basename>_1.<ext> and <basename>_2.<ext>,
# where <ext> can be fastq, fq, fastq.gz, or fq.gz.
found_samples=0
for file in *_[1].fastq *_[1].fq *_[1].fastq.gz *_[1].fq.gz; do
  if [ ! -f "$file" ]; then
    continue
  fi
  found_samples=1
  
  # Extract basename by removing the trailing _1 and extension.
  basename="${file%_1.*}"
  ext="${file##*_1.}"
  
  read1="${basename}_1.${ext}"
  read2="${basename}_2.${ext}"
  
  if [[ ! -f "$read1" || ! -f "$read2" ]]; then
    echo "Skipping sample '$basename': one or both paired-end files are missing."
    continue
  fi
  
  echo "Processing sample: $basename"
  
  ############################################
  # Determine front trim values for each read.
  # Support multiple naming conventions for lookup keys
  ############################################
  sample_key_r1="${basename}_1"
  sample_key_r2="${basename}_2"
  sample_key_r1_alt="${basename}_R1"
  sample_key_r2_alt="${basename}_R2"
  
  sample_trim_r1=""
  sample_trim_r2=""
  
  # Try to find trim values using different key formats
  if [[ -n "${per_sample_trim_R1[$sample_key_r1]}" ]]; then
    sample_trim_r1="${per_sample_trim_R1[$sample_key_r1]}"
  elif [[ -n "${per_sample_trim_R1[$sample_key_r1_alt]}" ]]; then
    sample_trim_r1="${per_sample_trim_R1[$sample_key_r1_alt]}"
  elif [[ -n "$global_trim_R1" ]]; then
    sample_trim_r1="$global_trim_R1"
  fi
  
  if [[ -n "${per_sample_trim_R2[$sample_key_r2]}" ]]; then
    sample_trim_r2="${per_sample_trim_R2[$sample_key_r2]}"
  elif [[ -n "${per_sample_trim_R2[$sample_key_r2_alt]}" ]]; then
    sample_trim_r2="${per_sample_trim_R2[$sample_key_r2_alt]}"
  elif [[ -n "$global_trim_R2" ]]; then
    sample_trim_r2="$global_trim_R2"
  fi
  
  ############################################
  # Define adapter options for cutadapt
  ############################################
  adapter_fasta="$script_dir/adapter_list.txt"
  
  if [[ -n "$adapter1" ]]; then
    adapter_opts1=( -a "$adapter1" )
  elif [[ -f "$adapter_fasta" ]]; then
    adapter_opts1=( -a "file:$adapter_fasta" )
  else
    adapter_opts1=( -a "AGATCGGAAGAG" -a "AAAAAAAAAAAA" -a "GGGGGGGGGGGG" )
  fi
  
  if [[ -n "$adapter2" ]]; then
    adapter_opts2=( -A "$adapter2" )
  elif [[ -f "$adapter_fasta" ]]; then
    adapter_opts2=( -A "file:$adapter_fasta" )
  else
    adapter_opts2=( -A "AGATCGGAAGAG" -A "AAAAAAAAAAAA" -A "GGGGGGGGGGGG" )
  fi
  
  ############################################
  # Build the cutadapt command
  ############################################
  cmd=(cutadapt)
  cmd+=( "${adapter_opts1[@]}" "${adapter_opts2[@]}" )
  cmd+=( -j "$threads" -q "$quality1" -Q "$quality2" -m "$min_length" )
  
  # Add front trim options if trim values were determined.
  if [[ -n "$sample_trim_r1" ]]; then
    cmd+=( -u "$sample_trim_r1" )
  fi
  if [[ -n "$sample_trim_r2" ]]; then
    cmd+=( -U "$sample_trim_r2" )
  fi
  
  # Define output filenames (placed in the output directory)
  out1="${output_dir}/${basename}_1.cleaned.${ext}"
  out2="${output_dir}/${basename}_2.cleaned.${ext}"
  
  cmd+=( -o "$out1" -p "$out2" "$read1" "$read2" )
  
  echo "Running command: ${cmd[*]}"
  "${cmd[@]}"
done

if [ $found_samples -eq 0 ]; then
  echo "No input files matching supported extensions were found in $(pwd)"
fi

echo "All samples processed at $(date)"