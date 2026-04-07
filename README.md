# Digital-breeding-modules
This repository will be deposited with developing digital-breeding related tools. The datasets used to construct model will skipped here due to potential authority issues.

## Output Files and Directories Overview
The pipeline writes everything under `output_dir` (default: `DGBreeding`).

- `00_Preprocessed_DNA/`: fastp-cleaned FASTQs (`*_R1.cleaned.fastq.gz`, `*_R2.cleaned.fastq.gz`) and `fastp_reports/` HTML/JSON.
- `01_DNAseq_alignment/`: copied reference FASTA + BWA index + `.fai`, alignment logs, BAMs (`*.bam`/`*.bai`), duplicate stats, optional `bam_archive.tar.gz`.
- `02_Variant_Calling/`: DeepVariant per-sample `*.vcf.gz`/`*.g.vcf.gz`, `gvcf.list`, `cohort.merged.vcf.gz`, `cohort.filtered.vcf.gz`, `cohort.filtered.ldpruned.vcf.gz`, logs.
- `evaluation/`: FastQC/MultiQC reports, Kraken2 summaries, mapping-yield CSV/plots, variant stats, circos plots.
- `<prediction_output_dir>/`: prediction outputs (relative to `output_dir`; default in config: `prediction/`), including:
  - `evaluation/PCA_Scree_Plot/` (PCA plots/scree)
  - `GWAS_dir/` (per-fold inputs/outputs, BLINK logs, Manhattan/QQ plots)
  - `PRS_dir/` (per-fold PRS train/test scores and train-derived weights)
  - `Model_result/` (metrics CSVs, probabilities/test labels `.npy`)
  - `Shap_dir/` (SHAP plots/feature importance CSVs)
  - `AUC_dir/` (ROC/AUC comparison plots + `auc_summary.csv`)

## Conda Environments (Biotools + ML)
The Docker image builds two separate conda environments using micromamba:

- `DGbreeding-biotools`: alignment/variant calling/QC (fastp, bwa-mem2, samtools, bcftools, glnexus, plink2, etc.)
- `DGbreeding-ml`: model training/inference (numpy, pandas, scikit-learn, xgboost, shap, etc.)

RNA-seq note:
- `fastp --dedup` is DNA-only in this project. Do not use duplicate removal when preprocessing RNA-seq libraries; use `--lib-type RNA` without `--dedup`.

To construct them locally, copy the YAML specs from `Dockerfile` and run:

```bash
micromamba create -n DGbreeding-biotools -f biotools.yml
micromamba create -n DGbreeding-ml -f ml.yml
```

If you use conda instead of micromamba, replace `micromamba` with `conda`.


## Illustration of concept below steps 
Run the pipeline via:
`python run_prediction_pipeline.py --config_file path/to/configure.txt`

Core steps (`step=1-6` in config):

    raw FASTQ + ref FASTA (+ optional GFF)
        |
        | [1] fastp clean + reference index (bwa-mem2, samtools faidx)
        v
    00_Preprocessed_DNA/        01_DNAseq_alignment/
        |
        | [2] bwa-mem2 alignment -> BAM (+ sort/markdup/filter)
        v
    01_DNAseq_alignment/*.bam + *.bai
        |
        | [3] DeepVariant -> per-sample VCF/gVCF
        |     GLnexus -> cohort.merged.vcf.gz
        |     plink2 QC -> cohort.filtered.vcf.gz
        |     plink2 LD-prune -> cohort.filtered.ldpruned.vcf.gz
        v
    02_Variant_Calling/
        |
        | [4] BLINK GWAS + PCA (uses phenotype CSV + VCF + chromosome mapping)
        v
    <prediction_output_dir>/GWAS_dir + evaluation/PCA_Scree_Plot
        |
        | [5] optional PRS scoring (fold-specific train/test scores; train-derived weights only)
        | [6] prediction models + SHAP + ROC/AUC
        v
    <prediction_output_dir>/PRS_dir + Model_result + Shap_dir + AUC_dir

Evaluation steps (`eva_step=0-6` in config, can be run alongside or after core steps):

- `0`: FastQC + MultiQC on raw FASTQ -> `evaluation/fastqc/raw`, `evaluation/multiqc/raw`
- `1`: FastQC + MultiQC on processed FASTQ + Kraken2 -> `evaluation/fastqc/processed`, `evaluation/kraken`
- `2`: Mapping yield summary -> `evaluation/unique_reads_in_bam_vs_fastp.csv`, `evaluation/mapping_yield_boxplots.png`
- `3`: Variant stats + Circos -> `evaluation/variant_stats/*.stats`, `evaluation/variant_summary.csv`, `evaluation/variant_circos/*.png`
- `4`: Manhattan/QQ plots -> `<prediction_output_dir>/GWAS_dir/*/manhattan_plot_*.png`, `qq_plot_*.png`
- `5`: PRS evaluation (not implemented)
- `6`: ROC/AUC plots -> `<prediction_output_dir>/AUC_dir/auc_comparison_*.png`, `auc_summary.csv`

## Minimal Config Additions
New prediction feature modes keep the old SNP workflow as the default:

```ini
prediction_feature_mode=baseline_snp
```

- `baseline_snp`: unchanged current behavior.
- `prs_only`: step 6 uses only PRS columns from step 5.
- `snp_plus_prs`: step 6 appends PRS columns onto the existing SNP matrix.

Minimal DeepVariant backend/resource keys:

```ini
deepvariant_mode=auto
deepvariant_docker_mode=true
deepvariant_gpu_image=google/deepvariant:1.10.0-gpu
deepvariant_cpu_image=google/deepvariant:1.10.0-beta
gpu_devices=0
step3_max_concurrent_samples=1
step3_gpu_jobs=1
step3_cpu_jobs=1
deepvariant_num_shards=96
```

- Only DeepVariant uses GPU.
- GLnexus, plink2, BLINK, PRS, model training, and the rest of the pipeline remain CPU-based.

## Example Snippets
Baseline SNP-only, unchanged default:

```ini
prediction_feature_mode=baseline_snp
step=4,6
```

PRS-only experiment:

```ini
prediction_feature_mode=prs_only
step=4,5,6
prs_top_n_snps_list=1-10,11-20,21-50
```

SNP + PRS experiment:

```ini
prediction_feature_mode=snp_plus_prs
step=4,5,6
top_n_snps_list=10,50
prs_top_n_snps_list=1-10,11-20,21-50
```

DeepVariant GPU via Docker:

```ini
deepvariant_mode=gpu
deepvariant_docker_mode=true
gpu_devices=0
step3_gpu_jobs=1
deepvariant_num_shards=96
```

Auto mode with clean CPU fallback:

```ini
deepvariant_mode=auto
deepvariant_docker_mode=true
deepvariant_cpu_image=google/deepvariant:1.10.0-beta
step3_cpu_jobs=1
```


## Tools, Links, and Papers
Papers are listed when available; some tools only have preprints or software docs.

| Tool | Link | Paper |
| --- | --- | --- |
| fastp | https://github.com/OpenGene/fastp | Chen et al., 2018, Bioinformatics |
| bwa-mem2 | https://github.com/bwa-mem2/bwa-mem2 | Vasimuddin et al., 2019, bioRxiv (preprint) |
| SAMtools | https://github.com/samtools/samtools | Li et al., 2009, Bioinformatics |
| BCFtools | https://github.com/samtools/bcftools | Li et al., 2009, Bioinformatics |
| DeepVariant | https://github.com/google/deepvariant | Poplin et al., 2018, Nature Biotechnology |
| GLnexus | https://github.com/dnanexus-rnd/GLnexus | Yun et al., 2020, bioRxiv (preprint) |
| PLINK2 | https://www.cog-genomics.org/plink/ | Purcell et al., 2007, AJHG; Chang et al., 2015, GigaScience |
| BLINK (GWAS) | https://zzlab.net/ | Huang et al., 2019, Bioinformatics |
| FastQC | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ | N/A |
| MultiQC | https://github.com/ewels/MultiQC | Ewels et al., 2016, Bioinformatics |
| Kraken2 | https://github.com/DerrickWood/kraken2 | Wood et al., 2019, Genome Biology |
| STAR (RNA-seq) | https://github.com/alexdobin/STAR | Dobin et al., 2013, Bioinformatics |
| StringTie (RNA-seq) | https://github.com/gpertea/stringtie | Pertea et al., 2015, Nature Biotechnology |
| scikit-learn | https://github.com/scikit-learn/scikit-learn | Pedregosa et al., 2011, JMLR |
| XGBoost | https://github.com/dmlc/xgboost | Chen and Guestrin, 2016, KDD |
| SHAP | https://github.com/shap/shap | Lundberg and Lee, 2017, NeurIPS |

If you use a different distribution or repo for any tool (for example, a BLINK GitHub fork), swap in that link.
