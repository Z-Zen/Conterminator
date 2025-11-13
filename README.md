
<div style="text-align: center;">
  <img src="conterminator.png">
</div>


A comprehensive Nextflow pipeline for RNA-seq quality control and contamination detection.

## Overview

**Conterminator** is a production-ready bioinformatics pipeline designed to perform comprehensive quality control and contamination detection on RNA-seq data. The pipeline combines multiple tools to identify potential contamination sources in sequencing data, including cross-strain contamination, bacterial/fungal/viral contamination, and technical artifacts.

### Key Features

- **Flexible input handling**: Direct file paths for FASTQ, pre-aligned BAM, or unmapped reads
- **Strain-specific alignment** with STAR using custom pseudogenomes or reference genome
- **Multi-database contamination screening** via BLAST
- **Microbial contamination detection** using DecontaMiner
- **Comprehensive QC suite**: FastQC, FastQ Screen, Qualimap, DeepTools, Picard, BEDTools
- **Singularity container** with all tools bundled for reproducible execution
- **SLURM cluster support** with intelligent job scheduling and resource management
- **Automatic retry** with dynamic memory scaling on out-of-memory errors
- **Flexible subsampling** for performance optimization
- **Interactive visualizations** and MultiQC reporting
- **Per-sample strain assignment** with optional default (GRCm39) for mixed experiments
- **Highly configurable** with sensible defaults

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Format](#input-format)
- [Usage](#usage)
- [Configuration](#configuration)
- [Output Structure](#output-structure)
- [Troubleshooting](#troubleshooting)
- [Execution Profiles](#execution-profiles)
- [Singularity Container](#singularity-container)
- [Contact](#contact)

## Requirements

### Software Dependencies

- **Nextflow** ≥ 23.04.0
- **Java** ≥ 11
- **Singularity** ≥ 3.5 (optional, for containerized execution)


To install Nextflow, you can run:
```bash
curl -s https://get.nextflow.io | bash
```

If you get an error that your Java is not up-to-date, read below on how to install it.

To install Java 24, you can run the following commands:
```bash
cd
wget https://download.oracle.com/java/24/archive/jdk-24.0.2_linux-x64_bin.tar.gz
tar -zxf jdk-24.0.2_linux-x64_bin.tar.gz
echo 'export JAVA_CMD="/home/${USER%%@*}/jdk-24.0.2/bin/java"' >> ~/.bashrc
echo 'export JAVA_HOME="/home/${USER%%@*}/jdk-24.0.2"' >> ~/.bashrc
echo 'export PATH="$JAVA_HOME/bin:$PATH"' >> ~/.bashrc
source .bashrc
```

**Note:** When using Singularity, Python and R dependencies are bundled in the container.

### Required Tools

The pipeline uses the following bioinformatics tools (paths configurable in `nextflow.config`):

| Tool | Purpose | Version |
|------|---------|---------|
| STAR | RNA-seq alignment | 2.7.11b |
| Samtools | BAM manipulation | 1.22.1 |
| FastQC | Quality control | 0.12.1 |
| FastQ Screen | Contamination screening | 0.16.0 |
| Qualimap | BAM quality control | 2.3 |
| Picard | GC bias analysis | 3.4.0 |
| DeepTools | GC bias computation | 3.5.6 |
| BEDTools | Coverage analysis | 2.31.1 |
| BLAST+ | Sequence alignment | 2.17.0+ |
| DecontaMiner | Microbial detection | 1.4 |
| Seqtk | FASTQ subsampling | 1.5-r133 |
| BBMap | Format conversion | 39.37 |
| MultiQC | Report aggregation | 1.31 |

### Reference Data

- Strain-specific pseudogenomes and annotations
- STAR genome indices (auto-generated if missing)
- BLAST databases for contamination screening
- FastQ Screen configuration and indices

## Installation

### Clone the Repository

```bash
git clone git@gitlab.epfl.ch:abadreddine/conterminator.git
cd conterminator
```

If the above command didn't work, make sure you have added your SSH key to gitlab. You can do it using the following commands:
```bash
ssh-keygen -t ed25519 -C "aleksandar.mihaylov@epfl.ch" -f ~/.ssh/id_ed25519_epfl
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519_epfl
cat ~/.ssh/id_ed25519_epfl.pub
```

Copy the full output of the `cat` command and paste it in the [SSH keys](https://gitlab.epfl.ch/-/user_settings/ssh_keys) of Gitlab.
Click on `Add new key` and copy the key inside the box. Then in the expiration date, click on the year and another time to prolong the key life.
Once done, click on `Add key`. Now, you can try cloning again the repository.

### Tools

#### Option 1: download pre-built tools container (Recommended)

The container includes all required bioinformatics tools pre-configured.

```bash
# Download from Google Drive
pip install gdown
gdown 1a2CevfBkMUjSt5R4AnQyDX4ofAHBZmfZ

# Test the build
singularity test conterminator.sif
```

#### Option 2: build Singularity container (requires sudo)

Build the Singularity container with all dependencies bundled:

```bash
# Build the container
sudo singularity build conterminator.sif conterminator.def

# Test the build
singularity test conterminator.sif
```

#### Option 3: native installation 

Install all required tools manually and configure paths in `nextflow.config`.

Edit `nextflow.config` to specify paths to installed tools and reference data:

```groovy
params {
    // Tool paths
    star_bin = "/path/to/STAR"
    samtools_bin = "/path/to/samtools"
    fastqc_bin = "/path/to/fastqc"
    // ... (see nextflow.config for all options)
    
    // Reference directories
    strains_base_dir = "/path/to/pseudogenomes"          // Strain-specific pseudogenomes
    standard_references_dir = "/path/to/references/Mus"  // Standard reference genomes (GRCm39, etc.)
    star_index_dir = "/path/to/star/indices"
    contamination_blast_dbs = "/path/to/blast/databases"
}
```

#### Test Installation

```bash
# Native installation
nextflow run main.nf --help

# With Singularity
nextflow run main.nf --singularity_path /path/to/conterminator.sif -profile singularity --help

# With Singularity + Slurm
nextflow run main.nf --singularity_path /path/to/conterminator.sif -profile singularity,slurm --help
```

## Quick Start

### Step 1: Prepare Sample Sheet

Create a tab-separated file (e.g., `samples.tsv`):

```tsv
sample	strain	type	read1	read2
sample1	C57BL_6J	fastq	/data/sample1_R1.fq.gz	/data/sample1_R2.fq.gz
sample2	C57BL_6J	fastq	/data/sample2_R1.fq.gz	/data/sample2_R2.fq.gz
```

### Step 2: Run Pipeline

#### Minimal Example (Local)

```bash
nextflow run main.nf \
    --sample_sheet samples.tsv \
    --outdir results \
    -bg &> results.log
```

#### With Singularity Container

```bash
nextflow run main.nf \
    -profile singularity \
    --singularity_path /path/to/conterminator.sif \
    --sample_sheet samples.tsv \
    --outdir results \
    -bg &> results.log
```

#### On SLURM Cluster with Singularity (Recommended)

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --singularity_path /path/to/conterminator.sif \
    --sample_sheet samples.tsv \
    --outdir results \
    --max_parallel_samples 10 \
    -bg &> results.log
```

#### Complete Example with Performance Tuning

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --singularity_path /path/to/conterminator.sif \
    --sample_sheet samples.tsv \
    --outdir results_full \
    --subset_for_fastq_qc true \
    --subset_fastq_qc_reads 100000 \
    --subset_bam_for_qc true \
    --bam_qc_subset_mapped 200000 \
    --max_parallel_samples 10 \
    -bg &> results.log
```

## Input Format

### Sample Sheet (Required)

The pipeline requires a **tab-separated** sample sheet file specifying samples, input types, and file paths.

#### Format Specification

Create a TSV file (e.g., `samples.tsv`) with **5 columns**:

```tsv
sample	strain	type	read1	read2
sample1	C57BL_6J	fastq	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
sample2	C57BL_6J	bam	/path/to/sample2.bam
sample3	DBA_2J	unmapped_fastq	/path/to/sample3_unmapped_R1.fastq.gz	/path/to/sample3_unmapped_R2.fastq.gz
sample4		fastq	/path/to/sample4_R1.fastq.gz	/path/to/sample4_R2.fastq.gz
sample5		bam	/path/to/sample5.bam
```

#### Column Descriptions

| Column | Required | Description | Values/Examples |
|--------|----------|-------------|-----------------|
| **sample** | Yes | Unique sample identifier | `sample1`, `exp_001`, `ctrl-A` |
| **strain** | No | Reference strain name (defaults to `GRCm39` if empty) | `C57BL_6J`, `DBA_2J`, `BALB_cJ`, or leave empty |
| **type** | Yes | Input data type | `fastq`, `bam`, or `unmapped_fastq` |
| **read1** | Yes | Path to first read file or BAM file | `/path/to/sample_R1.fastq.gz` or `/path/to/sample.bam` |
| **read2** | Conditional | Path to second read file (required for `fastq` and `unmapped_fastq`, empty for `bam`) | `/path/to/sample_R2.fastq.gz` or empty |

#### Input Type Definitions

| Type | Description | Use Case | Required Columns |
|------|-------------|----------|------------------|
| **fastq** | Paired-end FASTQ files for alignment | Raw sequencing data to be aligned with STAR | `read1`, `read2` |
| **bam** | Pre-aligned BAM file | Already aligned data, skip STAR alignment | `read1` only |
| **unmapped_fastq** | Unmapped reads in FASTQ format | Already extracted unmapped reads for contamination analysis | `read1`, `read2` |

#### Requirements and Rules

1. **Header row is mandatory** - First line must be: `sample	strain	type	read1	read2`
2. **Tab-separated format** - Columns must be separated by tabs (not spaces)
3. **Five columns required** - All 5 columns must be present in header
4. **Unique sample IDs** - Each sample name must be unique across the sheet
5. **Valid type values** - Type must be one of: `fastq`, `bam`, or `unmapped_fastq`
6. **Strain handling**:
   - If strain is empty or not specified, defaults to `GRCm39`
   - If strain is specified, it must exist in `params.strains_base_dir`
7. **read2 rules**:
   - For `type=fastq`: read2 is **required** (paired-end reads)
   - For `type=unmapped_fastq`: read2 is **required** (paired-end reads)
   - For `type=bam`: read2 must be **empty** (BAM files don't have separate read files)
8. **File paths** - All file paths in `read1` and `read2` must exist and be accessible

#### Comprehensive Examples

**Example 1: Mixed input types**
```tsv
sample	strain	type	read1	read2
WT_rep1	C57BL_6J	fastq	/data/exp1/WT_rep1_R1.fq.gz	/data/exp1/WT_rep1_R2.fq.gz
WT_rep2	C57BL_6J	fastq	/data/exp1/WT_rep2_R1.fq.gz	/data/exp1/WT_rep2_R2.fq.gz
KO_rep1	DBA_2J	bam	/data/exp1/KO_rep1_aligned.bam
archived_sample		bam	/archive/old_alignment.bam
contaminated_sample	BALB_cJ	unmapped_fastq	/data/unmapped/sample_R1.fq.gz	/data/unmapped/sample_R2.fq.gz
```

**Example 2: All FASTQ inputs (standard RNA-seq workflow)**
```tsv
sample	strain	type	read1	read2
ctrl_1	C57BL_6J	fastq	/mnt/data/ctrl_1_R1.fastq.gz	/mnt/data/ctrl_1_R2.fastq.gz
ctrl_2	C57BL_6J	fastq	/mnt/data/ctrl_2_R1.fastq.gz	/mnt/data/ctrl_2_R2.fastq.gz
treat_1	C57BL_6J	fastq	/mnt/data/treat_1_R1.fastq.gz	/mnt/data/treat_1_R2.fastq.gz
treat_2	C57BL_6J	fastq	/mnt/data/treat_2_R1.fastq.gz	/mnt/data/treat_2_R2.fastq.gz
```

**Example 3: Pre-aligned BAM files (QC only)**
```tsv
sample	strain	type	read1	read2
sample1	C57BL_6J	bam	/alignments/sample1.bam
sample2	DBA_2J	bam	/alignments/sample2.bam
sample3		bam	/alignments/sample3.bam
```

**Example 4: Using default strain (GRCm39)**
```tsv
sample	strain	type	read1	read2
exp001		fastq	/data/exp001_R1.fq.gz	/data/exp001_R2.fq.gz
exp002		fastq	/data/exp002_R1.fq.gz	/data/exp002_R2.fq.gz
exp003		bam	/data/exp003.bam
```

#### Strain Reference Requirements

For each strain specified in the sample sheet (or the default `GRCm39`), the following files must exist in `params.strains_base_dir`:

```
strains_base_dir/
└── STRAIN_NAME/
    ├── *pseudogenome__strain_STRAIN_NAME.fa.gz    # FASTA reference
    └── *pseudogenome__strain_STRAIN_NAME.gtf.gz   # Gene annotations
```

**Example for C57BL_6J strain:**
```
/path/to/strains_base_dir/C57BL_6J/
├── mm39_pseudogenome__strain_C57BL_6J.fa.gz
└── mm39_pseudogenome__strain_C57BL_6J.gtf.gz
```

**Example for default GRCm39 strain:**
```
/path/to/strains_base_dir/GRCm39/
├── GRCm39.genome.fa.gz
└── gencode.vM35.primary_assembly.annotation.gtf.gz
```

**Note on Standard Reference Genomes (v1.2+):** The pipeline now supports storing standard reference genomes (like GRCm39, GRCm38) in a separate directory specified by `params.standard_references_dir` (default: `/mnt/sas/Data/References/Mus`). The pipeline will automatically search both `strains_base_dir` (for strain-specific pseudogenomes) and `standard_references_dir` (for standard references) when resolving strain names. This allows you to maintain your existing HDP pseudogenome collection while using standard Ensembl/GENCODE references without reorganizing your directory structure.

## Usage

### Key Parameters

#### Input/Output
```bash
--sample_sheet <file>      # Sample sheet with 5 columns: sample, strain, type, read1, read2 (required)
--outdir <path>            # Output directory (default: results)
```

#### Subsampling (Performance Tuning)
```bash
--subset_for_fastq_qc true          # Subsample for FastQ QC (default: true)
--subset_fastq_qc_reads 100000      # Reads for QC (default: 100k)

--subset_for_star false             # Subsample for STAR (default: false)
--subset_star_reads 100000          # Reads for STAR if enabled

--subset_bam_for_qc true            # Subsample BAMs for QC (default: true)
--bam_qc_subset_mapped 200000       # Mapped reads for BAM QC (default: 200k)

--subset_unmapped_for_blast true    # Subsample for BLAST (default: true)
--unmapped_subset_reads 100000      # Unmapped reads for BLAST (default: 100k)
```

#### Tool Toggles
```bash
--run_star_alignment true     # Run STAR alignment (default: true)
--run_fastqc true             # Run FastQC (default: true)
--run_fastq_screen true       # Run FastQ Screen (default: true)
--run_deeptools true          # Run DeepTools GC bias (default: true)
--run_picard_gc true          # Run Picard GC bias (default: true)
--run_bedtools_gc true        # Run BEDTools coverage (default: true)
--run_qualimap true           # Run Qualimap (default: true)
--run_mapinsights true        # Run MapInsights (default: true)
--run_decontaminer true       # Run DecontaMiner (default: true)
--run_contamination_check true # Run BLAST contamination (default: true)
--run_multiqc true            # Run MultiQC (default: true)
```

#### Resources
```bash
--max_cpus 16              # Maximum CPUs per process
--max_mem "32 GB"          # Maximum memory per process
--max_time "24h"           # Maximum time per process
--max_parallel_samples 4   # Max samples in parallel
```

## Configuration

### Qualimap Settings

```bash
--qualimap_mode "both"                    # "bamqc", "rnaseq", or "both"
--qualimap_genome "mm10"                  # Reference genome name
--qualimap_protocol "strand-specific-reverse"
--qualimap_threads 8
```

### Contamination Detection

```bash
--contamination_blast_dbs "/path/to/blast/dbs"  # BLAST database directory
--contamination_evalue "1e-10"                   # E-value threshold
--blast_max_parallel 10                          # Max parallel BLAST jobs
```

### DecontaMiner Settings

```bash
--decontaminer_config "/path/to/config.txt"
--decontaminer_organisms "bfv"          # b=bacteria, f=fungi, v=viruses
--decontaminer_pairing "P"              # P=paired, S=single
--decontaminer_quality_filter "yes"
--decontaminer_ribo_filter "yes"
```

## Output Structure

```
results/
├── Input/
│   ├── reference/              # Prepared references
│   ├── strain_references/      # Strain-specific files
│   ├── subsampled_bams/        # Subsampled BAM files
│   ├── subsampled_fastq/       # Subsampled FASTQ files
│   └── sample_sheet.tsv        # Copy of input sample sheet
├── Output/
│   ├── bedtools_gc/            # GC content coverage
│   ├── contamination_check/    # BLAST results and plots
│   ├── decontaminer/           # DecontaMiner reports
│   ├── deeptools_gc_bias/      # DeepTools GC bias
│   ├── fastqc/                 # FastQC reports
│   ├── fastq_screen/           # FastQ Screen results
│   ├── mapinsights/            # MapInsights reports
│   ├── multiqc/                # MultiQC aggregate report
│   ├── picard_gc_bias/         # Picard GC metrics
│   └── qualimap/               # Qualimap QC
├── Temporary/
│   ├── decontaminer/           # Intermediate files
│   └── star_alignment/         # STAR outputs
├── NextflowReports/            # Internal reports by Nextflow
├── pid.txt                     # Contains Nextflow PIDs
└── pipeline_info.txt           # Run metadata
```

## Troubleshooting

### Common Issues

**Issue: STAR index not found**
```bash
# The pipeline will auto-build missing indices
# Ensure strains_base_dir contains:
# strains_base_dir/STRAIN_NAME/*pseudogenome__strain_STRAIN_NAME.fa.gz
# strains_base_dir/STRAIN_NAME/*pseudogenome__strain_STRAIN_NAME.gtf.gz
```

**Issue: Out of memory errors**
```bash
# Increase memory allocation
--max_mem "64 GB"

# Enable more aggressive subsampling
--subset_bam_for_qc true
--bam_qc_subset_mapped 100000
```

**Issue: Pipeline hangs during BLAST**
```bash
# Reduce parallel BLAST jobs
--blast_max_parallel 5

# Enable BLAST subsampling
--subset_unmapped_for_blast true
--unmapped_subset_reads 50000
```

### Nextflow tips

#### Resuming failed runs

If your run has failed due to an error, you can resume it after fixing the problem by running the following:

```bash
# Get all the runs launched by nextflow
nextflow log

# Copy the session id of the run you want to resume. It looks like this: 2b29d621-8eff-400c-83a2-05146c4f6131

# Resume your pipeline by using exactly the same command as before but by adding -resume [session id]
# For example, if your command was
nextflow run main.nf --outdir myresults --sample_sheet samples.tsv -profile singularity,slurm
# It will become
nextflow run main.nf --outdir myresults --sample_sheet samples.tsv -profile singularity,slurm -resume 2b29d621-8eff-400c-83a2-05146c4f6131
```

#### Run Nextflow in the background

You can simply run nextflow in the background using `-bg`:

```bash
nextflow run main.nf --outdir myresults --sample_sheet samples.tsv -profile singularity,slurm -bg &> myrun.log
```

#### Stopping background run

You can find the latest running pid in the `pid.txt` in the root of output directory.

```bash
cat pid.txt
# ---------------------
# 4149804
# nextflow run main.nf --sample_sheet input_test/sample_sheet.tsv --outdir results_test --singularity_path /mnt/sas/Users/abadreddine/Projects/Conterminator/conterminator.sif -profile singularity -resume 2b29d621-8eff-400c-83a2-05146c4f6131
# Start: 07-Nov-2025 01:25:30

# Then you can use the PID number to kill the running job
kill 4149804
```

#### Cleaning the work directory

Nextflow works by running his processes in a folder called `work` in your current path. The `work` directory is important to resume your jobs. If this folder gets deleted, you cannot resume your jobs anymore.

However, once your analysis has finished, you can clean your working directory by running:

```bash
rm -rf work .nextflow.log*
# OR 
nextflow clean -f $(nextflow log -q)
```

## Contact

**Author:** Alaa Badreddine

**Project Repository:** https://gitlab.epfl.ch/abadreddine/conterminator

For bug reports and feature requests, please open an issue on GitLab.

---

## Version History

- **v1.2** (December 2025) - Enhanced reference genome support:
  - Dual directory support for reference genomes:
    - `strains_base_dir` for strain-specific pseudogenomes
    - `standard_references_dir` for standard reference genomes (GRCm39, GRCm38, etc.)
  - Automatic directory resolution - pipeline searches both locations for each strain
  - Flexible reference file pattern matching for standard genomes
  - Fixed WRITE_PIPELINE_INFO staging to avoid file collision errors
  - Pipeline configuration files (main.nf, nextflow.config, conterminator.def) now copied to output for reproducibility
  - Dynamic user path support using `${System.getProperty('user.home')}`

- **v1.1** (November 2025) - Major updates:
  - Full Singularity container support with all tools bundled
  - SLURM cluster execution with job scheduling
  - Per-process resource configuration with automatic retry on OOM
  - Dynamic memory scaling (2×/3× on retry)
  - Parallelization control with `maxForks` parameter
  - Enhanced error handling with intelligent retry strategy
  - Singularity information in pipeline reports
  - Improved channel scoping for complex workflows
  - New flexible sample sheet format supporting multiple input types:
  - Direct FASTQ file paths (no directory scanning required)
  - Pre-aligned BAM files for QC-only workflows
  - Unmapped FASTQ files for contamination-only analysis
  - Optional strain specification with GRCm39 default

- **v1.0** (2025) - Initial release with full QC and contamination detection suite
