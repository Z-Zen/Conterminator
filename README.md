
<div style="text-align: center;">
  <img src="conterminator.png">
</div>


A comprehensive Nextflow pipeline for RNA-seq quality control and contamination detection.

## Overview

**Conterminator** is a production-ready bioinformatics pipeline designed to perform comprehensive quality control and contamination detection on RNA-seq data. The pipeline combines multiple tools to identify potential contamination sources in sequencing data, including cross-strain contamination, bacterial/fungal/viral contamination, and technical artifacts.

### Key Features

- **Strain-specific alignment** with STAR using custom pseudogenomes or reference genome
- **Multi-database contamination screening** via BLAST
- **Microbial contamination detection** using DecontaMiner
- **Comprehensive QC suite**: FastQC, FastQ Screen, Qualimap, DeepTools, Picard, BEDTools
- **Flexible subsampling** for performance optimization
- **Interactive visualizations** and MultiQC reporting
- **Sample sheet support** for batch processing with per-sample strain assignment
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
- [To do](#todo)
- [Contact](#contact)

## Requirements

### Software Dependencies

- **Nextflow** ≥ 23.04.0
- **Java** ≥ 11
- **Python** ≥ 3.7
- **R** ≥ 4.0

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
git clone https://gitlab.epfl.ch/abadreddine/conterminator.git
cd conterminator
```

### Configure Tool Paths

Edit `nextflow.config` to specify paths to installed tools and reference data:

```groovy
params {
    // Tool paths
    star_bin = "/path/to/STAR"
    samtools_bin = "/path/to/samtools"
    fastqc_bin = "/path/to/fastqc"
    // ... (see nextflow.config for all options)
    
    // Reference directories
    strains_base_dir = "/path/to/pseudogenomes"
    star_index_dir = "/path/to/star/indices"
    contamination_blast_dbs = "/path/to/blast/databases"
}
```

### Test Installation

```bash
nextflow run main_full.nf --help
```

## Quick Start

### Minimal Example

```bash
nextflow run main_full.nf \
    --sample_sheet samples.tsv \
    --outdir results
```

### Complete Example with Options

```bash
nextflow run main_full.nf \
    --sample_sheet samples.tsv \
    --outdir results_full \
    --subset_for_fastq_qc true \
    --subset_fastq_qc_reads 100000 \
    --subset_bam_for_qc true \
    --bam_qc_subset_mapped 200000 \
    --max_parallel_samples 4
```

## Input Format

### Sample Sheet (Required)

Create a tab-separated file (`samples.tsv`) with the following format:

```tsv
sample	strain	path
sample1	C57BL_6J	/data/fastq/sample1/
sample2	DBA_2J	/data/fastq/sample2/
sample3	C57BL_6J	/data/fastq/sample3/
```

**Requirements:**
- Header row is mandatory
- Three columns: `sample`, `strain`, `path`
- `sample`: Unique sample identifier
- `strain`: Reference strain name (must exist in `strains_base_dir`)
- `path`: Directory containing paired-end FASTQ files

**Supported FASTQ naming patterns:**
- `*_{1,2}.fq.gz` / `*_{1,2}.fastq.gz`
- `*_{R1,R2}.fq.gz` / `*_{R1,R2}.fastq.gz`
- Pattern auto-detection or specify with `--fastq_pattern`

## Usage

### Basic Commands

```bash
# Show help
nextflow run main_full.nf --help

# Run with defaults
nextflow run main_full.nf --sample_sheet samples.tsv

# Resume failed run
nextflow run main_full.nf --sample_sheet samples.tsv -resume

# Run with specific profile
nextflow run main_full.nf --sample_sheet samples.tsv -profile cluster
```

### Key Parameters

#### Input/Output
```bash
--sample_sheet <file>      # Sample sheet with strain assignments (required)
--outdir <path>            # Output directory (default: results)
--fastq_pattern <pattern>  # FASTQ file pattern (auto-detected)
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

### Resume Failed Runs

```bash
# Nextflow automatically caches completed tasks
nextflow run main_full.nf --sample_sheet samples.tsv -resume
```

### Check Pipeline Version

```bash
nextflow run main_full.nf --version
```

## Execution Profiles

The pipeline supports multiple execution profiles (configure in `nextflow.config`):

```bash
# Local execution
nextflow run main_full.nf -profile local

# With Singularity (not ready yet)
nextflow run main_full.nf -profile singularity
```

## To do

- Create the singularity image.

## Contact

**Author:** Alaa Badreddine

**Project Repository:** https://gitlab.epfl.ch/abadreddine/conterminator

For bug reports and feature requests, please open an issue on GitLab.

---

## Version History

- **v1.0** (2025) - Initial release with full QC and contamination detection suite
