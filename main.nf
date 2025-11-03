#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def content = new File('name.txt').text

println content


def helpMessage() {
    log.info"""
    ======================
      Help Documentation
    ======================
    A pipeline for identifying contamination in RNA-seq data
    
    USAGE:
        nextflow run main.nf --sample_sheet <file.tsv> [options]
        nextflow run main.nf --help
    
    REQUIRED PARAMETERS:
      --sample_sheet <file>       Tab-separated file with columns: sample, strain, path
                                  Example: samples.tsv
    
    INPUT/OUTPUT:
      --outdir <path>             Output directory (default: results)
      --fastq_pattern <pattern>   FASTQ file pattern (auto-detected if not specified)
                                  Examples: "*_{1,2}.fq.gz", "*_{R1,R2}.fastq.gz"
    
    REFERENCE FILES:
      --strains_base_dir <path>   Directory with strain-specific references
                                  (default: /home/abadreddine/.../HDP_pseudogenomes)
      --star_index_dir <path>     Directory for STAR indices (default: /mnt/sas/Data/References/Mus)
    
    STAR ALIGNMENT:
      --run_star_alignment        Enable STAR alignment (default: true)
      --strain <strain(s)>        Comma-separated strain list (alternative to sample sheet)
                                  Example: "C57BL_6J,DBA_2J"
      --star_threads <int>        STAR alignment threads (default: 20)
      --star_index_threads <int>  STAR index building threads (default: 20)
      --star_index_mem <string>   STAR index memory (default: "32 GB")
    
    SUBSAMPLING OPTIONS (Performance Tuning):
      --subset_for_fastq_qc       Subsample FASTQs for QC (default: true)
      --subset_fastq_qc_reads     Number of reads for FASTQ QC (default: 100000)
      
      --subset_for_star           Subsample FASTQs for STAR alignment (default: false)
      --subset_star_reads         Number of reads for STAR (default: 200000)
      
      --subset_bam_for_qc         Subsample BAMs for QC tools (default: true)
      --bam_qc_subset_mapped      Number of mapped reads for BAM QC (default: 2000000)
      
      --subset_mapped_for_blast   Subsample mapped reads for BLAST (default: true)
      --mapped_subset_reads       Number of mapped reads for BLAST (default: 100000)
      
      --subset_unmapped_for_blast Subsample unmapped for BLAST (default: true)
      --unmapped_subset_reads     Number of unmapped reads for BLAST (default: 100000)
      
      --subset_unmapped_for_decontaminer  Subsample unmapped for DecontaMiner (default: true)
      --subset_seed <int>         Random seed for reproducibility (default: 100)
    
    QC TOOLS:
      --run_fastqc                Run FastQC on input FASTQs (default: true)
      --run_fastq_screen          Run FastQ Screen (default: true)
      --fastq_screen_conf         FastQ Screen config file (default: provided)
      --fastq_screen_threads      FastQ Screen threads (default: 8)
      
      --run_deeptools             Run DeepTools GC bias (default: true)
      --run_picard_gc             Run Picard GC bias (default: true)
      --run_bedtools_gc           Run BEDTools GC coverage (default: true)
      --windowsize <int>          BEDTools window size (default: 500)
      --window_step <int>         BEDTools window step (default: 250)
      
      --run_mapinsights           Run MapInsights (default: true)
      --mapinsights_opts          MapInsights options (default: "")
      
      --run_qualimap              Run Qualimap (default: true)
      --qualimap_mode             Qualimap mode: "bamqc", "rnaseq", or "both" (default: "both")
      --qualimap_genome           Reference genome name (default: "mm10")
      --qualimap_java_mem         Java memory allocation (default: "10G")
      --qualimap_threads          Qualimap threads (default: 8)
      --qualimap_protocol         Protocol type (default: "strand-specific-reverse")
      --qualimap_skip_dup         Skip duplicates (default: false)
    
    CONTAMINATION ANALYSIS:
      --run_contamination_check   Enable contamination checking (default: true)
      --contamination_blast_dbs   Path to BLAST database(s) (default: /mnt/sas/Tools/blast_databases/)
      --contamination_evalue      BLAST e-value threshold (default: "1e-10")
      --contamination_min_reads   Minimum reads to report (default: 100)
      --contamination_gc_min      GC content minimum (default: 0.60)
      --contamination_gc_max      GC content maximum (default: 1.0)
      --contamination_mapq        MAPQ threshold (default: 10)
      --blast_max_parallel        Max parallel BLAST jobs (default: 10)
      
      --run_decontaminer          Run DecontaMiner (default: true)
      --decontaminer_dir          DecontaMiner installation directory
      --decontaminer_config       DecontaMiner configuration file (default: provided)
      --decontaminer_pairing      Pairing mode: "P" or "S" (default: "P")
      --decontaminer_organisms    Organism codes (default: "bfv")
      --decontaminer_format       Input format (default: "bam")
      --decontaminer_quality_filter  Quality filter (default: "yes")
      --decontaminer_ribo_filter  Ribosomal filter (default: "yes")
    
    REPORTING:
      --run_multiqc               Run MultiQC (default: true)
      --run_blast_plots           Generate BLAST plots (default: true)
      --multiqc_config            MultiQC config file (default: null)
    
    RESOURCES:
      --max_cpus <int>            Maximum CPUs (default: 16)
      --max_mem <string>          Maximum memory (default: "32 GB")
      --max_time <string>         Maximum time (default: "24h")
      --max_parallel_samples      Maximum parallel samples (default: 4)
    
    DEFAULT RUN (Recommended Settings):
       
       nextflow run main_full.nf \\
       --sample_sheet input_test/sample_sheet.tsv \\
       --outdir results_test
       
       This runs with defaults:
       - FASTQ QC: 100k reads
       - STAR: full FASTQs
       - BAM QC: 200k reads
       - BLAST: 100k reads
       - All QC tools enabled

       Note: sample sheet must include strain column
       References are automatically prepared per-strain from the strains_base_dir
    
    ========================================
    SAMPLE SHEET FORMAT
    ========================================
    
    Tab-separated file (TSV) with header:
    
    sample	strain	path
    sample1	C57BL_6J	/data/fastq/sample1/
    sample2	DBA_2J	/data/fastq/sample2/
    sample3	C57BL_6J	/data/fastq/sample3/
    
    Required columns:
    - sample: Unique sample identifier
    - strain: Reference strain name (must exist in strains_base_dir)
    - path: Directory containing paired-end FASTQ files
    
    Notes:
    - Header row is mandatory
    - Each sample is aligned to its specified strain
    - FASTQ pattern will be auto-detected (or use --fastq_pattern)
    - Supported patterns: *_{1,2}.fq.gz, *_{R1,R2}.fastq.gz, etc.
    
    ========================================
    PERFORMANCE TIPS
    ========================================
    
    1. For testing/development:
       - Use heavy subsampling (10K-50K reads)
       - Disable slow tools if not needed
    
    2. For production with many samples:
       - Use default subsampling (balances speed/accuracy)
       - Adjust --max_parallel_samples based on resources
       - Consider --blast_max_parallel for BLAST jobs
    
    3. For comprehensive analysis:
       - Disable subsampling for critical samples
       - Increase resource allocations (--max_cpus, --max_mem)
       - Enable all QC tools
    
    4. For contamination-focused runs:
       - Keep BLAST/DecontaMiner enabled
       - Disable unnecessary BAM QC tools
       - Use moderate subsampling on unmapped reads
    
    ========================================
    TOOL PATHS
    ========================================
    
    All tool paths are pre-configured in nextflow.config.
    Override if needed using command line parameters:
    
      --samtools_bin <path>
      --star_bin <path>
      --blastn_bin <path>
      --fastqc_bin <path>
      (and more... see nextflow.config)
    
    ========================================
    EXECUTION PROFILES
    ========================================
    
    Run with different execution modes:
    
      nextflow run main.nf -profile local      # Local execution
      nextflow run main.nf -profile cluster    # SLURM cluster
      nextflow run main.nf -profile singularity # With containers
    
    ========================================
    OUTPUT STRUCTURE
    ========================================
    
    results/
    ├── Input
    │   ├── reference
    │   ├── strain_references
    │   ├── subsampled_bams
    │   ├── subsampled_fastq
    │   ├── subsampled_unmapped_blast
    │   └── subsampled_unmapped_decontaminer
    ├── Output
    │   ├── bedtools_gc
    │   ├── contamination_check
    │   ├── decontaminer
    │   ├── deeptools_gc_bias
    │   ├── fastqc
    │   ├── fastq_screen
    │   ├── mapinsights
    │   ├── multiqc
    │   ├── picard_gc_bias
    │   └── qualimap
    └── Temporary
        ├── decontaminer
        └── star_alignment
    
    ========================================
    MORE INFORMATION
    ========================================
    
    Project: Conterminator v1.0
    Author:  Alaa Badreddine
    Home:    https://gitlab.epfl.ch/abadreddine/conterminator
    
    For detailed documentation and issue reporting,
    visit the project repository on GitLab.
    
    ========================================
    """.stripIndent()
}

// Check for help parameter
if (params.help) {
    helpMessage()
    exit 0
}


// ====================
// Strain Management
// ====================

def detectFastqPattern(dir) {
    if (params.fastq_pattern) {
        return params.fastq_pattern
    }
    
    def patterns = [
        "*_{1,2}.fq.gz",
        "*_{1,2}.fq",
        "*_{1,2}.fastq.gz",
        "*_{1,2}.fastq",
        "*_{R1,R2}.fq.gz",
        "*_{R1,R2}.fq",
        "*_{R1,R2}.fastq.gz",
        "*_{R1,R2}.fastq",
        "*.{1,2}.fq.gz",
        "*.{1,2}.fq",
        "*.{1,2}.fastq.gz",
        "*.{1,2}.fastq"
    ]
    
    for (pattern in patterns) {
        def files = file("${dir}/${pattern}")
        if (files && files.size() >= 2) {
            println "Auto-detected FASTQ pattern: ${pattern}"
            return pattern
        }
    }
    
    return null
}

def getAvailableStrains() {
    def strainsDir = new File(params.strains_base_dir)
    if (!strainsDir.exists()) {
        return []
    }
    return strainsDir.listFiles()
        .findAll { it.isDirectory() }
        .collect { it.name }
        .sort()
}

def parseSampleSheet() {
    if (!params.sample_sheet) {
        return [:]
    }
    
    def sampleSheet = [:]
    def sheetFile = file(params.sample_sheet)
    
    if (!sheetFile.exists()) {
        println "ERROR: Sample sheet file not found: ${params.sample_sheet}"
        System.exit(1)
    }
    
    def lines = sheetFile.readLines()
    def header = lines[0].split('\t')
    
    // Check for 3 columns
    if (header.size() < 3 || header[0].toLowerCase() != 'sample' || 
        header[1].toLowerCase() != 'strain' || header[2].toLowerCase() != 'path') {
        println """
        ========================================
        ERROR: Invalid sample sheet format
        ========================================
        Expected header: sample\tstrain\tpath
        Found: ${header.join('\t')}
        
        Example format:
        sample\tstrain\tpath
        43\tC57BL_6J\t/path/to/fastq/
        ========================================
        """
        System.exit(1)
    }
    
    // Store both strain and path
    lines[1..-1].each { line ->
        if (line.trim() && !line.startsWith('#')) {
            def fields = line.split('\t')
            if (fields.size() >= 3) {
                def sample = fields[0].trim()
                def strain = fields[1].trim()
                def path = fields[2].trim()
                sampleSheet[sample] = [strain: strain, path: path]
            }
        }
    }
    
    return sampleSheet
}

def getStrainsFromSampleSheet(sampleSheet) {
    return sampleSheet.values().collect { it.strain }.unique().sort()
}

def parseStrains() {
    if (!params.strain) {
        return []
    }
    return params.strain.split(',').collect { it.trim() }
}

// Validate strains at startup
def availableStrains = getAvailableStrains()
def requestedStrains = []
def sampleToStrain = [:]

if (params.sample_sheet && params.run_star_alignment) {
    if (params.sample_sheet && params.strain) {
        println """
        ========================================
        ERROR: Conflicting parameters
        ========================================
        Please provide EITHER --sample_sheet OR --strain, not both.
        ========================================
        """
        System.exit(1)
    }
    
    if (params.sample_sheet) {
        sampleToStrain = parseSampleSheet()
        requestedStrains = getStrainsFromSampleSheet(sampleToStrain)
        
        println """
        ========================================
        Sample Sheet Configuration
        ========================================
        Sample sheet: ${params.sample_sheet}
        Sample-strain mappings:
        ${sampleToStrain.collect { k, v -> "  ${k} → ${v}" }.join('\n')}
        
        Unique strains: ${requestedStrains.join(', ')}
        ========================================
        """
    } else if (params.strain) {
        requestedStrains = parseStrains()
        
        println """
        ========================================
        Strain Configuration
        ========================================
        Requested strains: ${requestedStrains.join(', ')}
        Mode: All samples will be aligned to ALL specified strains
        ========================================
        """
    } else {
        println """
        ========================================
        ERROR: No strain specified
        ========================================
        Please specify strains using EITHER:
        
        Option 1: Sample sheet (recommended)
          --sample_sheet samples.tsv
        
        Option 2: Strain list
          --strain C57BL_6J
          --strain "C57BL_6J,DBA_2J"
        
        Available strains (${availableStrains.size()}):
        ========================================
        ${availableStrains.join('\n        ')}
        ========================================
        """
        System.exit(1)
    }
    
    def invalidStrains = requestedStrains.findAll { !availableStrains.contains(it) }
    if (invalidStrains) {
        println """
        ========================================
        ERROR: Invalid strain(s) specified
        ========================================
        The following strain(s) are not available:
        ${invalidStrains.join(', ')}
        
        Available strains (${availableStrains.size()}):
        ========================================
        ${availableStrains.join('\n        ')}
        ========================================
        """
        System.exit(1)
    }
}

def discoverBlastDatabases(blast_dir) {
    def databases = []
    def dir = new File(blast_dir)
    
    if (!dir.exists() || !dir.isDirectory()) {
        return databases
    }
    
    def nhrFiles = dir.listFiles().findAll { it.name.endsWith('.nhr') }
    
    nhrFiles.each { file ->
        def baseName = file.name.replaceAll(/\.\d+\.nhr$/, '').replaceAll(/\.nhr$/, '')
        if (!databases.contains(baseName)) {
            databases.add(baseName)
        }
    }
    
    return databases.unique().sort()
}

// ====================
// Input Validation
// ====================

if (!params.sample_sheet) {
    exit 1, "ERROR: Please provide input via --sample_sheet"
}

// QC tools now use per-strain references from sample sheet
def need_ref = params.run_deeptools || params.run_picard_gc || params.run_bedtools_gc || params.run_contamination_check || params.run_mapinsights

if (params.run_decontaminer) {
    if (!params.decontaminer_config) {
        exit 1, "ERROR: DecontaMiner requires --decontaminer_config"
    }
    if (!file(params.decontaminer_config).exists()) {
        exit 1, "ERROR: DecontaMiner config file not found: ${params.decontaminer_config}"
    }
}

// Print configuration
println """
========================================
FLEXIBLE SUBSAMPLING CONFIGURATION
========================================
FASTQ QC Subsampling:    ${params.subset_for_fastq_qc} (${params.subset_fastq_qc_reads} reads)
STAR Subsampling:        ${params.subset_for_star} (${params.subset_star_reads} reads)
BAM QC Subsampling:      ${params.subset_bam_for_qc} (${params.bam_qc_subset_mapped} reads)
Mapped BLAST Subset:     ${params.subset_mapped_for_blast} (${params.mapped_subset_reads} reads)
Unmapped BLAST Subset:   ${params.subset_unmapped_for_blast} (${params.unmapped_subset_reads} reads)
Unmapped Decon Subset:   ${params.subset_unmapped_for_decontaminer} (${params.unmapped_subset_reads} reads)
========================================
"""

// ====================
// Processes
// ====================

include { WRITE_PIPELINE_INFO } from './modules/01.WRITE_PIPELINE_INFO.nf'
include { COPY_SAMPLE_SHEET } from './modules/02.COPY_SAMPLE_SHEET.nf'
include { SUBSET_FASTQ_FOR_QC } from './modules/03.SUBSET_FASTQ_FOR_QC.nf'
include { SUBSET_FASTQ_FOR_STAR } from './modules/04.SUBSET_FASTQ_FOR_STAR.nf'
include { PREPARE_STRAIN_REFERENCE } from './modules/05.PREPARE_STRAIN_REFERENCE.nf'
include { BUILD_STAR_INDEX } from './modules/06.BUILD_STAR_INDEX.nf'
include { STAR_ALIGN } from './modules/07.STAR_ALIGN.nf'
include { PREP_STRAIN_REFERENCE_FOR_QC } from './modules/08.PREP_STRAIN_REFERENCE_FOR_QC.nf'
include { SUBSET_BAM_FOR_QC } from './modules/09.SUBSET_BAM_FOR_QC.nf'
include { SUBSET_UNMAPPED_FOR_BLAST } from './modules/10.SUBSET_UNMAPPED_FOR_BLAST.nf'
include { SUBSET_UNMAPPED_FOR_DECONTAMINER } from './modules/11.SUBSET_UNMAPPED_FOR_DECONTAMINER.nf'
include { DEEPTOOLS_GC } from './modules/12.DEEPTOOLS_GC.nf'
include { PICARD_GC_BIAS } from './modules/13.PICARD_GC_BIAS.nf'
include { BEDTOOLS_GC_COVERAGE } from './modules/14.BEDTOOLS_GC_COVERAGE.nf'
include { QUALIMAP_BAMQC } from './modules/15.QUALIMAP_BAMQC.nf'
include { QUALIMAP_RNASEQ } from './modules/16.QUALIMAP_RNASEQ.nf'
include { MAPINSIGHTS } from './modules/17.MAPINSIGHTS.nf'
include { FASTQ_SCREEN } from './modules/18.FASTQ_SCREEN.nf'
include { FASTQC_INPUT_FASTQ } from './modules/19.FASTQC_INPUT_FASTQ.nf'
include { DECONTAMINER_STEP1_STAR_MAPPED } from './modules/20.DECONTAMINER_STEP1_STAR_MAPPED.nf'
include { DECONTAMINER_STEP1_STAR_UNMAPPED } from './modules/21.DECONTAMINER_STEP1_STAR_UNMAPPED.nf'
include { DECONTAMINER_STEP2 } from './modules/22.DECONTAMINER_STEP2.nf'
include { DECONTAMINER_STEP3 } from './modules/23.DECONTAMINER_STEP3.nf'
include { FASTQC_MAPPED_BAM } from './modules/24.FASTQC_MAPPED_BAM.nf'
include { FASTQC_UNMAPPED_FASTQ } from './modules/25.FASTQC_UNMAPPED_FASTQ.nf'
include { BLAST_MAPPED_READS_MULTI } from './modules/26.BLAST_MAPPED_READS_MULTI.nf'
include { BLAST_UNMAPPED_READS_MULTI } from './modules/27.BLAST_UNMAPPED_READS_MULTI.nf'
include { BLAST_PLOT_CHARTS } from './modules/28.BLAST_PLOT_CHARTS.nf'
include { CONTAMINATION_FINAL_REPORT_MULTI } from './modules/29.CONTAMINATION_FINAL_REPORT_MULTI.nf'
include { MULTIQC } from './modules/30.MULTIQC.nf'


// ====================
// Main Workflow
// ====================
workflow {
    // ===================================
    // PIPELINE METADATA
    // ===================================

    // Write pipeline information to output directory
    WRITE_PIPELINE_INFO()

    // ===================================
    // PROCESS EXECUTION TRACKING FLAGS
    // ===================================
    def ran_fastqc_input = false
    def ran_fastq_screen = false
    def ran_star = false
    def ran_picard = false
    def ran_qualimap_bamqc = false
    def ran_qualimap_rnaseq = false
    def ran_fastqc_mapped = false
    def ran_fastqc_unmapped = false
    def ran_deeptools = false
    def ran_bedtools = false

    // ===================================
    // SECTION 1: INPUT HANDLING
    // ===================================
    
    if (params.sample_sheet) {
        // MODE 1: Sample sheet with individual paths per sample
        println "Using sample sheet with per-sample paths"
        
        def fastq_pattern = params.fastq_pattern
        if (!fastq_pattern) {
            def firstPath = sampleToStrain.values()[0].path
            fastq_pattern = detectFastqPattern(firstPath)
            if (!fastq_pattern) {
                println "ERROR: Could not detect FASTQ pattern. Please specify --fastq_pattern"
                System.exit(1)
            }
        }
        
        // Create channel from sample sheet paths
        def fastqList = []
        sampleToStrain.each { sample, info ->
            def samplePath = info.path
            def pattern = "${samplePath}/${fastq_pattern}"
            
            def files = file(pattern)
            if (files instanceof List && files.size() >= 2) {
                // Sort to ensure consistent R1, R2 pairing
                files = files.sort()
                fastqList.add([sample, files[0], files[1]])
            } else {
                println "WARNING: Could not find FASTQ pair for sample ${sample} at ${samplePath}"
                println "  Pattern: ${pattern}"
            }
        }
        
        if (fastqList.isEmpty()) {
            println "ERROR: No FASTQ files found for any samples in the sample sheet"
            System.exit(1)
        }
        
        Channel
            .fromList(fastqList)
            .set { raw_fastq_ch }
        
    } else {
        println """
        ========================================
        ERROR: No input provided
        ========================================
        Please provide input via:
          --sample_sheet <file.tsv> (with path column)
        ========================================
        """
        System.exit(1)
    }
    
    // ===================================
    // SECTION 2: FASTQ PROCESSING (if FASTQ input)
    // ===================================
    
    if (params.sample_sheet) {
        // Copy sample sheet to output directory
        COPY_SAMPLE_SHEET(file(params.sample_sheet))
        
        // SSubset for FastQ QC
        if (params.subset_for_fastq_qc) {
            SUBSET_FASTQ_FOR_QC(raw_fastq_ch)
            fastq_for_screen = SUBSET_FASTQ_FOR_QC.out.subsampled_fastq
        } else {
            fastq_for_screen = raw_fastq_ch
        }
        
        // Subset (or not) for STAR
        if (params.subset_for_star) {
            SUBSET_FASTQ_FOR_STAR(raw_fastq_ch)
            fastq_for_alignment = SUBSET_FASTQ_FOR_STAR.out.subsampled_fastq
        } else {
            fastq_for_alignment = raw_fastq_ch
        }
        
        // FastQC on input FASTQs
        if (params.run_fastqc) {
            FASTQC_INPUT_FASTQ(fastq_for_screen)
            ran_fastqc_input = true
        }
        
        // FastQ Screen
        if (params.run_fastq_screen) {
            FASTQ_SCREEN(fastq_for_screen)
            ran_fastq_screen = true
        }
        
        // ===================================
        // SECTION 3: STAR ALIGNMENT
        // ===================================
        
        if (params.run_star_alignment) {
            // Sample sheet mode - each sample to its specific strain
            strain_ch = Channel.from(requestedStrains)
            PREPARE_STRAIN_REFERENCE(strain_ch)

            // Prepare per-strain references for QC tools
            if (need_ref) {
                PREP_STRAIN_REFERENCE_FOR_QC(PREPARE_STRAIN_REFERENCE.out.references)
                strain_qc_refs = PREP_STRAIN_REFERENCE_FOR_QC.out.qc_references
            }

            strain_index_status = PREPARE_STRAIN_REFERENCE.out.references
                .map { strain, fasta, gtf ->
                    def indexDir = file("${params.star_index_dir}/${strain}")
                    def saFile = file("${params.star_index_dir}/${strain}/SA")
                    def genomeFile = file("${params.star_index_dir}/${strain}/Genome")
                    def indexExists = indexDir.exists() && saFile.exists() && genomeFile.exists()
                    
                    if (indexExists) {
                        println "OK: Found existing STAR index for ${strain}"
                    } else {
                        println "Warning: STAR index not found for ${strain}, building..."
                    }
                    
                    return tuple(strain, fasta, gtf, indexDir, indexExists)
                }
            
            strains_needing_index = strain_index_status
                .filter { strain, fasta, gtf, indexDir, exists -> !exists }
                .map { strain, fasta, gtf, indexDir, exists -> tuple(strain, fasta, gtf) }
            
            BUILD_STAR_INDEX(strains_needing_index)
            
            existing_indexes = strain_index_status
                .filter { strain, fasta, gtf, indexDir, exists -> exists }
                .map { strain, fasta, gtf, indexDir, exists -> tuple(strain, indexDir) }
            
            all_indexes = existing_indexes.mix(BUILD_STAR_INDEX.out.star_index)

            fastq_with_strain = fastq_for_alignment
                .map { sample, fq1, fq2 ->
                    def strain = sampleToStrain[sample]?.strain
                    if (!strain) {
                        println "WARNING: Sample ${sample} not found in sample sheet, skipping..."
                        return null
                    }
                    return tuple(strain, sample, fq1, fq2)
                }
                .filter { it != null }

            fastq_strain_combinations = fastq_with_strain
                .combine(all_indexes, by: 0)
                .map { strain, sample, fq1, fq2, index_dir ->
                    tuple(sample, fq1, fq2, strain, index_dir)
                }
            
            // Run STAR alignment
            STAR_ALIGN(fastq_strain_combinations)
            ran_star = true
            
            // Save unmapped FASTQ for later use
            star_unmapped_fastq = STAR_ALIGN.out.unmapped_fastq
        }
    }

    // ===================================
    // SECTION 4: BAM QC TOOLS
    // ===================================
    
    if (params.sample_sheet && params.run_star_alignment) {
        
        if (params.subset_bam_for_qc) {
            SUBSET_BAM_FOR_QC(STAR_ALIGN.out.aligned_bam)
            bams_for_qc = SUBSET_BAM_FOR_QC.out.subsampled_bam
            qc_bams = bams_for_qc.map { sample, strain, bam, bai -> 
                tuple(sample, strain, bam, bai, "subset") 
            }
        } else {
            bams_for_qc = STAR_ALIGN.out.aligned_bam
            qc_bams = bams_for_qc.map { sample, strain, bam, bai -> 
                tuple(sample, strain, bam, bai, "full") 
            }
        }

        // Join BAMs with their strain-specific references for QC tools
        if (need_ref) {
            qc_bams_with_refs = qc_bams
                .map { sample, strain, bam, bai, type -> tuple(strain, sample, bam, bai, type) }
                .combine(strain_qc_refs, by: 0)
                .map { strain, sample, bam, bai, type, fa, fai, dict, bed, eff_size, twobit ->
                    tuple(sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit)
                }

            // Run QC tools
            if (params.run_deeptools) {
                qc_bams_for_deeptools = qc_bams_with_refs.map { sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit ->
                    tuple("${sample}_${strain}", bam, bai, type, twobit, eff_size)
                }
                // Group by strain and run once per strain
                qc_bams_for_deeptools
                    .map { id, bam, bai, type, twobit, eff_size -> tuple(eff_size.name, id, bam, bai, type, twobit, eff_size) }
                    .groupTuple(by: 0)
                    .flatMap { key, ids, bams, bais, types, twobits, eff_sizes ->
                        // Zip all inputs together
                        [ids, bams, bais, types].transpose().collect { id, bam, bai, type ->
                            tuple(id, bam, bai, type, twobits[0], eff_sizes[0])
                        }
                    }
                    .set { deeptools_grouped }
                DEEPTOOLS_GC(deeptools_grouped.map { id, bam, bai, type, twobit, eff_size -> tuple(id, bam, bai, type) },
                            deeptools_grouped.first().map { id, bam, bai, type, twobit, eff_size -> twobit },
                            deeptools_grouped.first().map { id, bam, bai, type, twobit, eff_size -> eff_size })
                ran_deeptools = true
            }

            if (params.run_picard_gc) {
                qc_bams_for_picard = qc_bams_with_refs.map { sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit ->
                    tuple("${sample}_${strain}", bam, bai, type, fa)
                }
                // Group by reference (strain) and run once per strain
                qc_bams_for_picard
                    .map { id, bam, bai, type, fa -> tuple(fa.name, id, bam, bai, type, fa) }
                    .groupTuple(by: 0)
                    .flatMap { key, ids, bams, bais, types, fas ->
                        [ids, bams, bais, types].transpose().collect { id, bam, bai, type ->
                            tuple(id, bam, bai, type, fas[0])
                        }
                    }
                    .set { picard_grouped }
                PICARD_GC_BIAS(picard_grouped.map { id, bam, bai, type, fa -> tuple(id, bam, bai, type) },
                              picard_grouped.first().map { id, bam, bai, type, fa -> fa })
                ran_picard = true
            }

            if (params.run_bedtools_gc) {
                qc_bams_for_bedtools = qc_bams_with_refs.map { sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit ->
                    tuple("${sample}_${strain}", bam, bai, type, fa, bed)
                }
                // Group by reference (strain) and run once per strain
                qc_bams_for_bedtools
                    .map { id, bam, bai, type, fa, bed -> tuple(fa.name, id, bam, bai, type, fa, bed) }
                    .groupTuple(by: 0)
                    .flatMap { key, ids, bams, bais, types, fas, beds ->
                        [ids, bams, bais, types].transpose().collect { id, bam, bai, type ->
                            tuple(id, bam, bai, type, fas[0], beds[0])
                        }
                    }
                    .set { bedtools_grouped }
                BEDTOOLS_GC_COVERAGE(bedtools_grouped.map { id, bam, bai, type, fa, bed -> tuple(id, bam, bai, type) },
                                    bedtools_grouped.first().map { id, bam, bai, type, fa, bed -> fa },
                                    bedtools_grouped.first().map { id, bam, bai, type, fa, bed -> bed })
                ran_bedtools = true
            }

            if (params.run_mapinsights) {
                qc_bams_for_mapinsights = qc_bams_with_refs.map { sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit ->
                    tuple("${sample}_${strain}", bam, bai, type, fa)
                }
                // Group by reference (strain) and run once per strain
                qc_bams_for_mapinsights
                    .map { id, bam, bai, type, fa -> tuple(fa.name, id, bam, bai, type, fa) }
                    .groupTuple(by: 0)
                    .flatMap { key, ids, bams, bais, types, fas ->
                        [ids, bams, bais, types].transpose().collect { id, bam, bai, type ->
                            tuple(id, bam, bai, type, fas[0])
                        }
                    }
                    .set { mapinsights_grouped }
                MAPINSIGHTS(mapinsights_grouped.map { id, bam, bai, type, fa -> tuple(id, bam, bai, type) },
                           mapinsights_grouped.first().map { id, bam, bai, type, fa -> fa })
            }

            if (params.run_qualimap) {
                qc_bams_for_qualimap = qc_bams_with_refs.map { sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit ->
                    tuple("${sample}_${strain}", strain, bam, bai, type, fa)
                }

                if (params.qualimap_mode == "bamqc" || params.qualimap_mode == "both") {
                    // Group by reference (strain) and run once per strain
                    qc_bams_for_qualimap
                        .map { id, strain, bam, bai, type, fa -> tuple(fa.name, id, strain, bam, bai, type, fa) }
                        .groupTuple(by: 0)
                        .flatMap { key, ids, strains, bams, bais, types, fas ->
                            [ids, strains, bams, bais, types].transpose().collect { id, strain, bam, bai, type ->
                                tuple(id, strain, bam, bai, type, fas[0])
                            }
                        }
                        .set { qualimap_grouped }
                    QUALIMAP_BAMQC(qualimap_grouped.map { id, strain, bam, bai, type, fa -> tuple(id, strain, bam, bai, type) },
                                  qualimap_grouped.first().map { id, strain, bam, bai, type, fa -> fa })
                    ran_qualimap_bamqc = true
                }

                if (params.qualimap_mode == "rnaseq" || params.qualimap_mode == "both") {
                    QUALIMAP_RNASEQ(qc_bams_for_qualimap.map { id, strain, bam, bai, type, fa -> tuple(id, strain, bam, bai, type) })
                    ran_qualimap_rnaseq = true
                }
            }
        }
    }

    // ===================================
    // SECTION 5: DECONTAMINER
    // ===================================
    
    if (params.run_decontaminer) {
        // Initialize empty channel for decontaminer outputs
        decontaminer_all_outputs = Channel.empty()
        
        if (params.sample_sheet && params.run_star_alignment) {
            // Use STAR outputs
            
            // Run on mapped reads
            decontaminer_mapped_input = bams_for_qc
                .map { sample, strain, bam, bai -> 
                    tuple("${sample}_${strain}", bam, bai, "mapped") 
                }
            
            // Subset unmapped for DecontaMiner
            if (params.subset_unmapped_for_decontaminer) {
                SUBSET_UNMAPPED_FOR_DECONTAMINER(star_unmapped_fastq)
                unmapped_for_decontaminer = SUBSET_UNMAPPED_FOR_DECONTAMINER.out.subsampled_fastq
            } else {
                unmapped_for_decontaminer = star_unmapped_fastq
            }
            
            decontaminer_unmapped_input = unmapped_for_decontaminer
                .map { sample, strain, r1, r2 -> 
                    tuple("${sample}_${strain}", r1, r2, "unmapped") 
                }
            
            
            DECONTAMINER_STEP1_STAR_MAPPED(decontaminer_mapped_input)
            DECONTAMINER_STEP1_STAR_UNMAPPED(decontaminer_unmapped_input)
            
            // Combine outputs for steps 2 and 3
            decontaminer_all_outputs = DECONTAMINER_STEP1_STAR_MAPPED.out.output_dir
                .mix(DECONTAMINER_STEP1_STAR_UNMAPPED.out.output_dir)
            
        }
        
        // Run steps 2 and 3 - they will automatically wait for step 1 outputs
        DECONTAMINER_STEP2(decontaminer_all_outputs)
        DECONTAMINER_STEP3(DECONTAMINER_STEP2.out.filtered_dir)
    }

    // ===================================
    // SECTION 6: CONTAMINATION CHECK
    // ===================================
    
    if (params.run_contamination_check && (params.sample_sheet && params.run_star_alignment)) {
    
        // Use full BAMs, extract mapped stats
        extract_mapped_input = qc_bams.map { sample, strain, bam, bai, type ->
            tuple("${sample}", "${strain}", bam, bai)
        }
        
        // Use STAR unmapped directly
        unmapped_for_contam_check = star_unmapped_fastq
            .map { sample, strain, r1, r2 -> 
                tuple("${sample}", "${strain}", r1, r2)
            }

        // FastQC on mapped
        if (params.run_fastqc && params.run_star_alignment) {
            FASTQC_MAPPED_BAM(extract_mapped_input)
            ran_fastqc_mapped = true
        }
        
        // FastQC on unmapped
        if (params.run_fastqc && params.run_star_alignment) {
            FASTQC_UNMAPPED_FASTQ(unmapped_for_contam_check)
            ran_fastqc_unmapped = true
        }

        // ===================================
        // SECTION 7: BLAST ANALYSIS
        // ===================================
        
        if (params.run_star_alignment) {
            
            // Multi-database BLAST support
            def blast_source = params.contamination_blast_dbs

            if (blast_source) {
                def blast_dir = new File(blast_source)
                def blast_dbs = []
                
                if (blast_dir.isDirectory() && blast_dir.exists()) {
                    def dbFiles = blast_dir.listFiles().findAll { 
                        it.name.endsWith('.nsq') || it.name.endsWith('.nhr') || it.name.endsWith('.nin')
                    }
                    
                    // Extract base names, handling both single and split databases
                    blast_dbs = dbFiles.collect { 
                        def name = it.name
                        // For split databases (e.g., db.00.nsq), extract base name (db)
                        if (name =~ /\.\d+\.(nsq|nhr|nin)$/) {
                            name = name.replaceAll(/\.\d+\.(nsq|nhr|nin)$/, '')
                        } else {
                            // For single databases (e.g., db.nsq), just remove extension
                            name = name.replaceAll(/\.(nsq|nhr|nin)$/, '')
                        }
                        return name
                    }.unique()
                    
                    if (!blast_dbs.isEmpty()) {
                        println """
    ========================================
    Available BLAST Databases
    ========================================
    ${blast_dbs.collect { "${it} → ${blast_source}/${it}" }.join('\n    ')}
    ========================================
                        """
                    }
                } else if (blast_dir.exists()) {
                    def dbName = blast_source.split('/').last()
                    blast_dbs = [dbName]
                    println """
    ========================================
    Using Single BLAST Database
    ========================================
    Database: ${dbName}
    Path: ${blast_source}
    ========================================
                    """
                } else {
                    println "WARNING: BLAST database path does not exist: ${blast_source}"
                }
                
                if (!blast_dbs.isEmpty()) {
                    // Create database channel
                    if (blast_dir.isDirectory()) {
                        blast_db_ch = Channel.fromList(
                            blast_dbs.collect { db -> tuple(db, "${blast_source}/${db}") }
                        )
                    } else {
                        blast_db_ch = Channel.of(tuple(blast_dbs[0], blast_source))
                    }
                    
                    // ===================================
                    // BLAST ON MAPPED READS
                    // ===================================
                    if (params.subset_mapped_for_blast) {
                        // Combine BAMs with databases
                        // bams_for_qc structure: tuple(sample, strain, bam, bai, type)
                        // blast_db_ch structure: tuple(db_name, db_path)
                        blast_mapped_input_ch = bams_for_qc
                            .combine(blast_db_ch)
                            .map { combined -> 
                                // After combine, we get all elements in a single list
                                // Order: [sample, strain, bam, bai, type, db_name, db_path]
                                def sample = combined[0]
                                def strain = combined[1]
                                def bam = combined[2]
                                def bai = combined[3]
                                def type = combined[4]
                                def db_path = combined[5]
                                def db_name = db_path.split("/").last()
                                
                                // Return properly ordered tuple for the process
                                tuple("${sample}_${strain}", strain, bam, bai, db_name, db_path)
                            }
                        // Run BLAST on mapped reads
                        println "Launching BLAST jobs on mapped reads for ${blast_dbs.size()} databases..."
                        BLAST_MAPPED_READS_MULTI(blast_mapped_input_ch)
                    }

                    // ===================================
                    // BLAST ON UNMAPPED READS
                    // ===================================
                    if (params.subset_unmapped_for_blast) {
                        // Subset unmapped for BLAST
                        SUBSET_UNMAPPED_FOR_BLAST(star_unmapped_fastq)
                        unmapped_for_blast = SUBSET_UNMAPPED_FOR_BLAST.out.subsampled_fastq
                        
                        // Combine unmapped with databases
                        // unmapped_for_blast structure: tuple(sample, strain, r1, r2)
                        // blast_db_ch structure: tuple(db_name, db_path)
                        blast_unmapped_input_ch = unmapped_for_blast
                            .combine(blast_db_ch)
                            .map { combined -> 
                                // After combine: [sample, strain, r1, r2, db_name, db_path]
                                def sample = combined[0]
                                def strain = combined[1]
                                def r1 = combined[2]
                                def r2 = combined[3]
                                def db_name = combined[4]
                                def db_path = combined[5]
                                
                                // Return properly ordered tuple for the process
                                tuple("${sample}_${strain}", r1, r2, db_name, db_path)
                            }
                        
                        // Run BLAST on unmapped reads
                        println "Launching BLAST jobs on unmapped reads for ${blast_dbs.size()} databases..."
                        BLAST_UNMAPPED_READS_MULTI(blast_unmapped_input_ch)
                    }
                    
                    // ===================================
                    // GENERATE PLOTS FOR ALL BLAST RESULTS
                    // ===================================
                    if (params.run_blast_plots) {
                    all_blast_for_plotting = Channel.empty()
                    
                    if (params.subset_mapped_for_blast) {
                        mapped_labeled = BLAST_MAPPED_READS_MULTI.out.for_plotting
                            .map { it + ['mapped'] }
                        all_blast_for_plotting = all_blast_for_plotting.mix(mapped_labeled)
                    }
                    
                    if (params.subset_unmapped_for_blast) {
                        unmapped_labeled = BLAST_UNMAPPED_READS_MULTI.out.for_plotting
                            .map { it + ['unmapped'] }
                        all_blast_for_plotting = all_blast_for_plotting.mix(unmapped_labeled)
                    }
                    
                    BLAST_PLOT_CHARTS(all_blast_for_plotting)
                }
                }
            }
        }
    }

    // ===================================
    // SECTION 8: MULTIQC
    // ===================================

    if (params.run_multiqc) {
        multiqc_files = Channel.empty()
        
        // 1. FastQC on input FASTQs
        if (ran_fastqc_input) {
            multiqc_files = multiqc_files.mix(
                FASTQC_INPUT_FASTQ.out.results.flatten()
            )
        }
        
        // 2. FastQ Screen
        if (ran_fastq_screen) {
            multiqc_files = multiqc_files.mix(
                FASTQ_SCREEN.out.summary.map { it[1] }
            )
        }
        
        // 3. STAR alignment logs
        if (ran_star) {
            multiqc_files = multiqc_files.mix(
                STAR_ALIGN.out.log
            )
        }
        
        // 4. Picard GC metrics
        if (ran_picard) {
            multiqc_files = multiqc_files.mix(
                PICARD_GC_BIAS.out.metrics,
                PICARD_GC_BIAS.out.summary
            )
        }
        
        // 5. Qualimap reports
        if (ran_qualimap_bamqc) {
            multiqc_files = multiqc_files.mix(
                QUALIMAP_BAMQC.out.results.flatten().filter { it.name == 'genome_results.txt' }
            )
        }
        if (ran_qualimap_rnaseq) {
            multiqc_files = multiqc_files.mix(
                QUALIMAP_RNASEQ.out.results.flatten().filter { it.name == 'rnaseq_results.txt' }
            )
        }
        
        // 6. FastQC on mapped BAM
        if (ran_fastqc_mapped) {
            multiqc_files = multiqc_files.mix(
                FASTQC_MAPPED_BAM.out.results.flatten()
            )
        }
        
        // 7. FastQC on unmapped reads - USE FLAG
        if (ran_fastqc_unmapped) {
            multiqc_files = multiqc_files.mix(
                FASTQC_UNMAPPED_FASTQ.out.results_r1.flatten(),
                FASTQC_UNMAPPED_FASTQ.out.results_r2.flatten()
            )
        }
        
        // 8. DeepTools GC bias
        if (ran_deeptools) {
            multiqc_files = multiqc_files.mix(
                DEEPTOOLS_GC.out.freq
            )
        }
        
        // 9. BEDTools coverage
        if (ran_bedtools) {
            multiqc_files = multiqc_files.mix(
                BEDTOOLS_GC_COVERAGE.out.coverage
            )
        }
        
        // Collect all files and run MultiQC
        multiqc_files
            .collect()
            .map { files -> files.findAll { it != null && it.exists() } }
            .filter { it.size() > 0 }
            .set { multiqc_input }
        
        MULTIQC(multiqc_input)
    }
    
} // End of workflow block

// ===================================
// WORKFLOW COMPLETION HANDLERS
// ===================================

workflow.onComplete {
    println """
    ========================================
    QC Pipeline Completed
    ========================================
    Status:   ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    Output:   ${params.outdir}
    Input:    ${params.sample_sheet}
    ${params.sample_sheet ? "Sample Sheet: ${params.sample_sheet}" : ''}
    ${params.strain ? "Strains:  ${params.strain}" : ''}
    
    SUBSAMPLING SUMMARY:
    ${params.subset_for_fastq_qc ? "  ✓ FastQ QC: ${params.subset_fastq_qc_reads} reads" : "  X FastQ QC: Full FASTQs"}
    ${params.subset_for_star ? "  ✓ STAR: ${params.subset_star_reads} reads" : "  X STAR: Full FASTQs"}
    ${params.subset_bam_for_qc ? "  ✓ BAM QC: ${params.bam_qc_subset_mapped} reads" : "  X BAM QC: Full BAMs"}
    ${params.subset_unmapped_for_blast ? "  ✓ BLAST: ${params.unmapped_subset_reads} reads" : "  X BLAST: Full unmapped"}
    ${params.subset_unmapped_for_decontaminer ? "  ✓ DecontaMiner: ${params.unmapped_subset_reads} reads" : "  X DecontaMiner: Full unmapped"}
    
    ${params.run_multiqc ? "To generate MultiQC report:\n    cd ${params.outdir} && multiqc . --force" : ''}
    ========================================
    """
}

workflow.onError {
    println """
    ========================================
    Pipeline Error
    ========================================
    Error: ${workflow.errorMessage}
    ========================================
    """
}