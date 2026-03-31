#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ====================
// Pipeline Parameters
// ====================

// help
params.help = false

// Inputs
params.sample_sheet             = null
params.fastq_pattern            = null

// FLEXIBLE FastQ SUBSETTING
params.subset_for_fastq_qc      = false      // Subset for FastQ Screen
params.subset_fastq_qc_reads    = 100000    // Size for FastQ QC (100k reads)
params.subset_for_star          = false     // Use full FastQ for STAR alignment
params.subset_star_reads        = 100000    // Size for STAR if subsampling (100k reads)
params.subset_seed              = 100       // Random seed for reproducibility

// STAR alignment
params.run_star_alignment       = true
params.strain                   = null
params.strains_base_dir         = null      // Directory with strain-specific pseudogenomes (required for STAR)
params.standard_references_dir  = null      // Directory with standard reference genomes (e.g., GRCm39, GRCm38)
params.star_index_dir           = null      // Directory for STAR indices

// BAM SUBSETTING FOR QC TOOLS
params.subset_bam_for_qc                 = false      // Subset BAMs for QC tools
params.bam_qc_subset_mapped              = 200000   // Mapped reads for QC (200k reads)

// READS SUBSETTING
params.subset_unmapped_for_blast         = true   // Subset for BLAST
params.subset_mapped_for_blast           = true   // Subset for BLAST
params.unmapped_subset_reads             = 100000 // Unmapped reads to keep (100k reads)
params.mapped_subset_reads               = 100000 // Mapped reads to keep (100k reads)

// Tool toggles
params.run_deeptools           = true
params.run_picard_gc           = true
params.run_fastq_screen        = true
params.run_bedtools_gc         = true
params.run_mapinsights         = true
params.run_contamination_check = true
params.run_fastqc              = true

// Deeptools settings
params.blacklist_bed      = null

// BEDTools windowing
params.windowsize         = 500
params.window_step        = 250

// Contamination check
params.contamination_gc_min    = 0.60
params.contamination_gc_max    = 1.0
params.contamination_mapq      = 10
params.contamination_blast_dbs = null       // Path to BLAST database directory (required for contamination check)
params.contamination_evalue    = "1e-10"

// Mapinsights
params.mapinsights_opts   = ""

// Tool paths — set these when running WITHOUT a container (ignored with Singularity)
params.bedtools_bin       = null
params.mapinsights_bin    = null
params.picard_jar         = null
params.samtools_bin       = null
params.blastn_bin         = null
params.deeptools_gcbias   = null
params.facount_bin        = null
params.fa2bit_bin         = null
params.fastq_screen_bin   = null
params.qualimap_bin       = null
params.seqtk_bin          = null
params.star_bin           = null
params.fastqc_bin         = null
params.reformat_bin       = null

// Qualimap settings
params.run_qualimap       = true
params.qualimap_mode      = "both"
params.qualimap_genome    = "mm10"
params.qualimap_protocol  = "strand-specific-reverse"

// FastQ Screen conf
params.fastq_screen_conf   = null           // FastQ Screen config file (required for fastq_screen)

// Plotting and reporting
params.gc_plot_script         =  "${projectDir}/bin/R/plot_gc_content.R"
params.run_blast_plots        = true
params.blast_plot             = "${projectDir}/bin/py/plot_blast_pie.py"
params.blast_plot_interac     = "${projectDir}/bin/py/interactive_plot_blast_pie.py"
params.run_multiqc            = true
params.multiqc_config         = null

params.singularity_path = null

def bannerFile = new File("${projectDir}/name.txt")
if (bannerFile.exists()) {
    println bannerFile.text
}

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
      --sample_sheet <file>       Tab-separated file with columns: sample, strain, type, read1, read2
                                  strain: leave empty to use default "GRCm39"
                                  type values: "fastq", "bam", or "unmapped_fastq"
                                  For BAM files, read2 should be empty
                                  Example: samples.tsv
    
    INPUT/OUTPUT:
      --outdir <path>             Output directory (default: results)
      --fastq_pattern <pattern>   FastQ file pattern (auto-detected if not specified)
                                  Examples: "*_{1,2}.fq.gz", "*_{R1,R2}.fastq.gz"
    
    REFERENCE FILES:
      --strains_base_dir <path>        Directory with strain-specific pseudogenomes (required for STAR)
      --standard_references_dir <path> Directory with standard reference genomes (e.g., GRCm39, GRCm38)
      --star_index_dir <path>          Directory for STAR indices (required for STAR)
    
    STAR ALIGNMENT:
      --run_star_alignment        Enable STAR alignment (default: true)
      --strain <strain(s)>        Comma-separated strain list (alternative to sample sheet)
                                  Example: "C57BL_6J,DBA_2J"
    
    SUBSAMPLING OPTIONS (Performance Tuning):
      --subset_for_fastq_qc       Subsample FASTQs for QC (default: true)
      --subset_fastq_qc_reads     Number of reads for FastQ QC (default: 100000)
      
      --subset_for_star           Subsample FASTQs for STAR alignment (default: false)
      --subset_star_reads         Number of reads for STAR (default: 200000)
      
      --subset_bam_for_qc         Subsample BAMs for QC tools (default: true)
      --bam_qc_subset_mapped      Number of mapped reads for BAM QC (default: 2000000)
      
      --subset_mapped_for_blast   Subsample mapped reads for BLAST (default: true)
      --mapped_subset_reads       Number of mapped reads for BLAST (default: 100000)
      
      --subset_unmapped_for_blast Subsample unmapped for BLAST (default: true)
      --unmapped_subset_reads     Number of unmapped reads for BLAST (default: 100000)
      
      --subset_seed <int>         Random seed for reproducibility (default: 100)
    
    QC TOOLS:
      --run_fastqc                Run FastQC on input FASTQs (default: true)
      --run_fastq_screen          Run FastQ Screen (default: true)
      --fastq_screen_conf         FastQ Screen config file (default: provided)
      
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
      --qualimap_protocol         Protocol type (default: "strand-specific-reverse")
    
    CONTAMINATION ANALYSIS:
      --run_contamination_check   Enable contamination checking (default: true)
      --contamination_blast_dbs   Path to BLAST database(s) (required for contamination check)
      --contamination_evalue      BLAST e-value threshold (default: "1e-10")
      --contamination_min_reads   Minimum reads to report (default: 100)
      --contamination_gc_min      GC content minimum (default: 0.60)
      --contamination_gc_max      GC content maximum (default: 1.0)
      --contamination_mapq        MAPQ threshold (default: 10)
      
    REPORTING:
      --run_multiqc               Run MultiQC (default: true)
      --run_blast_plots           Generate BLAST plots (default: true)
      --multiqc_config            MultiQC config file (default: null)
    
    DEFAULT RUN (Recommended Settings):
       
       nextflow run main_full.nf \\
       --sample_sheet input_test/sample_sheet.tsv \\
       --outdir results_test
       
       This runs with defaults:
       - FastQ QC: 100k reads
       - STAR: full FASTQs
       - BAM QC: 200k reads
       - BLAST: 100k reads
       - All QC tools enabled

       Note: sample sheet must include strain column
       References are automatically prepared per-strain from the strains_base_dir
    
    ========================================
    MORE INFORMATION
    ========================================
    
    Project: Conterminator v1.3
    Author:  Alaa Badreddine
    Home:    https://github.com/Z-Zen/Conterminator
    
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
// Parameter Validation
// ====================
def errors = []

if (params.outdir == null) {
    errors << "Please specify the parameter '--outdir [output folder]'"
}

if (params.run_star_alignment) {
    if (!params.strains_base_dir && !params.standard_references_dir) {
        errors << "STAR alignment requires at least one of:\n  --strains_base_dir <path>  (strain-specific pseudogenomes)\n  --standard_references_dir <path>  (standard reference genomes)"
    }
    if (!params.star_index_dir) {
        errors << "STAR alignment requires --star_index_dir <path> (directory for STAR indices)"
    }
}

if (params.run_contamination_check && !params.contamination_blast_dbs) {
    errors << "Contamination check requires --contamination_blast_dbs <path> (BLAST database directory)"
}

if (params.run_fastq_screen && !params.fastq_screen_conf) {
    errors << "FastQ Screen requires --fastq_screen_conf <path> (FastQ Screen config file)"
}

if (!errors.isEmpty()) {
    println """
    ========================================
    ERROR: Missing required parameters
    ========================================
    ${errors.join('\n    ')}

    Run 'nextflow run main.nf --help' for more information
    ========================================
    """
    exit 1
}

if (workflow.profile.contains('singularity')) {
    if (params.singularity_path && !file(params.singularity_path).exists()) {
        println "ERROR: Singularity image '${params.singularity_path}' does not exist."
        println "Please provide a valid path with '--singularity_path /full/path/to/conterminator.sif'"
        println "Or check ${workflow.manifest.homePage} for instructions on how to build the image."
        exit 1
    }
}

///////////////////////////////////////////
// Copy Nextflow PID to output directory //
///////////////////////////////////////////

def nextflowPid = new File(System.getProperty("user.dir") + "/.nextflow.pid")
def target = new File(params.outdir + "/pid.txt")

// Ensure output directory exists
def outputDir = new File(params.outdir)
if (!outputDir.exists()) {
    outputDir.mkdirs()
}

// Read PID if file exists
def pidContent = nextflowPid.exists() ? nextflowPid.text : "PID file not found\n"

// Build pid content string
def pidString = "---------------------\n" +
            pidContent +
            "${workflow.commandLine}\n" +
            "Start: ${workflow.start.format('dd-MMM-yyyy HH:mm:ss')}\n"

// Write or append to target
if (target.exists()) {
    target.append(pidString)
} else {
    target.write(pidString)
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
            println "Auto-detected FastQ pattern: ${pattern}"
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
    
    // Read file and normalize line endings (handle both LF and CRLF)
    def content = sheetFile.text.replaceAll('\r\n', '\n').replaceAll('\r', '\n')
    def lines = content.split('\n').findAll { it.trim() }  // Split and remove empty lines
    
    if (lines.size() < 2) {
        println """
        ========================================
        ERROR: Sample sheet is empty or has no data rows
        ========================================
        The sample sheet must contain a header row and at least one data row.
        ========================================
        """
        System.exit(1)
    }
    
    def header = lines[0].split('\t')

    // Check for 5 columns: sample, strain, type, read1, read2
    if (header.size() < 5 || header[0].toLowerCase() != 'sample' ||
        header[1].toLowerCase() != 'strain' || header[2].toLowerCase() != 'type' ||
        header[3].toLowerCase() != 'read1' || header[4].toLowerCase() != 'read2') {
        println """
        ========================================
        ERROR: Invalid sample sheet format
        ========================================
        Expected header: sample\tstrain\ttype\tread1\tread2
        Found: ${header.join('\t')}

        Example format:
        sample\tstrain\ttype\tread1\tread2
        43\tC57BL_6J\tfastq\t/path/to/sample_R1.fq.gz\t/path/to/sample_R2.fq.gz
        44\tC57BL_6J\tbam\t/path/to/sample.bam\t
        45\tDBA_2J\tunmapped_fastq\t/path/to/unmapped_R1.fq.gz\t/path/to/unmapped_R2.fq.gz
        46\t\tfastq\t/path/to/sample_R1.fq.gz\t/path/to/sample_R2.fq.gz

        Valid type values: fastq, bam, unmapped_fastq
        For BAM files, read2 should be empty
        Leave strain column empty to use default: GRCm39
        ========================================
        """
        System.exit(1)
    }

    // Store strain, type, and file paths
    lines[1..-1].each { line ->
        // Skip empty lines and comments
        def trimmedLine = line.trim()
        if (trimmedLine && !trimmedLine.startsWith('#')) {
            def fields = line.split('\t')
            if (fields.size() >= 4) {
                def sample = fields[0].trim()
                def strain = fields[1].trim()
                def type = fields[2].trim().toLowerCase()
                def read1 = fields[3].trim()
                def read2 = fields.size() >= 5 ? fields[4].trim() : ''

                // Set default strain to GRCm39 if empty
                if (!strain || strain.isEmpty()) {
                    strain = 'GRCm39'
                    println "INFO: Sample ${sample} has no strain specified, using default: GRCm39"
                }

                // Validate type
                if (!['fastq', 'bam', 'unmapped_fastq'].contains(type)) {
                    println "ERROR: Invalid type '${type}' for sample ${sample}. Must be: fastq, bam, or unmapped_fastq"
                    System.exit(1)
                }

                // Validate BAM files don't have read2
                if (type == 'bam' && read2) {
                    println "WARNING: Sample ${sample} is type 'bam' but has read2 specified. Ignoring read2."
                    read2 = ''
                }

                // Validate fastq/unmapped_fastq have both reads
                if ((type == 'fastq' || type == 'unmapped_fastq') && !read2) {
                    println "ERROR: Sample ${sample} is type '${type}' but missing read2"
                    System.exit(1)
                }

                sampleSheet[sample] = [strain: strain, type: type, read1: read1, read2: read2]
            }
        }
    }
    
    if (sampleSheet.isEmpty()) {
        println """
        ========================================
        ERROR: No valid samples found in sample sheet
        ========================================
        Please check that your sample sheet has at least one valid data row.
        ========================================
        """
        System.exit(1)
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
}

if (params.run_star_alignment) {
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

    if (!params.sample_sheet && !params.strain) {
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

    if (params.strain) {
        requestedStrains = parseStrains()
        println """
        ========================================
        Strain Configuration
        ========================================
        Requested strains: ${requestedStrains.join(', ')}
        Mode: All samples will be aligned to ALL specified strains
        ========================================
        """
    }
    
    // Validate that each requested strain has a directory in either strains_base_dir or standard_references_dir
    def invalidStrains = requestedStrains.findAll { strain ->
        def pseudogenomeDir = new File("${params.strains_base_dir}/${strain}")
        def standardRefDir = new File("${params.standard_references_dir}/${strain}")

        boolean pseudogenomeExists = pseudogenomeDir.exists() && pseudogenomeDir.isDirectory()
        boolean standardRefExists = standardRefDir.exists() && standardRefDir.isDirectory()

        return !pseudogenomeExists && !standardRefExists
    }

    if (invalidStrains) {
        println """
        ========================================
        ERROR: Invalid strain(s) specified
        ========================================
        The following strain(s) do not have directories in either location:
        ${invalidStrains.join(', ')}

        Searched in:
          1. Strain-specific pseudogenomes: ${params.strains_base_dir}
          2. Standard reference genomes: ${params.standard_references_dir}

        Available strains in pseudogenomes directory (${availableStrains.size()}):
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

// Define tool paths based on execution profile
// When using Singularity, tools are on $PATH inside the container — use bare command names.
// When running locally, use explicit paths from params (user must provide them).
def usingSingularity = workflow.profile.contains('singularity')

def toolPath = { paramValue, toolName ->
    if (paramValue) return paramValue          // explicit param always wins
    return toolName                            // rely on $PATH (container or local)
}

// Core alignment tools
def STAR_BIN        = toolPath(params.star_bin, "STAR")
def SAMTOOLS_BIN    = toolPath(params.samtools_bin, "samtools")
def SEQTK_BIN       = toolPath(params.seqtk_bin, "seqtk")

// QC tools
def FASTQC_BIN       = toolPath(params.fastqc_bin, "fastqc")
def FASTQ_SCREEN_BIN = toolPath(params.fastq_screen_bin, "fastq_screen")
def FASTQ_SCREEN_CONF = params.fastq_screen_conf
def QUALIMAP_BIN     = toolPath(params.qualimap_bin, "qualimap")

// Analysis tools
def BEDTOOLS_BIN     = toolPath(params.bedtools_bin, "bedtools")
// Picard: if user provides a .jar path, prefix with "java -jar"; otherwise use bare command
def PICARD_CMD = params.picard_jar ? (params.picard_jar.endsWith('.jar') ? "java -jar ${params.picard_jar}" : params.picard_jar) : "picard"
def DEEPTOOLS_GCBIAS = toolPath(params.deeptools_gcbias, "computeGCBias")
def BLASTN_BIN       = toolPath(params.blastn_bin, "blastn")

// UCSC tools
def FACOUNT_BIN      = toolPath(params.facount_bin, "faCount")
def FA2BIT_BIN       = toolPath(params.fa2bit_bin, "faToTwoBit")

// Other tools
def REFORMAT_BIN     = toolPath(params.reformat_bin, "reformat.sh")
def MULTIQC_BIN      = toolPath(null, "multiqc")
def MAPINSIGHTS_BIN  = toolPath(params.mapinsights_bin, "mapinsights")

if (usingSingularity) {
    println "INFO: Using Singularity - tools resolved via container \$PATH"
}

// Print configuration
println """
========================================
FLEXIBLE SUBSAMPLING CONFIGURATION
========================================
FastQ QC Subsampling:    ${params.subset_for_fastq_qc} (${params.subset_fastq_qc_reads} reads)
STAR Subsampling:        ${params.subset_for_star} (${params.subset_star_reads} reads)
BAM QC Subsampling:      ${params.subset_bam_for_qc} (${params.bam_qc_subset_mapped} reads)
Mapped BLAST Subset:     ${params.subset_mapped_for_blast} (${params.mapped_subset_reads} reads)
Unmapped BLAST Subset:   ${params.subset_unmapped_for_blast} (${params.unmapped_subset_reads} reads)
========================================
"""

// ====================
// Processes
// ====================

// Write pipeline manifest info to output directory
process WRITE_PIPELINE_INFO {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path main_script, stageAs: "main_script_copy.nf"
    path nextflow_config, stageAs: "nextflow_config_copy.config"
    path singularity_def, stageAs: "singularity_def_copy.def"

    output:
    path "pipeline_info.txt"
    path "main.nf"
    path "nextflow.config"
    path "conterminator.def"

    script:
    """
    cat > pipeline_info.txt << EOF
========================================
CONTERMINATOR PIPELINE INFORMATION
========================================

Pipeline:     ${workflow.manifest.description}
Version:      ${workflow.manifest.version}
Author:       ${workflow.manifest.author}
Homepage:     ${workflow.manifest.homePage}
Main Script:  ${workflow.manifest.mainScript}

Run Information:
----------------
Run Name:     ${workflow.runName}
Session ID:   ${workflow.sessionId}
Started:      ${workflow.start}
Profile:      ${workflow.profile}
${workflow.profile.contains('singularity') ? "Singularity:  ${params.singularity_path ?: 'enabled'}" : ''}
Work Dir:     ${workflow.workDir}
Output Dir:   ${params.outdir}
Command Line: ${workflow.commandLine}

Nextflow:
---------
Version:      ${workflow.nextflow.version}
Build:        ${workflow.nextflow.build}

Parameters:
-----------
Sample Sheet: ${params.sample_sheet}
Strains Dir:  ${params.strains_base_dir}
STAR Index:   ${params.star_index_dir}

QC Tools:
---------
FastQC:       ${params.run_fastqc}
FastQ Screen: ${params.run_fastq_screen}
STAR Align:   ${params.run_star_alignment}
Deeptools:    ${params.run_deeptools}
Picard GC:    ${params.run_picard_gc}
BEDtools GC:  ${params.run_bedtools_gc}
MapInsights:  ${params.run_mapinsights}
Qualimap:     ${params.run_qualimap}

Contamination:
--------------
Check:        ${params.run_contamination_check}
BLAST DBs:    ${params.contamination_blast_dbs ?: 'Not specified'}

Subsampling:
------------
FastQ QC:     ${params.subset_for_fastq_qc} (${params.subset_fastq_qc_reads} reads)
STAR:         ${params.subset_for_star} (${params.subset_star_reads} reads)
BAM QC:       ${params.subset_bam_for_qc} (${params.bam_qc_subset_mapped} reads)
Mapped BLAST: ${params.subset_mapped_for_blast} (${params.mapped_subset_reads} reads)
Unmapped BLAST: ${params.subset_unmapped_for_blast} (${params.unmapped_subset_reads} reads)

Pipeline Files:
---------------
Copies of main.nf, nextflow.config, and conterminator.def are saved
alongside this file for reproducibility.

========================================
EOF

    # Copy pipeline configuration files
    cp main_script_copy.nf main.nf
    cp nextflow_config_copy.config nextflow.config
    cp singularity_def_copy.def conterminator.def
    """
}

// Simple process to copy sample sheet to output
process COPY_SAMPLE_SHEET {
    publishDir "${params.outdir}/Input/", mode: 'copy'

    input:
    path sample_sheet

    output:
    path sample_sheet

    script:
    """
    # Just pass through - publishDir will copy it
    """
}

process SUBSET_FASTQ_FOR_QC {
    tag "${sample}"
    publishDir "${params.outdir}/Input/subsampled_fastq/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fastq1), path(fastq2)

    output:
    tuple val(sample), path("${sample}_subsampled_1.fastq.gz"), path("${sample}_subsampled_2.fastq.gz"), emit: subsampled_fastq
    path "${sample}_subsample_stats.txt", emit: stats

    when:
    params.subset_for_fastq_qc

    script:
    """
    set -euo pipefail
    
    echo "Subsampling FastQ for QC: ${sample}" > ${sample}_subsample_stats.txt
    echo "======================================" >> ${sample}_subsample_stats.txt
    echo "Target: ${params.subset_fastq_qc_reads} reads per file" >> ${sample}_subsample_stats.txt
    echo "Random seed: ${params.subset_seed}" >> ${sample}_subsample_stats.txt
    echo "" >> ${sample}_subsample_stats.txt
    
    ORIG_READS_R1=\$(zcat ${fastq1} | wc -l | awk '{print \$1/4}')
    ORIG_READS_R2=\$(zcat ${fastq2} | wc -l | awk '{print \$1/4}')
    
    echo "Original reads R1: \$ORIG_READS_R1" >> ${sample}_subsample_stats.txt
    echo "Original reads R2: \$ORIG_READS_R2" >> ${sample}_subsample_stats.txt
    
    ${SEQTK_BIN} sample -s${params.subset_seed} ${fastq1} ${params.subset_fastq_qc_reads} | gzip > ${sample}_subsampled_1.fastq.gz
    ${SEQTK_BIN} sample -s${params.subset_seed} ${fastq2} ${params.subset_fastq_qc_reads} | gzip > ${sample}_subsampled_2.fastq.gz
    
    SUB_READS_R1=\$(zcat ${sample}_subsampled_1.fastq.gz | wc -l | awk '{print \$1/4}')
    SUB_READS_R2=\$(zcat ${sample}_subsampled_2.fastq.gz | wc -l | awk '{print \$1/4}')
    
    echo "Subsampled reads R1: \$SUB_READS_R1" >> ${sample}_subsample_stats.txt
    echo "Subsampled reads R2: \$SUB_READS_R2" >> ${sample}_subsample_stats.txt
    echo "Purpose: FastQ Screen and FastQ-based QC" >> ${sample}_subsample_stats.txt
    """
}

process SUBSET_FASTQ_FOR_STAR {
    tag "${sample}"
    publishDir "${params.outdir}/Input/subsampled_fastq_star/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fastq1), path(fastq2)

    output:
    tuple val(sample), path("${sample}_star_subsampled_1.fastq.gz"), path("${sample}_star_subsampled_2.fastq.gz"), emit: subsampled_fastq
    path "${sample}_star_subsample_stats.txt", emit: stats

    when:
    params.subset_for_star

    script:
    """
    set -euo pipefail
    
    echo "Subsampling FastQ for STAR: ${sample}" > ${sample}_star_subsample_stats.txt
    echo "======================================" >> ${sample}_star_subsample_stats.txt
    echo "Target: ${params.subset_star_reads} reads per file" >> ${sample}_star_subsample_stats.txt
    echo "Random seed: ${params.subset_seed}" >> ${sample}_star_subsample_stats.txt
    echo "" >> ${sample}_star_subsample_stats.txt
    
    ORIG_READS_R1=\$(zcat ${fastq1} | wc -l | awk '{print \$1/4}')
    ORIG_READS_R2=\$(zcat ${fastq2} | wc -l | awk '{print \$1/4}')
    
    echo "Original reads R1: \$ORIG_READS_R1" >> ${sample}_star_subsample_stats.txt
    echo "Original reads R2: \$ORIG_READS_R2" >> ${sample}_star_subsample_stats.txt
    
    ${SEQTK_BIN} sample -s${params.subset_seed} ${fastq1} ${params.subset_star_reads} | gzip > ${sample}_star_subsampled_1.fastq.gz
    ${SEQTK_BIN} sample -s${params.subset_seed} ${fastq2} ${params.subset_star_reads} | gzip > ${sample}_star_subsampled_2.fastq.gz
    
    SUB_READS_R1=\$(zcat ${sample}_star_subsampled_1.fastq.gz | wc -l | awk '{print \$1/4}')
    SUB_READS_R2=\$(zcat ${sample}_star_subsampled_2.fastq.gz | wc -l | awk '{print \$1/4}')
    
    echo "Subsampled reads R1: \$SUB_READS_R1" >> ${sample}_star_subsample_stats.txt
    echo "Subsampled reads R2: \$SUB_READS_R2" >> ${sample}_star_subsample_stats.txt
    echo "Purpose: STAR alignment" >> ${sample}_star_subsample_stats.txt
    """
}

process PREPARE_STRAIN_REFERENCE {
    tag "${strain}"
    publishDir "${params.star_index_dir}", mode: 'copy', pattern: "${strain}"

    input:
    val strain

    output:
    tuple val(strain), path("${strain}_genome.fa{,.gz}"), path("${strain}_annotation.gtf{,.gz}"), emit: references

    script:
    """
    set -euo pipefail

    # Check both possible locations for the strain
    PSEUDOGENOME_DIR="${params.strains_base_dir}/${strain}"
    STANDARD_REF_DIR="${params.standard_references_dir}/${strain}"

    # For standard reference genomes like GRCm39, GRCm38, etc., check standard_references_dir FIRST
    # For strain-specific pseudogenomes, check strains_base_dir first
    case "${strain}" in
        GRCm39|GRCm38|GRCh38|GRCh37|mm10|mm39|hg38|hg19)
            # This is a standard reference genome - prioritize standard_references_dir
            if [ -d "\${STANDARD_REF_DIR}" ]; then
                STRAIN_DIR="\${STANDARD_REF_DIR}"
                echo "INFO: Using standard reference directory: \${STRAIN_DIR}" >&2
            elif [ -d "\${PSEUDOGENOME_DIR}" ]; then
                STRAIN_DIR="\${PSEUDOGENOME_DIR}"
                echo "INFO: Using strain-specific pseudogenome directory: \${STRAIN_DIR}" >&2
            else
                echo "ERROR: Could not find directory for ${strain}" >&2
                echo "Searched in:" >&2
                echo "  - \${STANDARD_REF_DIR} (standard reference)" >&2
                echo "  - \${PSEUDOGENOME_DIR} (pseudogenome)" >&2
                exit 1
            fi
            ;;
        *)
            # This is a strain-specific pseudogenome - prioritize strains_base_dir
            if [ -d "\${PSEUDOGENOME_DIR}" ]; then
                STRAIN_DIR="\${PSEUDOGENOME_DIR}"
                echo "INFO: Using strain-specific pseudogenome directory: \${STRAIN_DIR}" >&2
            elif [ -d "\${STANDARD_REF_DIR}" ]; then
                STRAIN_DIR="\${STANDARD_REF_DIR}"
                echo "INFO: Using standard reference directory: \${STRAIN_DIR}" >&2
            else
                echo "ERROR: Could not find directory for ${strain}" >&2
                echo "Searched in:" >&2
                echo "  - \${PSEUDOGENOME_DIR} (pseudogenome)" >&2
                echo "  - \${STANDARD_REF_DIR} (standard reference)" >&2
                exit 1
            fi
            ;;
    esac

    # First, try to find strain-specific pseudogenome files
    FASTA=\$(find "\${STRAIN_DIR}" -name "*pseudogenome__strain_${strain}.fa.gz" 2>/dev/null | head -1)

    # If not found, look for any FASTA file (for standard reference genomes like GRCm39)
    if [ -z "\${FASTA}" ]; then
        echo "INFO: Strain-specific pseudogenome not found, searching for standard reference genome..." >&2
        # Priority 1: Primary assembly or toplevel DNA (typical for Ensembl/UCSC)
        FASTA=\$(find "\${STRAIN_DIR}" -maxdepth 1 -type f \\( -name "*dna.primary_assembly.fa*" -o -name "*dna.toplevel.fa*" \\) 2>/dev/null | head -1)
        
        # Priority 2: Any DNA FASTA
        if [ -z "\${FASTA}" ]; then
            FASTA=\$(find "\${STRAIN_DIR}" -maxdepth 1 -type f -name "*dna.fa*" 2>/dev/null | head -1)
        fi
        
        # Priority 3: Any FASTA excluding cDNA/ncRNA/peptide (to avoid transcriptomes)
        if [ -z "\${FASTA}" ]; then
            FASTA=\$(find "\${STRAIN_DIR}" -maxdepth 1 -type f \\( -name "*.fa*" -o -name "*.fasta*" \\) ! -name "*.cdna.*" ! -name "*.ncrna.*" ! -name "*.pep.*" 2>/dev/null | head -1)
        fi

        # Fallback: Original behavior
        if [ -z "\${FASTA}" ]; then
            FASTA=\$(find "\${STRAIN_DIR}" -maxdepth 1 -type f \\( -name "*.fa.gz" -o -name "*.fasta.gz" -o -name "*.fa" -o -name "*.fasta" \\) 2>/dev/null | head -1)
        fi
    fi

    if [ -z "\${FASTA}" ]; then
        echo "ERROR: Could not find FASTA file for ${strain} in \${STRAIN_DIR}" >&2
        echo "Searched for:" >&2
        echo "  - *pseudogenome__strain_${strain}.fa.gz" >&2
        echo "  - *.dna.primary_assembly.fa*, *.dna.toplevel.fa*, *.dna.fa*" >&2
        echo "  - *.fa.gz, *.fasta.gz, *.fa, *.fasta (excluding cdna/ncrna/pep)" >&2
        exit 1
    fi

    # First, try to find strain-specific GTF
    GTF=\$(find "\${STRAIN_DIR}" -name "*pseudogenome__strain_${strain}.gtf.gz" 2>/dev/null | head -1)

    # If not found, look for any GTF file (for standard reference genomes)
    if [ -z "\${GTF}" ]; then
        echo "INFO: Strain-specific GTF not found, searching for standard annotation..." >&2
        # Priority 1: Standard GTF (avoiding abinitio, etc.)
        GTF=\$(find "\${STRAIN_DIR}" -maxdepth 1 -type f \\( -name "*.gtf.gz" -o -name "*.gtf" \\) ! -name "*.abinitio.*" 2>/dev/null | head -1)
        
        # Fallback: Original behavior
        if [ -z "\${GTF}" ]; then
            GTF=\$(find "\${STRAIN_DIR}" -maxdepth 1 -type f \\( -name "*.gtf.gz" -o -name "*.gtf" -o -name "*.gff.gz" -o -name "*.gff3.gz" -o -name "*.gff" -o -name "*.gff3" \\) 2>/dev/null | head -1)
        fi
    fi

    if [ -z "\${GTF}" ]; then
        echo "ERROR: Could not find GTF/GFF file for ${strain} in \${STRAIN_DIR}" >&2
        echo "Searched for:" >&2
        echo "  - *pseudogenome__strain_${strain}.gtf.gz" >&2
        echo "  - *.gtf.gz, *.gtf, *.gff.gz, *.gff3.gz, *.gff, *.gff3" >&2
        exit 1
    fi

    # Link with appropriate extension
    if [[ "\${FASTA}" == *.gz ]]; then
        ln -s "\${FASTA}" ${strain}_genome.fa.gz
    else
        ln -s "\${FASTA}" ${strain}_genome.fa
    fi

    if [[ "\${GTF}" == *.gz ]]; then
        ln -s "\${GTF}" ${strain}_annotation.gtf.gz
    else
        ln -s "\${GTF}" ${strain}_annotation.gtf
    fi

    echo "Prepared references for ${strain}"
    echo "FASTA: \${FASTA}"
    echo "GTF: \${GTF}"
    """
}

process BUILD_STAR_INDEX {
    tag "${strain}"
    publishDir "${params.star_index_dir}", mode: 'copy', pattern: "${strain}"

    input:
    tuple val(strain), path(fasta), path(gtf)

    output:
    tuple val(strain), path("${strain}"), emit: star_index

    script:
    """
    set -euo pipefail
    
    echo "Building STAR index for strain: ${strain}"
    echo "This may take 30-60 minutes..."
    
    mkdir -p ${strain}
    
    if [[ ${fasta} == *.gz ]]; then
        zcat ${fasta} > ${strain}_genome.fa
        FASTA_FILE="${strain}_genome.fa"
    else
        FASTA_FILE="${fasta}"
    fi
    
    if [[ ${gtf} == *.gz ]]; then
        zcat ${gtf} > ${strain}_annotation.gtf
        GTF_FILE="${strain}_annotation.gtf"
    else
        GTF_FILE="${gtf}"
    fi
    
    ${STAR_BIN} \
        --runMode genomeGenerate \
        --runThreadN ${task.cpus} \
        --genomeDir ${strain} \
        --genomeFastaFiles \${FASTA_FILE} \
        --sjdbGTFfile \${GTF_FILE} \
        --sjdbOverhang 149

    echo "STAR index built successfully for ${strain}"
    """
}

process STAR_ALIGN {
    tag "${sample}_${strain}"
    publishDir "${params.outdir}/Temporary/star_alignment/${sample}/${strain}", mode: 'copy'

    input:
    tuple val(sample), path(fastq1), path(fastq2), val(strain), path(genome_dir)

    output:
    tuple val(sample), val(strain), path("${sample}_${strain}_Aligned.sortedByCoord.out.bam"), path("${sample}_${strain}_Aligned.sortedByCoord.out.bam.bai"), emit: aligned_bam
    tuple val(sample), val(strain), path("${sample}_${strain}_Unmapped.out.mate1"), path("${sample}_${strain}_Unmapped.out.mate2"), emit: unmapped_fastq
    path "${sample}_${strain}_Log.final.out", emit: log
    path "${sample}_${strain}_Log.out", emit: log_out
    path "${sample}_${strain}_Log.progress.out", emit: log_progress
    path "${sample}_${strain}_ReadsPerGene.out.tab", optional: true, emit: gene_counts
    path "${sample}_${strain}_Aligned.toTranscriptome.out.bam", optional: true, emit: transcriptome_bam
    path "${sample}_${strain}_star_alignment_stats.txt", emit: stats

    when:
    params.run_star_alignment

    script:
    """
    set -euo pipefail
    
    echo "STAR Alignment for ${sample} against ${strain}" > ${sample}_${strain}_star_alignment_stats.txt
    echo "===========================================" >> ${sample}_${strain}_star_alignment_stats.txt
    echo "Genome directory: ${genome_dir}" >> ${sample}_${strain}_star_alignment_stats.txt
    echo "Threads: ${task.cpus}" >> ${sample}_${strain}_star_alignment_stats.txt
    echo "" >> ${sample}_${strain}_star_alignment_stats.txt
    
    ${STAR_BIN} \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${genome_dir} \\
        --quantMode GeneCounts TranscriptomeSAM \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesIn ${fastq1} ${fastq2} \\
        --outFileNamePrefix ${sample}_${strain}_ \\
        --outReadsUnmapped Fastx
    
    ${SAMTOOLS_BIN} index -@ ${task.cpus} ${sample}_${strain}_Aligned.sortedByCoord.out.bam

    echo "" >> ${sample}_${strain}_star_alignment_stats.txt
    echo "ALIGNMENT SUMMARY:" >> ${sample}_${strain}_star_alignment_stats.txt
    cat ${sample}_${strain}_Log.final.out >> ${sample}_${strain}_star_alignment_stats.txt
    """
}

// Reference preparation for QC tools
process PREP_STRAIN_REFERENCE_FOR_QC {
    tag "${strain}"
    publishDir "${params.outdir}/References/${strain}", mode: 'copy'
    
    input:
    tuple val(strain), path(fasta), path(gtf)
    
    output:
    tuple val(strain), 
          path("${strain}.fa"), 
          path("${strain}.fa.fai"), 
          path("${strain}.dict"), 
          path("${strain}.bed"), 
          path("${strain}.txt"), 
          path("${strain}.2bit"),
          path("${strain}.gtf.gz"), emit: qc_references
    
    script:
    """
    set -euo pipefail
    
    # Decompress FASTA if needed (samtools can't index gzipped files with regular gzip)
    echo "Processing FASTA file: ${fasta}"
    if [[ ${fasta} == *.gz ]]; then
        echo "Decompressing FASTA..."
        gunzip -c ${fasta} > ${strain}.fa
    else
        echo "Linking uncompressed FASTA..."
        ln -s ${fasta} ${strain}.fa
    fi
    
    # Index FASTA
    echo "Indexing FASTA..."
    ${SAMTOOLS_BIN} faidx -@ ${task.cpus} ${strain}.fa
    
    # Create sequence dictionary
    echo "Creating sequence dictionary..."
    ${SAMTOOLS_BIN} dict ${strain}.fa > ${strain}.dict
    
    # Create BED file from FAI
    echo "Creating BED file..."
    awk '{print \$1"\\t0\\t"\$2}' ${strain}.fa.fai > ${strain}.bed
    
    # Calculate effective genome size
    echo "Calculating effective genome size..."
    ${FACOUNT_BIN} ${strain}.fa > ${strain}.faCount.txt
    awk 'NR>1 {total+=\$2; n+=\$3+\$4+\$5+\$6} END {print total-n}' ${strain}.faCount.txt > ${strain}.txt
    
    # Create 2bit file
    echo "Creating 2bit file..."
    ${FA2BIT_BIN} ${strain}.fa ${strain}.2bit
    
    # Handle GTF file (already staged as input)
    echo "Processing GTF file: ${gtf}"
    if [[ ${gtf} == *.gz ]]; then
        # Already compressed
        if [[ ${gtf} == ${strain}.gtf.gz ]]; then
            echo "GTF already has correct name"
        else
            echo "Copying compressed GTF..."
            cp ${gtf} ${strain}.gtf.gz
        fi
    else
        # Need to compress it
        echo "Compressing GTF..."
        gzip -c ${gtf} > ${strain}.gtf.gz
    fi
    
    echo "Reference preparation complete for ${strain}"
    """
}

process INDEX_INPUT_BAM {
    tag "${sample}_${strain}"
    publishDir "${params.outdir}/Input/indexed_bams/${sample}/${strain}", mode: 'copy'

    input:
    tuple val(sample), val(strain), path(bam)

    output:
    tuple val(sample), val(strain), path(bam), path("${bam}.bai"), emit: bam_with_index

    script:
    """
    ${SAMTOOLS_BIN} index -@ ${task.cpus} ${bam}
    """
}

process SUBSET_BAM_FOR_QC {
    tag "${sample}_${strain}"
    publishDir "${params.outdir}/Input/subsampled_bams/${sample}/${strain}", mode: 'copy'

    input:
    tuple val(sample), val(strain), path(bam), path(bai)

    output:
    tuple val(sample), val(strain), path("${sample}_${strain}_subset.bam"), path("${sample}_${strain}_subset.bam.bai"), emit: subsampled_bam
    path "${sample}_${strain}_bam_subset_stats.txt", emit: stats

    when:
    params.subset_bam_for_qc

    script:
    """
    set -euo pipefail
    
    echo "BAM Subsampling for QC: ${sample}_${strain}" > ${sample}_${strain}_bam_subset_stats.txt
    echo "==========================================" >> ${sample}_${strain}_bam_subset_stats.txt
    echo "Target: ${params.bam_qc_subset_mapped} mapped reads" >> ${sample}_${strain}_bam_subset_stats.txt
    echo "" >> ${sample}_${strain}_bam_subset_stats.txt
    
    # Count total mapped reads
    TOTAL_MAPPED=\$(${SAMTOOLS_BIN} view -@ ${task.cpus} -c -F 4 ${bam})
    echo "Total mapped reads: \$TOTAL_MAPPED" >> ${sample}_${strain}_bam_subset_stats.txt

    if [ "\$TOTAL_MAPPED" -gt "${params.bam_qc_subset_mapped}" ]; then
        FRACTION=\$(awk -v seed=${params.subset_seed} -v target=${params.bam_qc_subset_mapped} -v total=\$TOTAL_MAPPED 'BEGIN {frac=target/total*10000; printf "%d.%04d", seed, frac}')
        echo "Target reads: ${params.bam_qc_subset_mapped}" >> ${sample}_${strain}_bam_subset_stats.txt
        echo "Sampling fraction: \$FRACTION" >> ${sample}_${strain}_bam_subset_stats.txt
        
        # Subsample mapped reads only
        ${SAMTOOLS_BIN} view -@ ${task.cpus} -b -s \${FRACTION} -F 4 ${bam} > temp_mapped.bam
        
        # Verify actual count
        ACTUAL_SAMPLED=\$(${SAMTOOLS_BIN} view -@ ${task.cpus} -c temp_mapped.bam)
        echo "Actually sampled: \$ACTUAL_SAMPLED reads" >> ${sample}_${strain}_bam_subset_stats.txt
    else
        echo "Using all mapped reads (total: \$TOTAL_MAPPED < target: ${params.bam_qc_subset_mapped})" >> ${sample}_${strain}_bam_subset_stats.txt
        ${SAMTOOLS_BIN} view -@ ${task.cpus} -b -F 4 ${bam} > temp_mapped.bam
    fi
    
    # Sort and index
    ${SAMTOOLS_BIN} sort -@ ${task.cpus} -o ${sample}_${strain}_subset.bam temp_mapped.bam
    ${SAMTOOLS_BIN} index -@ ${task.cpus} ${sample}_${strain}_subset.bam
    
    # Report final count
    SUBSET_COUNT=\$(${SAMTOOLS_BIN} view -@ ${task.cpus} -c ${sample}_${strain}_subset.bam)
    echo "Subsampled reads: \$SUBSET_COUNT" >> ${sample}_${strain}_bam_subset_stats.txt
    echo "Purpose: QC tools (DeepTools, Picard, BEDTools, Qualimap, Mapinsights)" >> ${sample}_${strain}_bam_subset_stats.txt
    
    rm -f temp_mapped.bam
    """
}

process SUBSET_UNMAPPED_FOR_BLAST {
    tag "${sample}_${strain}"
    publishDir "${params.outdir}/Input/subsampled_unmapped_blast/${sample}/${strain}", mode: 'copy'

    input:
    tuple val(sample), val(strain), path(unmapped_r1), path(unmapped_r2)

    output:
    tuple val(sample), val(strain), path("${sample}_${strain}_unmapped_blast_R1.fastq.gz"), path("${sample}_${strain}_unmapped_blast_R2.fastq.gz"), emit: subsampled_fastq
    path "${sample}_${strain}_unmapped_blast_subset_stats.txt", emit: stats

    when:
    params.subset_unmapped_for_blast

    script:
    """
    set -euo pipefail
    
    echo "Subsampling Unmapped for BLAST: ${sample}_${strain}" > ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "===================================================" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "Target: ${params.unmapped_subset_reads} reads per file" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    
    TOTAL_R1=\$(cat ${unmapped_r1} | wc -l | awk '{print \$1/4}')
    TOTAL_R2=\$(cat ${unmapped_r2} | wc -l | awk '{print \$1/4}')
    
    echo "Total unmapped R1: \$TOTAL_R1" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "Total unmapped R2: \$TOTAL_R2" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    
    # Subsample R1
    if [ "\$TOTAL_R1" -gt "${params.unmapped_subset_reads}" ]; then
        ${SEQTK_BIN} sample -s${params.subset_seed} ${unmapped_r1} ${params.unmapped_subset_reads} | gzip > ${sample}_${strain}_unmapped_blast_R1.fastq.gz
        echo "R1 subsampled to ${params.unmapped_subset_reads} reads" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    else
        cat ${unmapped_r1} | gzip > ${sample}_${strain}_unmapped_blast_R1.fastq.gz
        echo "R1 using all \$TOTAL_R1 reads (less than target)" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    fi
    
    # Subsample R2
    if [ "\$TOTAL_R2" -gt "${params.unmapped_subset_reads}" ]; then
        ${SEQTK_BIN} sample -s${params.subset_seed} ${unmapped_r2} ${params.unmapped_subset_reads} | gzip > ${sample}_${strain}_unmapped_blast_R2.fastq.gz
        echo "R2 subsampled to ${params.unmapped_subset_reads} reads" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    else
        cat ${unmapped_r2} | gzip > ${sample}_${strain}_unmapped_blast_R2.fastq.gz
        echo "R2 using all \$TOTAL_R2 reads (less than target)" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    fi
    
    FINAL_R1=\$(zcat ${sample}_${strain}_unmapped_blast_R1.fastq.gz | wc -l | awk '{print \$1/4}')
    FINAL_R2=\$(zcat ${sample}_${strain}_unmapped_blast_R2.fastq.gz | wc -l | awk '{print \$1/4}')
    TOTAL_FINAL=\$(echo "\$FINAL_R1 + \$FINAL_R2" | bc)
    
    echo "Final R1 reads for BLAST: \$FINAL_R1" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "Final R2 reads for BLAST: \$FINAL_R2" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "Total final reads for BLAST: \$TOTAL_FINAL" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "Purpose: BLAST contamination screening" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    """
}

process DEEPTOOLS_GC {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/deeptools_gc_bias/${sample}/", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai), val(bam_type)
    path ref2bit
    path effective_size_file

    output:
    path "${sample}.gcBias.freq.txt", emit: freq
    path "${sample}.gcBias.plot.pdf", emit: plot

    when:
    params.run_deeptools

    script:
    def blacklist = params.blacklist_bed ? "--blackListFileName ${params.blacklist_bed}" : ""
    """
    set -euo pipefail
    
    EFFECTIVE_SIZE=\$(cat ${effective_size_file})
    
    ${DEEPTOOLS_GCBIAS} \\
        --bamfile ${bam} \\
        --genome ${ref2bit} \\
        --effectiveGenomeSize \$EFFECTIVE_SIZE \\
        --GCbiasFrequenciesFile ${sample}.gcBias.freq.txt \\
        --biasPlot ${sample}.gcBias.plot.pdf \\
        --numberOfProcessors ${task.cpus} \\
        ${blacklist}
    """
}

process PICARD_GC_BIAS {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/picard_gc_bias/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai), val(bam_type)
    path ref

    output:
    path "${sample}.gc_bias_metrics.txt", emit: metrics
    path "${sample}.gc_bias_summary.txt", emit: summary
    path "${sample}.gc_bias.pdf", emit: plot
    path "${sample}_allmetrics.*", emit: allmetrics

    when:
    params.run_picard_gc

    script:
    """
    set -euo pipefail
    export PATH="/opt/R/4.5.2/bin:\$PATH"

    ${PICARD_CMD} CollectGcBiasMetrics \\
        I=${bam} \\
        O=${sample}.gc_bias_metrics.txt \\
        SUMMARY_OUTPUT=${sample}.gc_bias_summary.txt \\
        CHART_OUTPUT=${sample}.gc_bias.pdf \\
        R=${ref} \\
        VALIDATION_STRINGENCY=LENIENT

    ${PICARD_CMD} CollectMultipleMetrics \\
        I=${bam} \\
        O=${sample}_allmetrics \\
        R=${ref}
    """
}

process BEDTOOLS_GC_COVERAGE {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/bedtools_gc/${sample}/", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai), val(bam_type)
    path ref
    path genome_bed

    output:
    path "${sample}.windows.bed.gz", emit: windows
    path "${sample}.windows.gc.tsv.gz", emit: gc
    path "${sample}.gc_coverage.tsv.gz", emit: coverage
    path "${sample}_gc_plots.png", emit: combined_plot, optional: true
    path "${sample}_gc_distribution.png", emit: gc_dist, optional: true
    path "${sample}_coverage_vs_gc.png", emit: cov_vs_gc, optional: true
    path "${sample}_coverage_distribution.png", emit: cov_dist, optional: true
    path "${sample}_gc_by_coverage_quartiles.png", emit: gc_quartiles, optional: true
    path "${sample}_gc_summary.txt", emit: summary, optional: true

    when:
    params.run_bedtools_gc

    script:
    """
    set -euo pipefail

    export PATH="/opt/R/4.5.2/bin:\$PATH"

    ${BEDTOOLS_BIN} makewindows -b ${genome_bed} -w ${params.windowsize} -s ${params.window_step} > ${sample}.windows.bed

    ${BEDTOOLS_BIN} nuc -fi ${ref} -bed ${sample}.windows.bed | \\
      awk 'BEGIN{OFS="\\t"} NR==1{print "chrom","start","end","pct_at","pct_gc"} NR>1{print \$1,\$2,\$3,\$4,\$5}' \\
      > ${sample}.windows.gc.tsv

    ${BEDTOOLS_BIN} coverage -a ${sample}.windows.bed -b ${bam} -mean | \\
      paste ${sample}.windows.gc.tsv - | \\
      awk 'BEGIN{OFS="\\t"} NR==1{print \$0,"mean_coverage"} NR>1{print}' \\
      > ${sample}.gc_coverage.tsv

    # Generate GC content plots
    Rscript ${params.gc_plot_script} ${sample}.gc_coverage.tsv ${sample}

    gzip ${sample}.gc_coverage.tsv ${sample}.windows.gc.tsv ${sample}.windows.bed
    """
}

process QUALIMAP_BAMQC {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/qualimap/${sample}/bamqc", mode: 'copy'

    input:
    tuple val(sample), val(strain), path(bam), path(bai), val(bam_type), path(gtf)  // GTF as input
    path ref

    output:
    path "qualimap_bamqc/**", emit: reports
    path "qualimap_bamqc/qualimapReport.html", optional: true, emit: html
    path "qualimap_bamqc/qualimapReport.pdf", optional: true, emit: pdf
    path "qualimap_bamqc/genome_results.txt", optional: true, emit: results

    when:
    params.run_qualimap && (params.qualimap_mode == "bamqc" || params.qualimap_mode == "both")

    script:
    """
    set -euo pipefail

    # Decompress GTF if needed
    if [[ ${gtf} == *.gz ]]; then
        zcat ${gtf} > ${strain}_annotation.gtf
        GTF_FILE="${strain}_annotation.gtf"
    else
        cp ${gtf} ${strain}_annotation.gtf
        GTF_FILE="${strain}_annotation.gtf"
    fi
    
    ${QUALIMAP_BIN} bamqc \\
        -bam ${bam} \\
        -c \\
        -gd ${params.qualimap_genome} \\
        -gff \${GTF_FILE} \\
        -outdir qualimap_bamqc \\
        -outformat PDF:HTML \\
        -p ${params.qualimap_protocol} \\
        -nt ${task.cpus}
    """
}

process QUALIMAP_RNASEQ {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/qualimap/${sample}/rnaseq", mode: 'copy'

    input:
    tuple val(sample), val(strain), path(bam), path(bai), val(bam_type), path(gtf)

    output:
    path "qualimap_rnaseq_unique_mapped_reads/**", emit: reports
    path "qualimap_rnaseq_unique_mapped_reads/qualimapReport.html", optional: true, emit: html
    path "qualimap_rnaseq_unique_mapped_reads/qualimapReport.pdf", optional: true, emit: pdf
    path "qualimap_rnaseq_unique_mapped_reads/rnaseq_results.txt", optional: true, emit: results
    path "qualimap_rnaseq_proportional/**", emit: reports_prop
    path "qualimap_rnaseq_proportional/qualimapReport.html", optional: true, emit: html_prop
    path "qualimap_rnaseq_proportional/qualimapReport.pdf", optional: true, emit: pdf_prop
    path "qualimap_rnaseq_proportional/rnaseq_results.txt", optional: true, emit: results_prop

    when:
    params.run_qualimap && (params.qualimap_mode == "rnaseq" || params.qualimap_mode == "both")

    script:
    """
    set -euo pipefail

    if [[ ${gtf} == *.gz ]]; then
        zcat ${gtf} > ${strain}_annotation.gtf
        GTF_FILE="${strain}_annotation.gtf"
    else
        GTF_FILE="${gtf}"
    fi
    
    ${QUALIMAP_BIN} rnaseq \\
        --algorithm uniquely-mapped-reads \\
        -bam ${bam} \\
        -gtf \${GTF_FILE} \\
        -outdir qualimap_rnaseq_unique_mapped_reads \\
        -outformat PDF:HTML \\
        -p ${params.qualimap_protocol} \\
        -pe \\
        -s

    ${QUALIMAP_BIN} rnaseq \\
        --algorithm proportional \\
        -bam ${bam} \\
        -gtf \${GTF_FILE} \\
        -outdir qualimap_rnaseq_proportional \\
        -outformat PDF:HTML \\
        -p ${params.qualimap_protocol} \\
        -pe \\
        -s 
    """
}

process MAPINSIGHTS {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/mapinsights/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai), val(bam_type)
    path ref

    output:
    path "${sample}_mapinsights/**", emit: reports
    path "${sample}_mapinsights/report.html", optional: true, emit: html

    when:
    params.run_mapinsights

    script:
    """
    set -euo pipefail
    
    export PATH="/opt/R/4.5.2/bin:\$PATH"

    mkdir ${sample}_mapinsights
    
    if [ -f "${MAPINSIGHTS_BIN}" ] && [ -x "${MAPINSIGHTS_BIN}" ]; then
        ${MAPINSIGHTS_BIN} bamqc \\
            -r ${ref} \\
            -i ${bam} \\
            -o ${sample}_mapinsights \\
            ${params.mapinsights_opts}
    else
        mkdir -p ${sample}_mapinsights
        ${SAMTOOLS_BIN} flagstat -@ ${task.cpus} ${bam} > ${sample}_mapinsights/${sample}.flagstat
        ${SAMTOOLS_BIN} stats -@ ${task.cpus} ${bam} > ${sample}_mapinsights/${sample}.stats
        ${SAMTOOLS_BIN} idxstats -@ ${task.cpus} ${bam} > ${sample}_mapinsights/${sample}.idxstats
    fi
    """
}

process FASTQ_SCREEN {
    tag "${sample}"
    publishDir "${params.outdir}/Output/fastq_screen/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fastq1), path(fastq2)

    output:
    tuple val(sample), path("*_screen.txt"), emit: summary
    path "*_screen.html", optional: true, emit: html
    path "*_screen.png", optional: true, emit: plot

    when:
    params.run_fastq_screen

    script:
    """
    set -euo pipefail
    ${FASTQ_SCREEN_BIN} \\
        --conf ${FASTQ_SCREEN_CONF} \\
        --threads ${task.cpus} \\
        --outdir . \\
        --force \\
        ${fastq1} ${fastq2}
    """
}

process FASTQC_INPUT_FASTQ {
    tag "${sample}"
    publishDir "${params.outdir}/Output/fastqc/${sample}/input_fastq", mode: 'copy'

    input:
    tuple val(sample), path(fastq1), path(fastq2)

    output:
    tuple path("*_fastqc.html"), path("*_fastqc.zip"), emit: results

    when:
    params.run_fastqc

    script:
    """
    set -euo pipefail
    
    ${FASTQC_BIN} \\
        --outdir . \\
        --nogroup \\
        --threads ${task.cpus} \\
        ${fastq1} ${fastq2}
    """
}

process FASTQC_MAPPED_BAM {
    tag "${sample} ${strain}"
    publishDir "${params.outdir}/Output/fastqc/${sample}/mapped_bam", mode: 'copy'

    input:
    tuple val(sample), val(strain), path(bam), path(bai)

    output:
    tuple path("*_fastqc.html"), path("*_fastqc.zip"), emit: results

    when:
    params.run_fastqc && params.run_contamination_check

    script:
    """
    set -euo pipefail
    
    ${FASTQC_BIN} \\
        --outdir . \\
        --nogroup \\
        --threads ${task.cpus} \\
        --format bam \\
        ${bam}
    """
}

process FASTQC_UNMAPPED_FASTQ {
    tag "${sample} ${strain}"
    publishDir "${params.outdir}/Output/fastqc/${sample}/unmapped_fastq", mode: 'copy'

    input:
    tuple val(sample), val(strain), path(fastq_r1), path(fastq_r2)

    output:
    tuple path("*1_fastqc.html"), path("*1_fastqc.zip"), emit: results_r1
    tuple path("*2_fastqc.html"), path("*2_fastqc.zip"), emit: results_r2

    when:
    params.run_fastqc && params.run_contamination_check

    script:
    """
    set -euo pipefail
    
    # Run FastQC on R1
    ${FASTQC_BIN} \\
        --outdir . \\
        --nogroup \\
        --threads ${task.cpus} \\
        ${fastq_r1}
    
    # Run FastQC on R2
    ${FASTQC_BIN} \\
        --outdir . \\
        --nogroup \\
        --threads ${task.cpus} \\
        ${fastq_r2}
    """
}

process BLAST_MAPPED_READS_MULTI {
    tag "${sample}_${db_name}"
    publishDir "${params.outdir}/Output/contamination_check/${sample}/mapped/${db_name}", mode: 'copy'

    input:
    tuple val(sample), val(strain), path(bam), path(bai), val(db_name), val(db_path)

    output:
    tuple val(sample), val(db_name), path("${sample}_${db_name}.blast.tsv"), emit: blast_results
    tuple val(sample), val(db_name), path("${sample}_${db_name}.blast_summary.txt"), emit: summary
    tuple val(sample), val(db_name), path("${sample}_${db_name}.top_contaminants.txt"), emit: top_hits
    tuple val(sample), val(db_name), path("${sample}_${db_name}.blast.tsv"), path("${sample}_${db_name}.blast_summary.txt"), emit: for_plotting

    when:
    params.run_contamination_check

    script:
    """
    export BLASTDB=${params.contamination_blast_dbs}
    
    echo "Starting BLAST analysis: ${sample} vs ${db_name}" >&2
    echo "Extracting mapped reads from BAM: ${bam}" >&2
    
    # Extract mapped reads from BAM to FastQ
    ${SAMTOOLS_BIN} fastq \\
        -@ ${task.cpus} \\
        -1 ${sample}_mapped_R1.fastq \\
        -2 ${sample}_mapped_R2.fastq \\
        -0 /dev/null \\
        -s /dev/null \\
        -N \\
        ${bam}
    
    ${REFORMAT_BIN} in=${sample}_mapped_R1.fastq out=${sample}_R1.fasta
    ${REFORMAT_BIN} in=${sample}_mapped_R2.fastq out=${sample}_R2.fasta

    echo "Running BLAST on R1 against ${db_name}..." >&2
    # BLAST R1
    ${BLASTN_BIN} \\
        -query ${sample}_R1.fasta \\
        -db ${db_path} \\
        -out ${sample}_${db_name}_R1.blast.tsv \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
        -evalue ${params.contamination_evalue} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 5 \\
        -max_hsps 1
    
    echo "Running BLAST on R2 against ${db_name}..." >&2
    # BLAST R2
    ${BLASTN_BIN} \\
        -query ${sample}_R2.fasta \\
        -db ${db_path} \\
        -out ${sample}_${db_name}_R2.blast.tsv \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
        -evalue ${params.contamination_evalue} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 5 \\
        -max_hsps 1

    # Combine results
    cat ${sample}_${db_name}_R1.blast.tsv ${sample}_${db_name}_R2.blast.tsv > ${sample}_${db_name}.blast.tsv
    
    BLAST_HITS=\$(wc -l < ${sample}_${db_name}.blast.tsv)
    echo "Total BLAST hits (R1+R2): \$BLAST_HITS" >> ${sample}_${db_name}.blast_summary.txt
    
    if [ "\$BLAST_HITS" -gt 0 ]; then
        echo "" >> ${sample}_${db_name}.blast_summary.txt
        echo "TOP CONTAMINATING SEQUENCES/STRAINS:" >> ${sample}_${db_name}.blast_summary.txt
        
        if [[ "${db_name}" == *"mus_strain"* ]]; then
            cut -f2 ${sample}_${db_name}.blast.tsv | sed 's/gnl|//; s/|.*//' | sort | uniq -c | sort -rn | head -20 >> ${sample}_${db_name}.blast_summary.txt
            cut -f2 ${sample}_${db_name}.blast.tsv | sed 's/gnl|//; s/|.*//' | sort | uniq -c | sort -rn | head -50 > ${sample}_${db_name}.top_contaminants.txt
        else
            cut -f13 ${sample}_${db_name}.blast.tsv | sort | uniq -c | sort -rn | head -20 >> ${sample}_${db_name}.blast_summary.txt
            cut -f13 ${sample}_${db_name}.blast.tsv | sort | uniq -c | sort -rn | head -50 > ${sample}_${db_name}.top_contaminants.txt
        fi
    else
        echo "No contamination detected" > ${sample}_${db_name}.top_contaminants.txt
    fi
    
    echo "BLAST analysis complete for ${sample} vs ${db_name}" >&2

    rm -f ${sample}_R1.fasta ${sample}_R2.fasta ${sample}_mapped_R1.fastq ${sample}_mapped_R2.fastq
    """
}

process BLAST_UNMAPPED_READS_MULTI {
    tag "${sample}_${db_name}"
    publishDir "${params.outdir}/Output/contamination_check/${sample}/unmapped/${db_name}", mode: 'copy'

    input:
    tuple val(sample), path(r1_fastq), path(r2_fastq), val(db_name), val(db_path)

    output:
    tuple val(sample), val(db_name), path("${sample}_${db_name}.blast.tsv"), emit: blast_results
    tuple val(sample), val(db_name), path("${sample}_${db_name}.blast_summary.txt"), emit: summary
    tuple val(sample), val(db_name), path("${sample}_${db_name}.top_contaminants.txt"), emit: top_hits
    tuple val(sample), val(db_name), path("${sample}_${db_name}.blast.tsv"), path("${sample}_${db_name}.blast_summary.txt"), emit: for_plotting

    when:
    params.run_contamination_check

    script:
    """
    
    export BLASTDB=${params.contamination_blast_dbs}
    
    echo "Starting BLAST analysis: ${sample} vs ${db_name}" >&2
    
    # Count reads from both R1 and R2
    READ_COUNT_R1=\$(zcat ${r1_fastq} | wc -l | awk '{print \$1/4}')
    READ_COUNT_R2=\$(zcat ${r2_fastq} | wc -l | awk '{print \$1/4}')
    TOTAL_READ_COUNT=\$(echo "\$READ_COUNT_R1 + \$READ_COUNT_R2" | bc)
    
    echo "BLAST Contamination Analysis for ${sample} against ${db_name}" > ${sample}_${db_name}.blast_summary.txt
    echo "==========================================================" >> ${sample}_${db_name}.blast_summary.txt
    echo "Database: ${db_path}" >> ${sample}_${db_name}.blast_summary.txt
    echo "Unmapped R1 reads: \$READ_COUNT_R1" >> ${sample}_${db_name}.blast_summary.txt
    echo "Unmapped R2 reads: \$READ_COUNT_R2" >> ${sample}_${db_name}.blast_summary.txt
    echo "Total unmapped reads to analyze: \$TOTAL_READ_COUNT" >> ${sample}_${db_name}.blast_summary.txt
    
    # Subsample if needed (sample from EACH file to maintain pairing info)
    MAX_READS_PER_FILE=${params.unmapped_subset_reads}
    if [ "\$READ_COUNT_R1" -gt "\$MAX_READS_PER_FILE" ]; then
        ${SEQTK_BIN} sample -s${params.subset_seed} ${sample}_mapped_R1.fastq.gz ${params.mapped_subset_reads} | gzip > ${sample}_R1.sampled.fastq.gz
        BLAST_INPUT_R1="${sample}_R1.sampled.fastq.gz"
    else
        BLAST_INPUT_R1="${r1_fastq}"
    fi
    
    if [ "\$READ_COUNT_R2" -gt "\$MAX_READS_PER_FILE" ]; then
        ${SEQTK_BIN} sample -s${params.subset_seed} ${sample}_mapped_R2.fastq.gz ${params.mapped_subset_reads} | gzip > ${sample}_R2.sampled.fastq.gz
        BLAST_INPUT_R2="${sample}_R2.sampled.fastq.gz"
    else
        BLAST_INPUT_R2="${r2_fastq}"
    fi
    
    # Convert R1 to FASTA
    ${REFORMAT_BIN} in=\$BLAST_INPUT_R1 out=${sample}_R1.fasta
    
    # Convert R2 to FASTA
    ${REFORMAT_BIN} in=\$BLAST_INPUT_R2 out=${sample}_R2.fasta
    
    echo "Running BLAST on R1 against ${db_name}..." >&2
    # BLAST R1
    ${BLASTN_BIN} \\
        -query ${sample}_R1.fasta \\
        -db ${db_path} \\
        -out ${sample}_${db_name}_R1.blast.tsv \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
        -evalue ${params.contamination_evalue} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 5 \\
        -max_hsps 1
    
    echo "Running BLAST on R2 against ${db_name}..." >&2
    # BLAST R2
    ${BLASTN_BIN} \\
        -query ${sample}_R2.fasta \\
        -db ${db_path} \\
        -out ${sample}_${db_name}_R2.blast.tsv \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
        -evalue ${params.contamination_evalue} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 5 \\
        -max_hsps 1
    
    # Combine results
    cat ${sample}_${db_name}_R1.blast.tsv ${sample}_${db_name}_R2.blast.tsv > ${sample}_${db_name}.blast.tsv
    
    BLAST_HITS=\$(wc -l < ${sample}_${db_name}.blast.tsv)
    echo "Total BLAST hits (R1+R2): \$BLAST_HITS" >> ${sample}_${db_name}.blast_summary.txt
    
    if [ "\$BLAST_HITS" -gt 0 ]; then
        echo "" >> ${sample}_${db_name}.blast_summary.txt
        echo "TOP CONTAMINATING SEQUENCES/STRAINS:" >> ${sample}_${db_name}.blast_summary.txt
        
        if [[ "${db_name}" == *"mus_strain"* ]]; then
            cut -f2 ${sample}_${db_name}.blast.tsv | sed 's/gnl|//; s/|.*//' | sort | uniq -c | sort -rn | head -20 >> ${sample}_${db_name}.blast_summary.txt
            cut -f2 ${sample}_${db_name}.blast.tsv | sed 's/gnl|//; s/|.*//' | sort | uniq -c | sort -rn | head -50 > ${sample}_${db_name}.top_contaminants.txt
        else
            cut -f13 ${sample}_${db_name}.blast.tsv | sort | uniq -c | sort -rn | head -20 >> ${sample}_${db_name}.blast_summary.txt
            cut -f13 ${sample}_${db_name}.blast.tsv | sort | uniq -c | sort -rn | head -50 > ${sample}_${db_name}.top_contaminants.txt
        fi
    else
        echo "No contamination detected" > ${sample}_${db_name}.top_contaminants.txt
    fi
    
    echo "BLAST analysis complete for ${sample} vs ${db_name}" >&2
    rm -f ${sample}_R1.fasta ${sample}_R2.fasta ${sample}_R1.sampled.fastq.gz ${sample}_R2.sampled.fastq.gz
    """
}

process BLAST_PLOT_CHARTS {
    tag "${sample}_${db_name}"
    publishDir "${params.outdir}/Output/contamination_check/${sample}/${reads_type}/${db_name}/plots", mode: 'copy'

    input:
    tuple val(sample), val(db_name), path(blast_tsv), path(blast_summary), val(reads_type)

    output:
    path "*.png", optional: true, emit: plots
    path "*.html", optional: true, emit: html

    when:
    params.run_contamination_check

    script:
    """
    set -euo pipefail
    
    # Check if BLAST results have hits
    if [ ! -s ${blast_tsv} ]; then
        echo "No BLAST hits found, skipping plots"
        exit 0
    fi
    
    # Count number of hits
    HITS=\$(wc -l < ${blast_tsv})
    
    if [ "\$HITS" -lt 1 ]; then
        echo "No BLAST hits found, skipping plots"
        exit 0
    fi
    
    echo "Generating plots for ${sample} (${db_name}): \$HITS hits"
    
    # Generate matplotlib pie chart (overview)
    python3 ${params.blast_plot} ${blast_tsv} ${sample}_${db_name}
    
    # Generate interactive plotly chart
    python3 ${params.blast_plot_interac} ${blast_tsv} ${sample}_${db_name}
    
    # Auto-detect top genus and create breakdown
    TOP_GENUS=\$(python3 -c "
import sys
from collections import Counter
import re

def extract_genus(desc):
    clean = re.sub(r'^PREDICTED:\\s+', '', desc)
    clean = re.sub(r'^\\S+:\\s+', '', clean)
    if clean.lower().startswith('mouse'):
        return 'Mus'
    match = re.search(r'^([A-Z][a-z]+)\\s+[a-z]+', clean)
    return match.group(1) if match else clean.split()[0]

genus_counts = Counter()
with open('${blast_tsv}') as f:
    for line in f:
        fields = line.strip().split('\\t')
        if len(fields) >= 13:
            genus_counts[extract_genus(fields[12])] += 1

if genus_counts:
    print(genus_counts.most_common(1)[0][0])
" || echo "")
    
    # Generate genus breakdown if top genus found and has enough hits
    if [ -n "\$TOP_GENUS" ] && [ "\$HITS" -ge 10 ]; then
        echo "Generating breakdown for genus: \$TOP_GENUS"
        python3 ${params.blast_plot} ${blast_tsv} ${sample}_${db_name} --genus "\$TOP_GENUS" || true
    fi
    
    echo "Plot generation complete"
    ls -lh *.png *.html 2>/dev/null || echo "No plot files found"
    """
}

process CONTAMINATION_FINAL_REPORT_MULTI {
    tag "${sample}"
    publishDir "${params.outdir}/Output/contamination_check/${sample}/final_report", mode: 'copy'

    input:
    tuple val(sample), path(mapped_stats), path(unmapped_stats), path(blast_summaries), path(top_hits_files)

    output:
    path "${sample}.contamination_final_summary.txt", emit: summary

    when:
    params.run_contamination_check

    script:
    """
    set -euo pipefail
    
    cat > ${sample}.contamination_final_summary.txt << 'EOF'
COMPREHENSIVE CONTAMINATION REPORT: ${sample}
========================================
EOF
    
    echo "" >> ${sample}.contamination_final_summary.txt
    echo "MAPPED READS ANALYSIS:" >> ${sample}.contamination_final_summary.txt
    cat ${mapped_stats} >> ${sample}.contamination_final_summary.txt
    
    echo "" >> ${sample}.contamination_final_summary.txt
    echo "UNMAPPED READS ANALYSIS:" >> ${sample}.contamination_final_summary.txt
    cat ${unmapped_stats} >> ${sample}.contamination_final_summary.txt
    
    echo "" >> ${sample}.contamination_final_summary.txt
    echo "BLAST CONTAMINATION SCREENING (ALL DATABASES):" >> ${sample}.contamination_final_summary.txt
    echo "===========================================" >> ${sample}.contamination_final_summary.txt
    
    for summary in ${blast_summaries}; do
        echo "" >> ${sample}.contamination_final_summary.txt
        cat \$summary >> ${sample}.contamination_final_summary.txt
        echo "" >> ${sample}.contamination_final_summary.txt
    done
    """
}

process MULTIQC {
    publishDir "${params.outdir}/Output/multiqc", mode: 'copy'

    input:
    path(reports, stageAs: "?/*")  // Stage with subdirectories

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data, optional: true
    path "multiqc_plots", optional: true, emit: plots

    when:
    params.run_multiqc

    script:
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    """
    set -euo pipefail
    
    # Count input files
    FILE_COUNT=\$(find . -type f 2>/dev/null | wc -l)
    echo "MultiQC found \$FILE_COUNT input files"
    
    if [ "\$FILE_COUNT" -eq 0 ]; then
        echo "ERROR: No input files for MultiQC"
        exit 1
    fi
    
    # Run MultiQC
    ${MULTIQC_BIN} . \\
        --force \\
        --title "Conterminator QC Report" \\
        --comment "RNA-seq Quality Control & Contamination Detection Pipeline" \\
        --filename multiqc_report.html \\
        --dirs \\
        --dirs-depth 2 \\
        ${config} \\
        --verbose
    
    echo "MultiQC report generated successfully"
    echo "Found \$FILE_COUNT QC output files"
    """
}

// ====================
// Main Workflow
// ====================
workflow {
    // ===================================
    // PIPELINE METADATA
    // ===================================

    // Write pipeline information to output directory
    WRITE_PIPELINE_INFO(
        file("${projectDir}/main.nf"),
        file("${projectDir}/nextflow.config"),
        file("${projectDir}/conterminator.def")
    )

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
        // MODE 1: Sample sheet with individual file paths per sample
        println "Using sample sheet with per-sample file paths"

        // Create separate channels for different input types
        def fastqList = []
        def bamList = []
        def unmappedFastqList = []

        sampleToStrain.each { sample, info ->
            def type = info.type
            def read1File = file(info.read1)
            def read2File = info.read2 ? file(info.read2) : null

            // Validate files exist
            if (!read1File.exists()) {
                println "ERROR: read1 file not found for sample ${sample}: ${info.read1}"
                System.exit(1)
            }

            if (type == 'fastq') {
                if (!read2File || !read2File.exists()) {
                    println "ERROR: read2 file not found for sample ${sample}: ${info.read2}"
                    System.exit(1)
                }
                fastqList.add([sample, read1File, read2File])
            } else if (type == 'bam') {
                bamList.add([sample, read1File])
            } else if (type == 'unmapped_fastq') {
                if (!read2File || !read2File.exists()) {
                    println "ERROR: read2 file not found for sample ${sample}: ${info.read2}"
                    System.exit(1)
                }
                unmappedFastqList.add([sample, read1File, read2File])
            }
        }

        // Create channels for each input type
        if (!fastqList.isEmpty()) {
            Channel.fromList(fastqList).set { raw_fastq_ch }
            println "Found ${fastqList.size()} FastQ sample(s)"
        } else {
            Channel.empty().set { raw_fastq_ch }
        }

        if (!bamList.isEmpty()) {
            Channel.fromList(bamList).set { raw_bam_ch }
            println "Found ${bamList.size()} BAM sample(s)"
        } else {
            Channel.empty().set { raw_bam_ch }
        }

        if (!unmappedFastqList.isEmpty()) {
            Channel.fromList(unmappedFastqList).set { raw_unmapped_fastq_ch }
            println "Found ${unmappedFastqList.size()} unmapped FastQ sample(s)"
        } else {
            Channel.empty().set { raw_unmapped_fastq_ch }
        }

        // Check that at least one input type was provided
        if (fastqList.isEmpty() && bamList.isEmpty() && unmappedFastqList.isEmpty()) {
            println "ERROR: No valid input files found in sample sheet"
            System.exit(1)
        }

    } else {
        println """
        ========================================
        ERROR: No input provided
        ========================================
        Please provide input via:
          --sample_sheet <file.tsv>

        Sample sheet format:
          sample\tstrain\ttype\tread1\tread2

        Where type is: fastq, bam, or unmapped_fastq
        ========================================
        """
        System.exit(1)
    }
    
    // ===================================
    // SECTION 2: FastQ PROCESSING (if FastQ input)
    // ===================================

    if (params.sample_sheet) {
        // Copy sample sheet to output directory
        COPY_SAMPLE_SHEET(file(params.sample_sheet))
    }

    // Initialize channels that may be defined conditionally
    star_unmapped_fastq = Channel.empty()
    star_aligned_bams = Channel.empty()
    bams_for_qc = Channel.empty()
    all_bams_ch = Channel.empty()
    input_bams = Channel.empty()

    // Process FASTQs (will be skipped if channel is empty)
    // Subset for FastQ QC
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
        // SECTION 3A: STAR ALIGNMENT (for FastQ inputs)
        // ===================================

        if (params.run_star_alignment) {
            // Get strains needed for FastQ samples
            def fastqStrains = []
            sampleToStrain.each { sample, info ->
                if (info.type == 'fastq') {
                    if (!fastqStrains.contains(info.strain)) {
                        fastqStrains.add(info.strain)
                    }
                }
            }

            if (!fastqStrains.isEmpty()) {
                // Sample sheet
                strain_ch = Channel.from(fastqStrains)
                PREPARE_STRAIN_REFERENCE(strain_ch)

                // Prepare per-strain references for QC tools
                if (need_ref) {
                    PREP_STRAIN_REFERENCE_FOR_QC(PREPARE_STRAIN_REFERENCE.out.references)
                    strain_qc_refs = PREP_STRAIN_REFERENCE_FOR_QC.out.qc_references
                }

                strain_index_status = PREPARE_STRAIN_REFERENCE.out.references
                    .map { strain, fasta, gtf ->
                        // Check both <star_index_dir>/<strain>/SA and <star_index_dir>/<strain>/star_index/SA
                        def indexDir = file("${params.star_index_dir}/${strain}")
                        def indexSubDir = file("${params.star_index_dir}/${strain}/star_index")
                        def saFile = file("${params.star_index_dir}/${strain}/SA")
                        def saFileSub = file("${params.star_index_dir}/${strain}/star_index/SA")

                        def indexExists = false
                        if (indexSubDir.exists() && saFileSub.exists()) {
                            indexDir = indexSubDir
                            indexExists = true
                            println "OK: Found existing STAR index for ${strain} in star_index/"
                        } else if (indexDir.exists() && saFile.exists()) {
                            indexExists = true
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

                // Save unmapped FastQ for later use
                star_unmapped_fastq = STAR_ALIGN.out.unmapped_fastq
                star_aligned_bams = STAR_ALIGN.out.aligned_bam
            }
        }

    // ===================================
    // SECTION 3B: PROCESS BAM INPUTS
    // ===================================

    // For BAM inputs, we need to prepare them for QC processing
    // BAMs need to be paired with their strain info
    // Get strains needed for BAM samples
    def bamStrains = []
    sampleToStrain.each { sample, info ->
        if (info.type == 'bam') {
            if (!bamStrains.contains(info.strain)) {
                bamStrains.add(info.strain)
            }
        }
    }

    if (!bamStrains.isEmpty()) {
        println "Processing BAM inputs (skipping FastQ QC and alignment steps)"

        if (!bamStrains.isEmpty() && need_ref) {
            // Prepare references for BAM samples (if not already prepared)
            if (!binding.hasVariable('strain_qc_refs')) {
                bam_strain_ch = Channel.from(bamStrains)
                PREPARE_STRAIN_REFERENCE(bam_strain_ch)
                PREP_STRAIN_REFERENCE_FOR_QC(PREPARE_STRAIN_REFERENCE.out.references)
                strain_qc_refs = PREP_STRAIN_REFERENCE_FOR_QC.out.qc_references
            }
        }

        // Add strain info and create index for BAM inputs
        bam_inputs_with_strain = raw_bam_ch
            .map { sample, bam ->
                def strain = sampleToStrain[sample]?.strain
                if (!strain) {
                    println "WARNING: Sample ${sample} not found in sample sheet, skipping..."
                    return null
                }
                return tuple(sample, strain, bam)
            }
            .filter { it != null }

        // Index the BAM files
        INDEX_INPUT_BAM(bam_inputs_with_strain)
        input_bams = INDEX_INPUT_BAM.out.bam_with_index
    }

    // ===================================
    // SECTION 4: BAM QC TOOLS
    // ===================================

    // Combine all BAM sources (Nextflow handles empty channels automatically)
    all_bams_ch = all_bams_ch.mix(star_aligned_bams).mix(input_bams)

    // Process all BAMs for QC
    if (params.subset_bam_for_qc) {
        SUBSET_BAM_FOR_QC(all_bams_ch)
        bams_for_qc = SUBSET_BAM_FOR_QC.out.subsampled_bam
        qc_bams = bams_for_qc.map { sample, strain, bam, bai ->
            tuple(sample, strain, bam, bai, "subset")
        }
    } else {
        bams_for_qc = all_bams_ch
        qc_bams = all_bams_ch.map { sample, strain, bam, bai ->
            tuple(sample, strain, bam, bai, "full")
        }
    }

    // Join BAMs with their strain-specific references for QC tools
    if (need_ref) {
        qc_bams_with_refs = qc_bams
                .map { sample, strain, bam, bai, type -> tuple(strain, sample, bam, bai, type) }
                .combine(strain_qc_refs, by: 0)
                .map { strain, sample, bam, bai, type, fa, fai, dict, bed, eff_size, twobit, gtf ->
                    tuple(sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit, gtf)
                }

        // Run QC tools
        if (params.run_deeptools) {
            qc_bams_for_deeptools = qc_bams_with_refs.map { sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit, gtf ->
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
            qc_bams_for_picard = qc_bams_with_refs.map { sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit, gtf ->
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
            qc_bams_for_bedtools = qc_bams_with_refs.map { sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit, gtf ->
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
            qc_bams_for_mapinsights = qc_bams_with_refs.map { sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit, gtf ->
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
            qc_bams_for_qualimap = qc_bams_with_refs.map { sample, strain, bam, bai, type, fa, fai, dict, bed, eff_size, twobit, gtf ->
                tuple("${sample}_${strain}", strain, bam, bai, type, fa, gtf)
            }

            if (params.qualimap_mode == "bamqc" || params.qualimap_mode == "both") {
                // Group by reference (strain) and run once per strain
                qc_bams_for_qualimap
                    .map { id, strain, bam, bai, type, fa, gtf -> tuple(fa.name, id, strain, bam, bai, type, fa, gtf) }
                    .groupTuple(by: 0)
                    .flatMap { key, ids, strains, bams, bais, types, fas, gtfs ->
                        [ids, strains, bams, bais, types].transpose().collect { id, strain, bam, bai, type ->
                            tuple(id, strain, bam, bai, type, fas[0], gtfs[0])  // Include GTF here
                        }
                    }
                    .set { qualimap_grouped }
                QUALIMAP_BAMQC(
                    qualimap_grouped.map { id, strain, bam, bai, type, fa, gtf -> tuple(id, strain, bam, bai, type, gtf) },  // Pass GTF
                    qualimap_grouped.first().map { id, strain, bam, bai, type, fa, gtf -> fa }
                )
                ran_qualimap_bamqc = true
            }

            if (params.qualimap_mode == "rnaseq" || params.qualimap_mode == "both") {
                QUALIMAP_RNASEQ(qc_bams_for_qualimap.map { id, strain, bam, bai, type, fa, gtf -> tuple(id, strain, bam, bai, type, gtf) })
                ran_qualimap_rnaseq = true
            }
        }
    }

    // ===================================
    // SECTION 5: CONTAMINATION CHECK
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
    def singularityUsed = workflow.profile.contains('singularity')
    def singularityInfo = singularityUsed ? "\n    Singularity: ${params.singularity_path ?: 'enabled'}" : ''

    println """
    ========================================
    QC Pipeline Completed
    ========================================
    Status:   ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    Output:   ${params.outdir}
    Input:    ${params.sample_sheet}
    ${params.sample_sheet ? "Sample Sheet: ${params.sample_sheet}" : ''}
    ${params.strain ? "Strains:  ${params.strain}" : ''}${singularityInfo}

    SUBSAMPLING SUMMARY:
    ${params.subset_for_fastq_qc ? "  FastQ QC: ${params.subset_fastq_qc_reads} reads" : "  X FastQ QC: Full FASTQs"}
    ${params.subset_for_star ? "  STAR: ${params.subset_star_reads} reads" : "  X STAR: Full FASTQs"}
    ${params.subset_bam_for_qc ? "  BAM QC: ${params.bam_qc_subset_mapped} reads" : "  X BAM QC: Full BAMs"}
    ${params.subset_unmapped_for_blast ? "  BLAST: ${params.unmapped_subset_reads} reads" : "  X BLAST: Full unmapped"}

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
