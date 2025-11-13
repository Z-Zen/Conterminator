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
params.subset_for_fastq_qc      = true      // Subset for FastQ Screen
params.subset_fastq_qc_reads    = 100000    // Size for FastQ QC (100k reads)
params.subset_for_star          = false     // Use full FastQ for STAR alignment
params.subset_star_reads        = 100000    // Size for STAR if subsampling (100k reads)
params.subset_seed              = 100       // Random seed for reproducibility

// STAR alignment
params.run_star_alignment       = true
params.star_threads             = 20
params.strain                   = null
params.strains_base_dir         = "/home/${params.user_home_dir}/rcp_storage/common/Users/vonalven/HDP_pseudogenomes_construction/Data/HPC_results/HDP_pseudogenomes"
params.standard_references_dir  = "/mnt/sas/Data/References/Mus"  // For standard reference genomes like GRCm39, GRCm38
params.star_index_dir           = "/mnt/sas/Data/References/Mus"
params.star_index_threads       = 20
params.star_index_mem           = "32 GB"

// BAM SUBSETTING FOR QC TOOLS
params.subset_bam_for_qc                 = true      // Subset BAMs for QC tools
params.bam_qc_subset_mapped              = 200000   // Mapped reads for QC (200k reads)

// READS SUBSETTING
params.subset_unmapped_for_blast         = true   // Subset for BLAST
params.subset_mapped_for_blast           = true   // Subset for BLAST
params.subset_unmapped_for_decontaminer  = true   // Subset for DecontaMiner
params.subset_mapped_for_decontaminer    = true   // Subset for DecontaMiner
params.unmapped_subset_reads             = 100000 // Unmapped reads to keep (100k reads)
params.mapped_subset_reads               = 100000 // Unmapped reads to keep (100k reads)

// Tool toggles
params.run_deeptools           = true
params.run_picard_gc           = true
params.run_fastq_screen        = true
params.run_bedtools_gc         = true
params.run_decontaminer        = true
params.run_mapinsights         = true
params.run_contamination_check = true
params.run_fastqc              = true

// Deeptools settings
params.blacklist_bed      = null

// BEDTools windowing
params.windowsize         = 500
params.window_step        = 250

// DecontaMiner settings
// Use container paths when running with Singularity, host paths otherwise
params.decontaminer_dir             = "/mnt/sas/Tools/decontaminer-RNAseq/decontaMiner_1.4"
params.decontaminer_config          = "/mnt/sas/Tools/decontaminer-RNAseq/decontaMiner_1.4/config_files/configure.txt"
params.decontaminer_pairing         = "P"
params.decontaminer_organisms       = "bfv"
params.decontaminer_format          = "bam"
params.decontaminer_quality_filter  = "yes"
params.decontaminer_ribo_filter     = "yes"
params.decontaminer_gap             = 5
params.decontaminer_mismatch        = 10
params.decontaminer_match_len       = 50
params.decontaminer_match_threshold = 5
params.decontaminer_generate_plots  = true

// Contamination check
params.contamination_gc_min    = 0.60
params.contamination_gc_max    = 1.0
params.contamination_mapq      = 10
params.contamination_blast_dbs = "/mnt/sas/Tools/blast_databases/"
params.contamination_evalue    = "1e-10"

// Mapinsights
params.mapinsights_opts   = ""

// Tool paths
params.bedtools_bin       = "/mnt/sas/Tools/bedtools2/bin/bedtools"
params.mapinsights_bin    = "/mnt/sas/Tools/mapinsights/mapinsights"
params.picard_jar         = "/mnt/sas/Tools/Picard/picard.jar"
params.samtools_bin       = "/mnt/sas/Tools/samtools-1.22.1/samtools"
params.blastn_bin         = "/mnt/sas/Tools/ncbi-blast-2.17.0+/bin/blastn"
params.deeptools_gcbias   = "/mnt/sas/Tools/deepTools/deeptools/computeGCBias.py"
params.facount_bin        = "/mnt/sas/Tools/faCount/faCount"
params.fa2bit_bin         = "/mnt/sas/Tools/faToTwoBit/faToTwoBit"
params.fastq_screen       = "/mnt/sas/Tools/FastQ-Screen-0.16.0/fastq_screen"
params.qualimap_bin       = "/mnt/sas/Tools/qualimap_v2.3/qualimap"
params.seqtk_bin          = "/mnt/sas/Tools/seqtk/seqtk"
params.star_bin           = "/mnt/sas/Tools/STAR_2.7.11b/Linux_x86_64_static/STAR"
params.fastqc_bin         = "/mnt/sas/Tools/FastQC/fastqc"
params.reformat_bin       = "/mnt/sas/Tools/bbmap/reformat.sh"

// Qualimap settings
params.run_qualimap       = true
params.qualimap_mode      = "both"
params.qualimap_genome    = "mm10"
params.qualimap_java_mem  = "10G"
params.qualimap_threads   = 8
params.qualimap_skip_dup  = false
params.qualimap_protocol  = "strand-specific-reverse"

// Resources
params.max_cpus               = 16
params.max_mem                = "32 GB"
params.max_time               = "24h"
params.max_parallel_samples   = 4

// FastQ Screen conf
params.fastq_screen_conf   = "/mnt/sas/Tools/FastQ-Screen-0.16.0/FastQ_Screen_Genomes/fastq_screen.conf"
params.fastq_screen_threads = 8

// Plotting and reporting
params.gc_plot_script         =  "${projectDir}/bin/R/plot_gc_content.R"
params.run_blast_plots        = true
params.blast_plot             = "${projectDir}/bin/py/plot_blast_pie.py"
params.blast_plot_interac     = "${projectDir}/bin/py/interactive_plot_blast_pie.py"
params.run_multiqc            = true
params.multiqc_bin            = "multiqc"
params.multiqc_config         = null

params.singularity_path = '/absolute/path/to/conterminator.sif'

def content = new File('name.txt').text

println content

// def welcomeMessage(){
//     log.info content
// }

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
      --strains_base_dir <path>        Directory with strain-specific pseudogenomes
                                       (default: /home/${params.user_home_dir}/.../HDP_pseudogenomes)
      --standard_references_dir <path> Directory with standard reference genomes (e.g., GRCm39, GRCm38)
                                       (default: /mnt/sas/Data/References/Mus)
      --star_index_dir <path>          Directory for STAR indices (default: /mnt/sas/Data/References/Mus)
    
    STAR ALIGNMENT:
      --run_star_alignment        Enable STAR alignment (default: true)
      --strain <strain(s)>        Comma-separated strain list (alternative to sample sheet)
                                  Example: "C57BL_6J,DBA_2J"
      --star_threads <int>        STAR alignment threads (default: 20)
      --star_index_threads <int>  STAR index building threads (default: 20)
      --star_index_mem <string>   STAR index memory (default: "32 GB")
    
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
       - FastQ QC: 100k reads
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
    - path: Directory containing paired-end FastQ files
    
    Notes:
    - Header row is mandatory
    - Each sample is aligned to its specified strain
    - FastQ pattern will be auto-detected (or use --fastq_pattern)
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
    
    Project: Conterminator v1.1
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


if (params.outdir == null){
    println "Please specify the parameter '--outdir [output folder]'"
    println "Run '~/nextflow run main.nf --help' for more information"
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
    
    def lines = sheetFile.readLines()
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
        if (line.trim() && !line.startsWith('#')) {
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
// When using Singularity, use container paths; otherwise use host paths
def usingSingularity = workflow.profile.contains('singularity')

// Core alignment tools
def STAR_BIN = usingSingularity ? "/opt/tools/STAR/bin/STAR" : params.star_bin
def SAMTOOLS_BIN = usingSingularity ? "/opt/tools/bin/samtools" : params.samtools_bin
def SEQTK_BIN = usingSingularity ? "/opt/tools/bin/seqtk" : params.seqtk_bin

// QC tools
def FASTQC_BIN = usingSingularity ? "/opt/tools/bin/fastqc" : params.fastqc_bin
def FASTQ_SCREEN_BIN = usingSingularity ? "/opt/tools/bin/fastq_screen" : params.fastq_screen
def QUALIMAP_BIN = usingSingularity ? "/opt/tools/bin/qualimap" : params.qualimap_bin

// Analysis tools
def BEDTOOLS_BIN = usingSingularity ? "/opt/tools/bin/bedtools" : params.bedtools_bin
def PICARD_JAR = usingSingularity ? "/opt/tools/picard/picard.jar" : params.picard_jar
def DEEPTOOLS_GCBIAS = usingSingularity ? "computeGCBias" : params.deeptools_gcbias
def BLASTN_BIN = usingSingularity ? "/opt/tools/blast/bin/blastn" : params.blastn_bin

// UCSC tools
def FACOUNT_BIN = usingSingularity ? "/opt/tools/bin/faCount" : params.facount_bin
def FA2BIT_BIN = usingSingularity ? "/opt/tools/bin/faToTwoBit" : params.fa2bit_bin

// Other tools
def REFORMAT_BIN = usingSingularity ? "/opt/tools/bbmap/reformat.sh" : params.reformat_bin
def MULTIQC_BIN = usingSingularity ? "multiqc" : params.multiqc_bin
def MAPINSIGHTS_BIN = usingSingularity ? "/opt/tools/bin/mapinsights" : params.mapinsights_bin

// DecontaMiner paths
def DECONTAMINER_DIR = usingSingularity ? "/opt/tools/decontaminer" : params.decontaminer_dir
def DECONTAMINER_CONFIG = usingSingularity ? "/opt/tools/decontaminer/config_files/configure.txt" : params.decontaminer_config

if (usingSingularity) {
    println "INFO: Using Singularity - Tool paths configured for container execution"
}

if (params.run_decontaminer) {
    if (!DECONTAMINER_CONFIG) {
        exit 1, "ERROR: DecontaMiner requires --decontaminer_config"
    }
    // Skip file existence check when using Singularity (file is inside container)
    if (!usingSingularity && !file(DECONTAMINER_CONFIG).exists()) {
        exit 1, "ERROR: DecontaMiner config file not found: ${DECONTAMINER_CONFIG}"
    }
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
Unmapped Decon Subset:   ${params.subset_unmapped_for_decontaminer} (${params.unmapped_subset_reads} reads)
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
DecontaMiner: ${params.run_decontaminer}
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
    tuple val(strain), path("${strain}_genome.fa.gz"), path("${strain}_annotation.gtf.gz"), emit: references

    script:
    """
    set -euo pipefail

    # Check both possible locations for the strain
    PSEUDOGENOME_DIR="${params.strains_base_dir}/${strain}"
    STANDARD_REF_DIR="${params.standard_references_dir}/${strain}"

    # Determine which directory to use
    if [ -d "\${PSEUDOGENOME_DIR}" ]; then
        STRAIN_DIR="\${PSEUDOGENOME_DIR}"
        echo "INFO: Using strain-specific pseudogenome directory: \${STRAIN_DIR}" >&2
    elif [ -d "\${STANDARD_REF_DIR}" ]; then
        STRAIN_DIR="\${STANDARD_REF_DIR}"
        echo "INFO: Using standard reference directory: \${STRAIN_DIR}" >&2
    else
        echo "ERROR: Could not find directory for ${strain}" >&2
        echo "Searched in:" >&2
        echo "  - \${PSEUDOGENOME_DIR}" >&2
        echo "  - \${STANDARD_REF_DIR}" >&2
        exit 1
    fi

    # First, try to find strain-specific pseudogenome files
    FASTA=\$(find "\${STRAIN_DIR}" -name "*pseudogenome__strain_${strain}.fa.gz" 2>/dev/null | head -1)

    # If not found, look for any FASTA file
    if [ -z "\${FASTA}" ]; then
        echo "INFO: Strain-specific pseudogenome not found, searching for standard reference genome..." >&2
        FASTA=\$(find "\${STRAIN_DIR}" -type f \\( -name "*.fa.gz" -o -name "*.fasta.gz" -o -name "*.genome.fa.gz" \\) 2>/dev/null | head -1)
    fi

    if [ -z "\${FASTA}" ]; then
        echo "ERROR: Could not find FASTA file for ${strain} in \${STRAIN_DIR}" >&2
        echo "Searched for:" >&2
        echo "  - *pseudogenome__strain_${strain}.fa.gz" >&2
        echo "  - *.fa.gz, *.fasta.gz, *.genome.fa.gz" >&2
        exit 1
    fi

    # First, try to find strain-specific GTF
    GTF=\$(find "\${STRAIN_DIR}" -name "*pseudogenome__strain_${strain}.gtf.gz" 2>/dev/null | head -1)

    # If not found, look for any GTF file
    if [ -z "\${GTF}" ]; then
        echo "INFO: Strain-specific GTF not found, searching for standard annotation..." >&2
        GTF=\$(find "\${STRAIN_DIR}" -type f \\( -name "*.gtf.gz" -o -name "*.gff.gz" -o -name "*.gff3.gz" \\) 2>/dev/null | head -1)
    fi

    if [ -z "\${GTF}" ]; then
        echo "ERROR: Could not find GTF/GFF file for ${strain} in \${STRAIN_DIR}" >&2
        echo "Searched for:" >&2
        echo "  - *pseudogenome__strain_${strain}.gtf.gz" >&2
        echo "  - *.gtf.gz, *.gff.gz, *.gff3.gz" >&2
        exit 1
    fi

    ln -s "\${FASTA}" ${strain}_genome.fa.gz
    ln -s "\${GTF}" ${strain}_annotation.gtf.gz

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
    
    ${STAR_BIN} \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${strain} \\
        --genomeFastaFiles \${FASTA_FILE} \\
        --sjdbGTFfile \${GTF_FILE} \\
        --sjdbOverhang 149
    
    sed 's///g'

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
    publishDir "${params.outdir}/Input/strain_references/${strain}", mode: 'copy'

    input:
    tuple val(strain), path(fasta_gz), path(gtf)

    output:
    tuple val(strain), path("${strain}.fa"), path("${strain}.fa.fai"), path("${strain}.dict"), path("${strain}.bed"), path("${strain}_effective_size.txt"), path("${strain}.2bit"), emit: qc_references

    script:
    """
    set -euo pipefail

    echo "Preparing QC references for strain: ${strain}"

    # Decompress fasta
    zcat ${fasta_gz} > ${strain}.fa

    # Create fasta index
    ${SAMTOOLS_BIN} faidx -@ ${task.cpus} ${strain}.fa

    # Create sequence dictionary
    ${SAMTOOLS_BIN} dict ${strain}.fa > ${strain}.dict

    # Create BED file from fai
    awk 'BEGIN{OFS="\\t"} {print \$1, 0, \$2}' ${strain}.fa.fai > ${strain}.bed

    # Calculate effective genome size
    ${FACOUNT_BIN} ${strain}.fa > ${strain}.facount.txt
    awk 'BEGIN{len=0; n=0} \$1=="total"{len=\$2; n=\$7} END{print len-n}' ${strain}.facount.txt > ${strain}_effective_size.txt
    echo "Effective genome size for ${strain}: \$(cat ${strain}_effective_size.txt)" >&2

    # Create 2bit file if deeptools is enabled
    if [[ "${params.run_deeptools}" == "true" ]]; then
        ${FA2BIT_BIN} ${strain}.fa ${strain}.2bit
    else
        # Create empty placeholder
        touch ${strain}.2bit
    fi

    echo "QC references prepared for ${strain}"
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

process SUBSET_UNMAPPED_FOR_DECONTAMINER {
    tag "${sample}_${strain}"
    publishDir "${params.outdir}/Input/subsampled_unmapped_decontaminer/${sample}/${strain}", mode: 'copy'

    input:
    tuple val(sample), val(strain), path(unmapped_r1), path(unmapped_r2)

    output:
    tuple val(sample), val(strain), path("${sample}_${strain}_unmapped_decon_subset_1.fq"), path("${sample}_${strain}_unmapped_decon_subset_2.fq"), emit: subsampled_fastq
    path "${sample}_${strain}_unmapped_decon_subset_stats.txt", emit: stats

    when:
    params.subset_unmapped_for_decontaminer

    script:
    """
    set -euo pipefail
    
    echo "Subsampling Unmapped for DecontaMiner: ${sample}_${strain}" > ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "=========================================================" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "Target: ${params.unmapped_subset_reads} reads per file" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    
    TOTAL_R1=\$(cat ${unmapped_r1} | wc -l | awk '{print \$1/4}')
    TOTAL_R2=\$(cat ${unmapped_r2} | wc -l | awk '{print \$1/4}')
    
    echo "Total unmapped R1: \$TOTAL_R1" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "Total unmapped R2: \$TOTAL_R2" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    
    # Subsample and add /1 suffix to R1 reads
    if [ "\$TOTAL_R1" -gt "${params.unmapped_subset_reads}" ]; then
        ${SEQTK_BIN} sample -s${params.subset_seed} ${unmapped_r1} ${params.unmapped_subset_reads} | \\
            sed '1~4 s/\$/\\/1/' > ${sample}_${strain}_unmapped_decon_subset_1.fq
        echo "Subsampled R1 to ${params.unmapped_subset_reads} reads" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    else
        sed '1~4 s/\$/\\/1/' ${unmapped_r1} > ${sample}_${strain}_unmapped_decon_subset_1.fq
        echo "Using all \$TOTAL_R1 unmapped R1 reads" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    fi
    
    # Subsample and add /2 suffix to R2 reads
    if [ "\$TOTAL_R2" -gt "${params.unmapped_subset_reads}" ]; then
        ${SEQTK_BIN} sample -s${params.subset_seed} ${unmapped_r2} ${params.unmapped_subset_reads} | \\
            sed '1~4 s/\$/\\/2/' > ${sample}_${strain}_unmapped_decon_subset_2.fq
        echo "Subsampled R2 to ${params.unmapped_subset_reads} reads" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    else
        sed '1~4 s/\$/\\/2/' ${unmapped_r2} > ${sample}_${strain}_unmapped_decon_subset_2.fq
        echo "Using all \$TOTAL_R2 unmapped R2 reads" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    fi
    
    FINAL_R1=\$(cat ${sample}_${strain}_unmapped_decon_subset_1.fq | wc -l | awk '{print \$1/4}')
    FINAL_R2=\$(cat ${sample}_${strain}_unmapped_decon_subset_2.fq | wc -l | awk '{print \$1/4}')
    echo "Final R1 reads: \$FINAL_R1" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "Final R2 reads: \$FINAL_R2" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "Purpose: DecontaMiner analysis" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
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
    
    python3 ${DEEPTOOLS_GCBIAS} \\
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
    export PATH="/opt/R/4.5.1/bin:\$PATH"

    java -jar ${PICARD_JAR} CollectGcBiasMetrics \\
        I=${bam} \\
        O=${sample}.gc_bias_metrics.txt \\
        SUMMARY_OUTPUT=${sample}.gc_bias_summary.txt \\
        CHART_OUTPUT=${sample}.gc_bias.pdf \\
        R=${ref} \\
        VALIDATION_STRINGENCY=LENIENT

    java -jar ${PICARD_JAR} CollectMultipleMetrics \\
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

    export PATH="/opt/R/4.5.1/bin:\$PATH"

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
    tuple val(sample), val(strain), path(bam), path(bai), val(bam_type)
    path ref

    output:
    path "qualimap_bamqc/**", emit: reports
    path "qualimap_bamqc/qualimapReport.html", optional: true, emit: html
    path "qualimap_bamqc/qualimapReport.pdf", optional: true, emit: pdf
    path "qualimap_bamqc/genome_results.txt", optional: true, emit: results

    when:
    params.run_qualimap && (params.qualimap_mode == "bamqc" || params.qualimap_mode == "both")

    script:
    def gtf = "/home/${params.user_home_dir}/rcp_storage/common/Users/vonalven/HDP_pseudogenomes_construction/Data/HPC_results/HDP_pseudogenomes/${strain}/HDP_merge_splitnorm_v1__pseudogenome__strain_${strain}.gtf.gz"
    """
    set -euo pipefail

    if [[ ${gtf} == *.gz ]]; then
        zcat ${gtf} > ${strain}_annotation.gtf
        GTF_FILE="${strain}_annotation.gtf"
    else
        GTF_FILE="${gtf}"
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
    tuple val(sample), val(strain), path(bam), path(bai), val(bam_type)

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
    def gtf = "/home/${params.user_home_dir}/rcp_storage/common/Users/vonalven/HDP_pseudogenomes_construction/Data/HPC_results/HDP_pseudogenomes/${strain}/HDP_merge_splitnorm_v1__pseudogenome__strain_${strain}.gtf.gz"
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
    
    export PATH="/opt/R/4.5.1/bin:\$PATH"

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
        --conf ${params.fastq_screen_conf} \\
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
        --threads ${task.cpus} \\
        ${fastq1} ${fastq2}
    """
}

process DECONTAMINER_STEP1_STAR_MAPPED {
    tag "${sample} ${reads_type}"
    publishDir "${params.outdir}/Temporary/decontaminer/${sample}/step1_decontaminer_mapped", mode: 'copy'

    input:
    tuple val(sample), path(mapped_bam), path(mapped_bai), val(reads_type)

    output:
    tuple val(sample), val(reads_type), path("decontaminer_output"), emit: output_dir
    path "decontaminer_output/RESULTS/**", optional: true

    when:
    params.run_decontaminer

    script:
    def organisms = ""
    if (params.decontaminer_organisms.contains('b')) organisms += " -b"
    if (params.decontaminer_organisms.contains('f')) organisms += " -f"
    if (params.decontaminer_organisms.contains('v')) organisms += " -v"
    
    def quality_filter = params.decontaminer_quality_filter ? "-Q ${params.decontaminer_quality_filter}" : ""
    def ribo_filter = params.decontaminer_ribo_filter ? "-R ${params.decontaminer_ribo_filter}" : ""
    """
    set -euo pipefail

    # Fix permissions on DecontaMiner scripts if needed
    if [ -f "${DECONTAMINER_DIR}/shell_scripts/decontaMiner.sh" ]; then
        chmod +x ${DECONTAMINER_DIR}/shell_scripts/decontaMiner.sh 2>/dev/null || true
        chmod +x ${DECONTAMINER_DIR}/*.pl ${DECONTAMINER_DIR}/*.sh 2>/dev/null || true
    fi
    
    SAMPLE="${sample}"
    mkdir -p input_fastq

    ${SAMTOOLS_BIN} fastq -@ ${task.cpus} -1 input_fastq/\${SAMPLE}_mapped_1.fq -2 input_fastq/\${SAMPLE}_mapped_2.fq -0 /dev/null -s /dev/null -N ${mapped_bam}
    
    echo "Running DecontaMiner on MAPPED reads...."
    bash ${DECONTAMINER_DIR}/shell_scripts/decontaMiner.sh \\
        -i \$(pwd)/input_fastq \\
        -o \$(pwd)/decontaminer_output \\
        -c ${DECONTAMINER_CONFIG} \\
        -F fastq \\
        -s ${params.decontaminer_pairing} \\
        ${quality_filter} \\
        ${ribo_filter} \\
        ${organisms} || echo "DecontaMiner step 1 (mapped) completed with warnings"
    
    mkdir -p decontaminer_output
    """
}

process DECONTAMINER_STEP1_STAR_UNMAPPED {
    tag "${sample} ${reads_type}"
    publishDir "${params.outdir}/Temporary/decontaminer/${sample}/step1_decontaminer_unmapped", mode: 'copy'

    input:
    tuple val(sample), path(unmapped_r1), path(unmapped_r2), val(reads_type)

    output:
    tuple val(sample), val(reads_type), path("decontaminer_output"), emit: output_dir
    path "decontaminer_output/RESULTS/**", optional: true

    when:
    params.run_decontaminer

    script:
    def organisms = ""
    if (params.decontaminer_organisms.contains('b')) organisms += " -b"
    if (params.decontaminer_organisms.contains('f')) organisms += " -f"
    if (params.decontaminer_organisms.contains('v')) organisms += " -v"
    
    def quality_filter = params.decontaminer_quality_filter ? "-Q ${params.decontaminer_quality_filter}" : ""
    def ribo_filter = params.decontaminer_ribo_filter ? "-R ${params.decontaminer_ribo_filter}" : ""
    
    """
    set -euo pipefail

    # Fix permissions on DecontaMiner scripts if needed
    if [ -f "${DECONTAMINER_DIR}/shell_scripts/decontaMiner.sh" ]; then
        chmod +x ${DECONTAMINER_DIR}/shell_scripts/decontaMiner.sh 2>/dev/null || true
        chmod +x ${DECONTAMINER_DIR}/*.pl ${DECONTAMINER_DIR}/*.sh 2>/dev/null || true
    fi
    
    SAMPLE="${sample}"
    mkdir -p input_fastq
    
    # Add /1 and /2 suffixes to read names (only to header lines)
    awk 'NR%4==1 {sub(/^@/, ""); print "@"\$0"/1"} NR%4!=1' ${unmapped_r1} > input_fastq/\${SAMPLE}_unmapped_1.fq
    awk 'NR%4==1 {sub(/^@/, ""); print "@"\$0"/2"} NR%4!=1' ${unmapped_r2} > input_fastq/\${SAMPLE}_unmapped_2.fq
    
    echo "Running DecontaMiner on UNMAPPED reads...."
    bash ${DECONTAMINER_DIR}/shell_scripts/decontaMiner.sh \\
        -i \$(pwd)/input_fastq \\
        -o \$(pwd)/decontaminer_output \\
        -c ${DECONTAMINER_CONFIG} \\
        -F fastq \\
        -s ${params.decontaminer_pairing} \\
        ${quality_filter} \\
        ${ribo_filter} \\
        ${organisms} || echo "DecontaMiner step 1 (unmapped) completed with warnings"
    """
}

process DECONTAMINER_STEP2 {
    tag "${sample} ${reads_type}"
    publishDir "${params.outdir}/Temporary/decontaminer/${sample}/step2_filtering", mode: 'copy'

    input:
    tuple val(sample), val(reads_type), path(output_dir)

    output:
    tuple val(sample), val(reads_type), path("${output_dir}"), emit: filtered_dir

    when:
    params.run_decontaminer

    script:
    """
    set -euo pipefail

    # Fix permissions on DecontaMiner scripts if needed
    if [ -f "${DECONTAMINER_DIR}/shell_scripts/filterBlastInfo.sh" ]; then
        chmod +x ${DECONTAMINER_DIR}/shell_scripts/*.sh 2>/dev/null || true
        chmod +x ${DECONTAMINER_DIR}/*.pl ${DECONTAMINER_DIR}/*.sh 2>/dev/null || true
    fi
    
    if [ -d "${output_dir}/RESULTS/BACTERIA" ]; then
        bash ${DECONTAMINER_DIR}/shell_scripts/filterBlastInfo.sh \\
            -i \$(pwd)/${output_dir}/RESULTS/BACTERIA/ \\
            -s ${params.decontaminer_pairing} \\
            -g ${params.decontaminer_gap} \\
            -m ${params.decontaminer_mismatch} \\
            -l ${params.decontaminer_match_len} \\
            -V O || echo "Bacteria/Fungi filtering completed with warnings"
    fi

    if [ -d "${output_dir}/RESULTS/FUNGI" ]; then
        bash ${DECONTAMINER_DIR}/shell_scripts/filterBlastInfo.sh \\
            -i \$(pwd)/${output_dir}/RESULTS/FUNGI/ \\
            -s ${params.decontaminer_pairing} \\
            -g ${params.decontaminer_gap} \\
            -m ${params.decontaminer_mismatch} \\
            -l ${params.decontaminer_match_len} \\
            -V O || echo "Bacteria/Fungi filtering completed with warnings"
    fi
    
    if [ -d "${output_dir}/RESULTS/VIRUSES" ]; then
        bash ${DECONTAMINER_DIR}/shell_scripts/filterBlastInfo.sh \\
            -i \$(pwd)/${output_dir}/RESULTS/VIRUSES \\
            -s ${params.decontaminer_pairing} \\
            -g ${params.decontaminer_gap} \\
            -m ${params.decontaminer_mismatch} \\
            -l ${params.decontaminer_match_len} \\
            -V V || echo "Virus filtering completed with warnings"
    fi
    
    echo "Filtering step completed for ${sample}"
    """
}

process DECONTAMINER_STEP3 {
    tag "${sample} ${reads_type}"
    publishDir "${params.outdir}/Output/decontaminer/${sample}/", mode: 'copy'

    input:
    tuple val(sample), val(reads_type), path(filtered_dir)

    output:
    path "bacteria_reports/**", optional: true, emit: bacteria_reports
    path "fungi_reports/**", optional: true, emit: fungi_reports
    path "virus_reports/**", optional: true, emit: virus_reports
    path "HTML_REPORTS/**", optional: true, emit: html_reports

    when:
    params.run_decontaminer

    script:
    def plots_flag = params.decontaminer_generate_plots ? "-P y" : ""
    """
    set -euo pipefail

    # Fix permissions on DecontaMiner scripts if needed
    if [ -f "${DECONTAMINER_DIR}/shell_scripts/collectInfo.sh" ]; then
        chmod +x ${DECONTAMINER_DIR}/shell_scripts/*.sh 2>/dev/null || true
        chmod +x ${DECONTAMINER_DIR}/*.pl ${DECONTAMINER_DIR}/*.sh 2>/dev/null || true
    fi

    export PATH="/opt/R/4.5.1/bin:\$PATH"

    if [ -d "${filtered_dir}/RESULTS/BACTERIA/COLLECTED_INFO" ]; then
        bash ${DECONTAMINER_DIR}/shell_scripts/collectInfo.sh \\
            -i \$(pwd)/${filtered_dir}/RESULTS/BACTERIA/COLLECTED_INFO \\
            -t ${params.decontaminer_match_threshold} \\
            -V O \\
            ${plots_flag} || echo "Bacteria collection completed with warnings"

        mkdir -p bacteria_reports
        # Copy COLLECTED_INFO outputs
        cp -r ${filtered_dir}/RESULTS/BACTERIA/COLLECTED_INFO/* bacteria_reports/ 2>/dev/null || true
    fi

    if [ -d "${filtered_dir}/RESULTS/FUNGI/COLLECTED_INFO" ]; then
        bash ${DECONTAMINER_DIR}/shell_scripts/collectInfo.sh \\
            -i \$(pwd)/${filtered_dir}/RESULTS/FUNGI/COLLECTED_INFO \\
            -t ${params.decontaminer_match_threshold} \\
            -V O \\
            ${plots_flag} || echo "Fungi collection completed with warnings"

        mkdir -p fungi_reports
        cp -r ${filtered_dir}/RESULTS/FUNGI/COLLECTED_INFO/* fungi_reports/ 2>/dev/null || true
    fi

    if [ -d "${filtered_dir}/RESULTS/VIRUSES/COLLECTED_INFO" ]; then
        bash ${DECONTAMINER_DIR}/shell_scripts/collectInfo.sh \\
            -i \$(pwd)/${filtered_dir}/RESULTS/VIRUSES/COLLECTED_INFO \\
            -t ${params.decontaminer_match_threshold} \\
            -V V \\
            ${plots_flag} || echo "Virus collection completed with warnings"

        mkdir -p virus_reports
        cp -r ${filtered_dir}/RESULTS/VIRUSES/COLLECTED_INFO/* virus_reports/ 2>/dev/null || true
    fi

    if [ -d "${filtered_dir}/RESULTS/HTML_REPORTS" ]; then
        echo "Copying shared HTML_REPORTS directory..."
        cp -r ${filtered_dir}/RESULTS/HTML_REPORTS . || echo "Warning: Could not copy HTML_REPORTS"
    else
        echo "No HTML_REPORTS directory found - this is normal if no valid contamination was detected"
    fi
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
        --threads ${task.cpus} \\
        ${fastq_r1}
    
    # Run FastQC on R2
    ${FASTQC_BIN} \\
        --outdir . \\
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
    set -euo pipefail
    
    export BLASTDB=/mnt/sas/Tools/blast_databases
    
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
    
    echo "1"

    # Combine results
    cat ${sample}_${db_name}_R1.blast.tsv ${sample}_${db_name}_R2.blast.tsv > ${sample}_${db_name}.blast.tsv

    echo "2"
    
    BLAST_HITS=\$(wc -l < ${sample}_${db_name}.blast.tsv)
    echo "Total BLAST hits (R1+R2): \$BLAST_HITS" >> ${sample}_${db_name}.blast_summary.txt

    echo "3"
    
    if [ "\$BLAST_HITS" -gt 0 ]; then
        echo "" >> ${sample}_${db_name}.blast_summary.txt
        echo "TOP CONTAMINATING SEQUENCES/STRAINS:" >> ${sample}_${db_name}.blast_summary.txt

        echo "33"
        
        if [[ "${db_name}" == *"mus_strain"* ]]; then
            echo "333"
            cut -f2 ${sample}_${db_name}.blast.tsv | sed 's/gnl|//; s/|.*//' | sort | uniq -c | sort -rn | head -20 >> ${sample}_${db_name}.blast_summary.txt
            echo "333"
            cut -f2 ${sample}_${db_name}.blast.tsv | sed 's/gnl|//; s/|.*//' | sort | uniq -c | sort -rn | head -50 > ${sample}_${db_name}.top_contaminants.txt
            echo "333"
        else
            echo "3333"
            cut -f13 ${sample}_${db_name}.blast.tsv | sort | uniq -c | sort -rn | head -20 >> ${sample}_${db_name}.blast_summary.txt
            echo "3333"
            cut -f13 ${sample}_${db_name}.blast.tsv | sort | uniq -c | sort -rn | head -50 > ${sample}_${db_name}.top_contaminants.txt
            echo "3333"
        fi
    else
        echo "No contamination detected" > ${sample}_${db_name}.top_contaminants.txt
    fi

    echo "4"
    
    echo "BLAST analysis complete for ${sample} vs ${db_name}" >&2

    echo "5"

    rm -f ${sample}_R1.fasta ${sample}_R2.fasta ${sample}_mapped_R1.fastq ${sample}_mapped_R2.fastq

    echo "6"
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
    set -euo pipefail
    
    export BLASTDB=/mnt/sas/Tools/blast_databases
    
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
    path('*')

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
    FILE_COUNT=\$(ls -1 2>/dev/null | wc -l)
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