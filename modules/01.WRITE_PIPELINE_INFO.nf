// Write pipeline manifest info to output directory
process WRITE_PIPELINE_INFO {
    publishDir "${params.outdir}/", mode: 'copy'

    output:
    path "pipeline_info.txt"

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
FASTQ QC:     ${params.subset_for_fastq_qc} (${params.subset_fastq_qc_reads} reads)
STAR:         ${params.subset_for_star} (${params.subset_star_reads} reads)
BAM QC:       ${params.subset_bam_for_qc} (${params.bam_qc_subset_mapped} reads)
Mapped BLAST: ${params.subset_mapped_for_blast} (${params.mapped_subset_reads} reads)
Unmapped BLAST: ${params.subset_unmapped_for_blast} (${params.unmapped_subset_reads} reads)

========================================
EOF
    """
}