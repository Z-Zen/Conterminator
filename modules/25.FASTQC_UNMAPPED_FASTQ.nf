process FASTQC_UNMAPPED_FASTQ {
    tag "${sample} ${strain}"
    publishDir "${params.outdir}/Output/fastqc/${sample}/unmapped_fastq", mode: 'copy'
    cpus { Math.min(4, params.max_cpus as int) }
    memory '4 GB'
    errorStrategy 'retry'
    maxRetries 2
    maxForks params.max_parallel_samples

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
    ${params.fastqc_bin} \\
        --outdir . \\
        --threads ${task.cpus} \\
        ${fastq_r1}
    
    # Run FastQC on R2
    ${params.fastqc_bin} \\
        --outdir . \\
        --threads ${task.cpus} \\
        ${fastq_r2}
    """
}