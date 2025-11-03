process FASTQC_INPUT_FASTQ {
    tag "${sample}"
    publishDir "${params.outdir}/Output/fastqc/${sample}/input_fastq", mode: 'copy'
    cpus { Math.min(4, params.max_cpus as int) }
    memory '4 GB'
    errorStrategy 'retry'
    maxRetries 2
    maxForks params.max_parallel_samples

    input:
    tuple val(sample), path(fastq1), path(fastq2)

    output:
    tuple path("*_fastqc.html"), path("*_fastqc.zip"), emit: results

    when:
    params.run_fastqc

    script:
    """
    set -euo pipefail
    
    ${params.fastqc_bin} \\
        --outdir . \\
        --threads ${task.cpus} \\
        ${fastq1} ${fastq2}
    """
}