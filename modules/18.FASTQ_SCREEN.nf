process FASTQ_SCREEN {
    tag "${sample}"
    publishDir "${params.outdir}/Output/fastq_screen/${sample}", mode: 'copy'
    cpus { Math.min(params.fastq_screen_threads as int, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'retry'
    maxRetries 2

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
    ${params.fastq_screen} \\
        --conf ${params.fastq_screen_conf} \\
        --threads ${task.cpus} \\
        --outdir . \\
        --force \\
        ${fastq1} ${fastq2}
    """
}