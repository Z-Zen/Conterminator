process FASTQC_MAPPED_BAM {
    tag "${sample} ${strain}"
    publishDir "${params.outdir}/Output/fastqc/${sample}/mapped_bam", mode: 'copy'
    cpus { Math.min(4, params.max_cpus as int) }
    memory '4 GB'
    errorStrategy 'retry'
    maxRetries 2
    maxForks params.max_parallel_samples

    input:
    tuple val(sample), val(strain), path(bam), path(bai)

    output:
    tuple path("*_fastqc.html"), path("*_fastqc.zip"), emit: results

    when:
    params.run_fastqc && params.run_contamination_check

    script:
    """
    set -euo pipefail
    
    ${params.fastqc_bin} \\
        --outdir . \\
        --threads ${task.cpus} \\
        --format bam \\
        ${bam}
    """
}