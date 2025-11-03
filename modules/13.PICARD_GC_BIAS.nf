process PICARD_GC_BIAS {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/picard_gc_bias/${sample}", mode: 'copy'
    memory { params.max_mem }
    errorStrategy 'retry'; maxRetries 2

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

    java -jar ${params.picard_jar} CollectGcBiasMetrics \\
        I=${bam} \\
        O=${sample}.gc_bias_metrics.txt \\
        SUMMARY_OUTPUT=${sample}.gc_bias_summary.txt \\
        CHART_OUTPUT=${sample}.gc_bias.pdf \\
        R=${ref} \\
        VALIDATION_STRINGENCY=LENIENT

    java -jar ${params.picard_jar} CollectMultipleMetrics \\
        I=${bam} \\
        O=${sample}_allmetrics \\
        R=${ref}
    """
}