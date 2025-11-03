process MAPINSIGHTS {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/mapinsights/${sample}", mode: 'copy'
    cpus { Math.min(8, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'ignore'

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
    
    if [ -f "${params.mapinsights_bin}" ] && [ -x "${params.mapinsights_bin}" ]; then
        ${params.mapinsights_bin} bamqc \\
            -r ${ref} \\
            -i ${bam} \\
            -o ${sample}_mapinsights \\
            ${params.mapinsights_opts}
    else
        mkdir -p ${sample}_mapinsights
        ${params.samtools_bin} flagstat ${bam} > ${sample}_mapinsights/${sample}.flagstat
        ${params.samtools_bin} stats ${bam} > ${sample}_mapinsights/${sample}.stats
        ${params.samtools_bin} idxstats ${bam} > ${sample}_mapinsights/${sample}.idxstats
    fi
    """
}