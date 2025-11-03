process DEEPTOOLS_GC {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/deeptools_gc_bias/${sample}/", mode: 'copy'
    cpus { Math.min(4, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'retry'; maxRetries 2

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
    
    python3 ${params.deeptools_gcbias} \\
        --bamfile ${bam} \\
        --genome ${ref2bit} \\
        --effectiveGenomeSize \$EFFECTIVE_SIZE \\
        --GCbiasFrequenciesFile ${sample}.gcBias.freq.txt \\
        --biasPlot ${sample}.gcBias.plot.pdf \\
        --numberOfProcessors ${task.cpus} \\
        ${blacklist}
    """
}