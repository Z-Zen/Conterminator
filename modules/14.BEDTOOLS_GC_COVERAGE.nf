process BEDTOOLS_GC_COVERAGE {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/bedtools_gc/${sample}/", mode: 'copy'
    cpus { Math.min(4, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'retry'; maxRetries 2

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

    ${params.bedtools_bin} makewindows -b ${genome_bed} -w ${params.windowsize} -s ${params.window_step} > ${sample}.windows.bed

    ${params.bedtools_bin} nuc -fi ${ref} -bed ${sample}.windows.bed | \\
      awk 'BEGIN{OFS="\\t"} NR==1{print "chrom","start","end","pct_at","pct_gc"} NR>1{print \$1,\$2,\$3,\$4,\$5}' \\
      > ${sample}.windows.gc.tsv

    ${params.bedtools_bin} coverage -a ${sample}.windows.bed -b ${bam} -mean | \\
      paste ${sample}.windows.gc.tsv - | \\
      awk 'BEGIN{OFS="\\t"} NR==1{print \$0,"mean_coverage"} NR>1{print}' \\
      > ${sample}.gc_coverage.tsv

    # Generate GC content plots
    Rscript ${params.gc_plot_script} ${sample}.gc_coverage.tsv ${sample}

    gzip ${sample}.gc_coverage.tsv ${sample}.windows.gc.tsv ${sample}.windows.bed
    """
}