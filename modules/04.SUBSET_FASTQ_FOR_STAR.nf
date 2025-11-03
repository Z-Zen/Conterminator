process SUBSET_FASTQ_FOR_STAR {
    tag "${sample}"
    publishDir "${params.outdir}/Input/subsampled_fastq_star/${sample}", mode: 'copy'
    cpus { Math.min(2, params.max_cpus as int) }
    memory '4 GB'
    errorStrategy 'retry'
    maxRetries 2
    maxForks params.max_parallel_samples

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
    
    echo "Subsampling FASTQ for STAR: ${sample}" > ${sample}_star_subsample_stats.txt
    echo "======================================" >> ${sample}_star_subsample_stats.txt
    echo "Target: ${params.subset_star_reads} reads per file" >> ${sample}_star_subsample_stats.txt
    echo "Random seed: ${params.subset_seed}" >> ${sample}_star_subsample_stats.txt
    echo "" >> ${sample}_star_subsample_stats.txt
    
    ORIG_READS_R1=\$(zcat ${fastq1} | wc -l | awk '{print \$1/4}')
    ORIG_READS_R2=\$(zcat ${fastq2} | wc -l | awk '{print \$1/4}')
    
    echo "Original reads R1: \$ORIG_READS_R1" >> ${sample}_star_subsample_stats.txt
    echo "Original reads R2: \$ORIG_READS_R2" >> ${sample}_star_subsample_stats.txt
    
    ${params.seqtk_bin} sample -s${params.subset_seed} ${fastq1} ${params.subset_star_reads} | gzip > ${sample}_star_subsampled_1.fastq.gz
    ${params.seqtk_bin} sample -s${params.subset_seed} ${fastq2} ${params.subset_star_reads} | gzip > ${sample}_star_subsampled_2.fastq.gz
    
    SUB_READS_R1=\$(zcat ${sample}_star_subsampled_1.fastq.gz | wc -l | awk '{print \$1/4}')
    SUB_READS_R2=\$(zcat ${sample}_star_subsampled_2.fastq.gz | wc -l | awk '{print \$1/4}')
    
    echo "Subsampled reads R1: \$SUB_READS_R1" >> ${sample}_star_subsample_stats.txt
    echo "Subsampled reads R2: \$SUB_READS_R2" >> ${sample}_star_subsample_stats.txt
    echo "Purpose: STAR alignment" >> ${sample}_star_subsample_stats.txt
    """
}