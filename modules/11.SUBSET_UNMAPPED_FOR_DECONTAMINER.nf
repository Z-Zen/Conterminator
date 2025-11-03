process SUBSET_UNMAPPED_FOR_DECONTAMINER {
    tag "${sample}_${strain}"
    publishDir "${params.outdir}/Input/subsampled_unmapped_decontaminer/${sample}/${strain}", mode: 'copy'
    cpus { Math.min(2, params.max_cpus as int) }
    memory '4 GB'
    errorStrategy 'retry'
    maxRetries 2
    maxForks params.max_parallel_samples

    input:
    tuple val(sample), val(strain), path(unmapped_r1), path(unmapped_r2)

    output:
    tuple val(sample), val(strain), path("${sample}_${strain}_unmapped_decon_subset_1.fq"), path("${sample}_${strain}_unmapped_decon_subset_2.fq"), emit: subsampled_fastq
    path "${sample}_${strain}_unmapped_decon_subset_stats.txt", emit: stats

    when:
    params.subset_unmapped_for_decontaminer

    script:
    """
    set -euo pipefail
    
    echo "Subsampling Unmapped for DecontaMiner: ${sample}_${strain}" > ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "=========================================================" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "Target: ${params.unmapped_subset_reads} reads per file" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    
    TOTAL_R1=\$(cat ${unmapped_r1} | wc -l | awk '{print \$1/4}')
    TOTAL_R2=\$(cat ${unmapped_r2} | wc -l | awk '{print \$1/4}')
    
    echo "Total unmapped R1: \$TOTAL_R1" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "Total unmapped R2: \$TOTAL_R2" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    
    # Subsample and add /1 suffix to R1 reads
    if [ "\$TOTAL_R1" -gt "${params.unmapped_subset_reads}" ]; then
        ${params.seqtk_bin} sample -s${params.subset_seed} ${unmapped_r1} ${params.unmapped_subset_reads} | \\
            sed '1~4 s/\$/\\/1/' > ${sample}_${strain}_unmapped_decon_subset_1.fq
        echo "Subsampled R1 to ${params.unmapped_subset_reads} reads" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    else
        sed '1~4 s/\$/\\/1/' ${unmapped_r1} > ${sample}_${strain}_unmapped_decon_subset_1.fq
        echo "Using all \$TOTAL_R1 unmapped R1 reads" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    fi
    
    # Subsample and add /2 suffix to R2 reads
    if [ "\$TOTAL_R2" -gt "${params.unmapped_subset_reads}" ]; then
        ${params.seqtk_bin} sample -s${params.subset_seed} ${unmapped_r2} ${params.unmapped_subset_reads} | \\
            sed '1~4 s/\$/\\/2/' > ${sample}_${strain}_unmapped_decon_subset_2.fq
        echo "Subsampled R2 to ${params.unmapped_subset_reads} reads" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    else
        sed '1~4 s/\$/\\/2/' ${unmapped_r2} > ${sample}_${strain}_unmapped_decon_subset_2.fq
        echo "Using all \$TOTAL_R2 unmapped R2 reads" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    fi
    
    FINAL_R1=\$(cat ${sample}_${strain}_unmapped_decon_subset_1.fq | wc -l | awk '{print \$1/4}')
    FINAL_R2=\$(cat ${sample}_${strain}_unmapped_decon_subset_2.fq | wc -l | awk '{print \$1/4}')
    echo "Final R1 reads: \$FINAL_R1" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "Final R2 reads: \$FINAL_R2" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    echo "Purpose: DecontaMiner analysis" >> ${sample}_${strain}_unmapped_decon_subset_stats.txt
    """
}