process SUBSET_BAM_FOR_QC {
    tag "${sample}_${strain}"
    publishDir "${params.outdir}/Input/subsampled_bams/${sample}/${strain}", mode: 'copy'
    cpus { Math.min(4, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'retry'
    maxRetries 2
    maxForks params.max_parallel_samples

    input:
    tuple val(sample), val(strain), path(bam), path(bai)

    output:
    tuple val(sample), val(strain), path("${sample}_${strain}_subset.bam"), path("${sample}_${strain}_subset.bam.bai"), emit: subsampled_bam
    path "${sample}_${strain}_bam_subset_stats.txt", emit: stats

    when:
    params.subset_bam_for_qc

    script:
    """
    set -euo pipefail
    
    echo "BAM Subsampling for QC: ${sample}_${strain}" > ${sample}_${strain}_bam_subset_stats.txt
    echo "==========================================" >> ${sample}_${strain}_bam_subset_stats.txt
    echo "Target: ${params.bam_qc_subset_mapped} mapped reads" >> ${sample}_${strain}_bam_subset_stats.txt
    echo "" >> ${sample}_${strain}_bam_subset_stats.txt
    
    # Count total mapped reads
    TOTAL_MAPPED=\$(${params.samtools_bin} view -c -F 4 ${bam})
    echo "Total mapped reads: \$TOTAL_MAPPED" >> ${sample}_${strain}_bam_subset_stats.txt

    if [ "\$TOTAL_MAPPED" -gt "${params.bam_qc_subset_mapped}" ]; then
        FRACTION=\$(awk -v seed=${params.subset_seed} -v target=${params.bam_qc_subset_mapped} -v total=\$TOTAL_MAPPED 'BEGIN {frac=target/total*10000; printf "%d.%04d", seed, frac}')
        echo "Target reads: ${params.bam_qc_subset_mapped}" >> ${sample}_${strain}_bam_subset_stats.txt
        echo "Sampling fraction: \$FRACTION" >> ${sample}_${strain}_bam_subset_stats.txt
        
        # Subsample mapped reads only
        ${params.samtools_bin} view -b -s \${FRACTION} -F 4 ${bam} > temp_mapped.bam
        
        # Verify actual count
        ACTUAL_SAMPLED=\$(${params.samtools_bin} view -c temp_mapped.bam)
        echo "Actually sampled: \$ACTUAL_SAMPLED reads" >> ${sample}_${strain}_bam_subset_stats.txt
    else
        echo "Using all mapped reads (total: \$TOTAL_MAPPED < target: ${params.bam_qc_subset_mapped})" >> ${sample}_${strain}_bam_subset_stats.txt
        ${params.samtools_bin} view -b -F 4 ${bam} > temp_mapped.bam
    fi
    
    # Sort and index
    ${params.samtools_bin} sort -@ ${task.cpus} -o ${sample}_${strain}_subset.bam temp_mapped.bam
    ${params.samtools_bin} index ${sample}_${strain}_subset.bam
    
    # Report final count
    SUBSET_COUNT=\$(${params.samtools_bin} view -c ${sample}_${strain}_subset.bam)
    echo "Subsampled reads: \$SUBSET_COUNT" >> ${sample}_${strain}_bam_subset_stats.txt
    echo "Purpose: QC tools (DeepTools, Picard, BEDTools, Qualimap, Mapinsights)" >> ${sample}_${strain}_bam_subset_stats.txt
    
    rm -f temp_mapped.bam
    """
}
