process SUBSET_UNMAPPED_FOR_BLAST {
    tag "${sample}_${strain}"
    publishDir "${params.outdir}/Input/subsampled_unmapped_blast/${sample}/${strain}", mode: 'copy'
    cpus { Math.min(2, params.max_cpus as int) }
    memory '4 GB'
    errorStrategy 'retry'
    maxRetries 2
    maxForks params.max_parallel_samples

    input:
    tuple val(sample), val(strain), path(unmapped_r1), path(unmapped_r2)

    output:
    tuple val(sample), val(strain), path("${sample}_${strain}_unmapped_blast_R1.fastq.gz"), path("${sample}_${strain}_unmapped_blast_R2.fastq.gz"), emit: subsampled_fastq
    path "${sample}_${strain}_unmapped_blast_subset_stats.txt", emit: stats

    when:
    params.subset_unmapped_for_blast

    script:
    """
    set -euo pipefail
    
    echo "Subsampling Unmapped for BLAST: ${sample}_${strain}" > ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "===================================================" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "Target: ${params.unmapped_subset_reads} reads per file" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    
    TOTAL_R1=\$(cat ${unmapped_r1} | wc -l | awk '{print \$1/4}')
    TOTAL_R2=\$(cat ${unmapped_r2} | wc -l | awk '{print \$1/4}')
    
    echo "Total unmapped R1: \$TOTAL_R1" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "Total unmapped R2: \$TOTAL_R2" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    
    # Subsample R1
    if [ "\$TOTAL_R1" -gt "${params.unmapped_subset_reads}" ]; then
        ${params.seqtk_bin} sample -s${params.subset_seed} ${unmapped_r1} ${params.unmapped_subset_reads} | gzip > ${sample}_${strain}_unmapped_blast_R1.fastq.gz
        echo "R1 subsampled to ${params.unmapped_subset_reads} reads" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    else
        cat ${unmapped_r1} | gzip > ${sample}_${strain}_unmapped_blast_R1.fastq.gz
        echo "R1 using all \$TOTAL_R1 reads (less than target)" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    fi
    
    # Subsample R2
    if [ "\$TOTAL_R2" -gt "${params.unmapped_subset_reads}" ]; then
        ${params.seqtk_bin} sample -s${params.subset_seed} ${unmapped_r2} ${params.unmapped_subset_reads} | gzip > ${sample}_${strain}_unmapped_blast_R2.fastq.gz
        echo "R2 subsampled to ${params.unmapped_subset_reads} reads" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    else
        cat ${unmapped_r2} | gzip > ${sample}_${strain}_unmapped_blast_R2.fastq.gz
        echo "R2 using all \$TOTAL_R2 reads (less than target)" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    fi
    
    FINAL_R1=\$(zcat ${sample}_${strain}_unmapped_blast_R1.fastq.gz | wc -l | awk '{print \$1/4}')
    FINAL_R2=\$(zcat ${sample}_${strain}_unmapped_blast_R2.fastq.gz | wc -l | awk '{print \$1/4}')
    TOTAL_FINAL=\$(echo "\$FINAL_R1 + \$FINAL_R2" | bc)
    
    echo "Final R1 reads for BLAST: \$FINAL_R1" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "Final R2 reads for BLAST: \$FINAL_R2" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "Total final reads for BLAST: \$TOTAL_FINAL" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    echo "Purpose: BLAST contamination screening" >> ${sample}_${strain}_unmapped_blast_subset_stats.txt
    """
}