process BLAST_UNMAPPED_READS_MULTI {
    tag "${sample}_${db_name}"
    publishDir "${params.outdir}/Output/contamination_check/${sample}/unmapped/${db_name}", mode: 'copy'
    cpus { Math.min(8, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'terminate'
    maxForks params.blast_max_parallel

    input:
    tuple val(sample), path(r1_fastq), path(r2_fastq), val(db_name), val(db_path)

    output:
    tuple val(sample), val(db_name), path("${sample}_${db_name}.blast.tsv"), emit: blast_results
    tuple val(sample), val(db_name), path("${sample}_${db_name}.blast_summary.txt"), emit: summary
    tuple val(sample), val(db_name), path("${sample}_${db_name}.top_contaminants.txt"), emit: top_hits
    tuple val(sample), val(db_name), path("${sample}_${db_name}.blast.tsv"), path("${sample}_${db_name}.blast_summary.txt"), emit: for_plotting

    when:
    params.run_contamination_check

    script:
    """
    set -euo pipefail
    
    export BLASTDB=/mnt/sas/Tools/blast_databases
    
    echo "Starting BLAST analysis: ${sample} vs ${db_name}" >&2
    
    # Count reads from both R1 and R2
    READ_COUNT_R1=\$(zcat ${r1_fastq} | wc -l | awk '{print \$1/4}')
    READ_COUNT_R2=\$(zcat ${r2_fastq} | wc -l | awk '{print \$1/4}')
    TOTAL_READ_COUNT=\$(echo "\$READ_COUNT_R1 + \$READ_COUNT_R2" | bc)
    
    echo "BLAST Contamination Analysis for ${sample} against ${db_name}" > ${sample}_${db_name}.blast_summary.txt
    echo "==========================================================" >> ${sample}_${db_name}.blast_summary.txt
    echo "Database: ${db_path}" >> ${sample}_${db_name}.blast_summary.txt
    echo "Unmapped R1 reads: \$READ_COUNT_R1" >> ${sample}_${db_name}.blast_summary.txt
    echo "Unmapped R2 reads: \$READ_COUNT_R2" >> ${sample}_${db_name}.blast_summary.txt
    echo "Total unmapped reads to analyze: \$TOTAL_READ_COUNT" >> ${sample}_${db_name}.blast_summary.txt
    
    # Subsample if needed (sample from EACH file to maintain pairing info)
    MAX_READS_PER_FILE=${params.unmapped_subset_reads}
    if [ "\$READ_COUNT_R1" -gt "\$MAX_READS_PER_FILE" ]; then
        ${params.seqtk_bin} sample -s${params.subset_seed} ${sample}_mapped_R1.fastq.gz ${params.mapped_subset_reads} | gzip > ${sample}_R1.sampled.fastq.gz
        BLAST_INPUT_R1="${sample}_R1.sampled.fastq.gz"
    else
        BLAST_INPUT_R1="${r1_fastq}"
    fi
    
    if [ "\$READ_COUNT_R2" -gt "\$MAX_READS_PER_FILE" ]; then
        ${params.seqtk_bin} sample -s${params.subset_seed} ${sample}_mapped_R2.fastq.gz ${params.mapped_subset_reads} | gzip > ${sample}_R2.sampled.fastq.gz
        BLAST_INPUT_R2="${sample}_R2.sampled.fastq.gz"
    else
        BLAST_INPUT_R2="${r2_fastq}"
    fi
    
    # Convert R1 to FASTA
    ${params.reformat_bin} in=\$BLAST_INPUT_R1 out=${sample}_R1.fasta
    
    # Convert R2 to FASTA
    ${params.reformat_bin} in=\$BLAST_INPUT_R2 out=${sample}_R2.fasta
    
    echo "Running BLAST on R1 against ${db_name}..." >&2
    # BLAST R1
    ${params.blastn_bin} \\
        -query ${sample}_R1.fasta \\
        -db ${db_path} \\
        -out ${sample}_${db_name}_R1.blast.tsv \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
        -evalue ${params.contamination_evalue} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 5 \\
        -max_hsps 1
    
    echo "Running BLAST on R2 against ${db_name}..." >&2
    # BLAST R2
    ${params.blastn_bin} \\
        -query ${sample}_R2.fasta \\
        -db ${db_path} \\
        -out ${sample}_${db_name}_R2.blast.tsv \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
        -evalue ${params.contamination_evalue} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 5 \\
        -max_hsps 1
    
    # Combine results
    cat ${sample}_${db_name}_R1.blast.tsv ${sample}_${db_name}_R2.blast.tsv > ${sample}_${db_name}.blast.tsv
    
    BLAST_HITS=\$(wc -l < ${sample}_${db_name}.blast.tsv)
    echo "Total BLAST hits (R1+R2): \$BLAST_HITS" >> ${sample}_${db_name}.blast_summary.txt
    
    if [ "\$BLAST_HITS" -gt 0 ]; then
        echo "" >> ${sample}_${db_name}.blast_summary.txt
        echo "TOP CONTAMINATING SEQUENCES/STRAINS:" >> ${sample}_${db_name}.blast_summary.txt
        
        if [[ "${db_name}" == *"mus_strain"* ]]; then
            cut -f2 ${sample}_${db_name}.blast.tsv | sed 's/gnl|//; s/|.*//' | sort | uniq -c | sort -rn | head -20 >> ${sample}_${db_name}.blast_summary.txt
            cut -f2 ${sample}_${db_name}.blast.tsv | sed 's/gnl|//; s/|.*//' | sort | uniq -c | sort -rn | head -50 > ${sample}_${db_name}.top_contaminants.txt
        else
            cut -f13 ${sample}_${db_name}.blast.tsv | sort | uniq -c | sort -rn | head -20 >> ${sample}_${db_name}.blast_summary.txt
            cut -f13 ${sample}_${db_name}.blast.tsv | sort | uniq -c | sort -rn | head -50 > ${sample}_${db_name}.top_contaminants.txt
        fi
    else
        echo "No contamination detected" > ${sample}_${db_name}.top_contaminants.txt
    fi
    
    echo "BLAST analysis complete for ${sample} vs ${db_name}" >&2
    rm -f ${sample}_R1.fasta ${sample}_R2.fasta ${sample}_R1.sampled.fastq.gz ${sample}_R2.sampled.fastq.gz
    """
}