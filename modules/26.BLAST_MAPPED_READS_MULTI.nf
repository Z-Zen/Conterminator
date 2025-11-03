process BLAST_MAPPED_READS_MULTI {
    tag "${sample}_${db_name}"
    publishDir "${params.outdir}/Output/contamination_check/${sample}/mapped/${db_name}", mode: 'copy'
    cpus { Math.min(8, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'terminate'
    maxForks params.blast_max_parallel

    input:
    tuple val(sample), val(strain), path(bam), path(bai), val(db_name), val(db_path)

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
    echo "Extracting mapped reads from BAM: ${bam}" >&2
    
    # Extract mapped reads from BAM to FASTQ
    ${params.samtools_bin} fastq \\
        -1 ${sample}_mapped_R1.fastq \\
        -2 ${sample}_mapped_R2.fastq \\
        -0 /dev/null \\
        -s /dev/null \\
        -N \\
        ${bam}
    
    ${params.reformat_bin} in=${sample}_mapped_R1.fastq out=${sample}_R1.fasta
    ${params.reformat_bin} in=${sample}_mapped_R2.fastq out=${sample}_R2.fasta

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
    
    echo "1"

    # Combine results
    cat ${sample}_${db_name}_R1.blast.tsv ${sample}_${db_name}_R2.blast.tsv > ${sample}_${db_name}.blast.tsv

    echo "2"
    
    BLAST_HITS=\$(wc -l < ${sample}_${db_name}.blast.tsv)
    echo "Total BLAST hits (R1+R2): \$BLAST_HITS" >> ${sample}_${db_name}.blast_summary.txt

    echo "3"
    
    if [ "\$BLAST_HITS" -gt 0 ]; then
        echo "" >> ${sample}_${db_name}.blast_summary.txt
        echo "TOP CONTAMINATING SEQUENCES/STRAINS:" >> ${sample}_${db_name}.blast_summary.txt

        echo "33"
        
        if [[ "${db_name}" == *"mus_strain"* ]]; then
            echo "333"
            cut -f2 ${sample}_${db_name}.blast.tsv | sed 's/gnl|//; s/|.*//' | sort | uniq -c | sort -rn | head -20 >> ${sample}_${db_name}.blast_summary.txt
            echo "333"
            cut -f2 ${sample}_${db_name}.blast.tsv | sed 's/gnl|//; s/|.*//' | sort | uniq -c | sort -rn | head -50 > ${sample}_${db_name}.top_contaminants.txt
            echo "333"
        else
            echo "3333"
            cut -f13 ${sample}_${db_name}.blast.tsv | sort | uniq -c | sort -rn | head -20 >> ${sample}_${db_name}.blast_summary.txt
            echo "3333"
            cut -f13 ${sample}_${db_name}.blast.tsv | sort | uniq -c | sort -rn | head -50 > ${sample}_${db_name}.top_contaminants.txt
            echo "3333"
        fi
    else
        echo "No contamination detected" > ${sample}_${db_name}.top_contaminants.txt
    fi

    echo "4"
    
    echo "BLAST analysis complete for ${sample} vs ${db_name}" >&2

    echo "5"

    rm -f ${sample}_R1.fasta ${sample}_R2.fasta ${sample}_mapped_R1.fastq ${sample}_mapped_R2.fastq

    echo "6"
    """
}