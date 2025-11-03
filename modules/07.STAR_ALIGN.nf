process STAR_ALIGN {
    tag "${sample}_${strain}"
    publishDir "${params.outdir}/Temporary/star_alignment/${sample}/${strain}", mode: 'copy'
    cpus { Math.min(params.star_threads, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'retry'
    maxRetries 2
    maxForks params.max_parallel_samples

    input:
    tuple val(sample), path(fastq1), path(fastq2), val(strain), path(genome_dir)

    output:
    tuple val(sample), val(strain), path("${sample}_${strain}_Aligned.sortedByCoord.out.bam"), path("${sample}_${strain}_Aligned.sortedByCoord.out.bam.bai"), emit: aligned_bam
    tuple val(sample), val(strain), path("${sample}_${strain}_Unmapped.out.mate1"), path("${sample}_${strain}_Unmapped.out.mate2"), emit: unmapped_fastq
    path "${sample}_${strain}_Log.final.out", emit: log
    path "${sample}_${strain}_Log.out", emit: log_out
    path "${sample}_${strain}_Log.progress.out", emit: log_progress
    path "${sample}_${strain}_ReadsPerGene.out.tab", optional: true, emit: gene_counts
    path "${sample}_${strain}_Aligned.toTranscriptome.out.bam", optional: true, emit: transcriptome_bam
    path "${sample}_${strain}_star_alignment_stats.txt", emit: stats

    when:
    params.run_star_alignment

    script:
    """
    set -euo pipefail
    
    echo "STAR Alignment for ${sample} against ${strain}" > ${sample}_${strain}_star_alignment_stats.txt
    echo "===========================================" >> ${sample}_${strain}_star_alignment_stats.txt
    echo "Genome directory: ${genome_dir}" >> ${sample}_${strain}_star_alignment_stats.txt
    echo "Threads: ${task.cpus}" >> ${sample}_${strain}_star_alignment_stats.txt
    echo "" >> ${sample}_${strain}_star_alignment_stats.txt
    
    ${params.star_bin} \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${genome_dir} \\
        --quantMode GeneCounts TranscriptomeSAM \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesIn ${fastq1} ${fastq2} \\
        --outFileNamePrefix ${sample}_${strain}_ \\
        --outReadsUnmapped Fastx
    
    ${params.samtools_bin} index ${sample}_${strain}_Aligned.sortedByCoord.out.bam

    echo "" >> ${sample}_${strain}_star_alignment_stats.txt
    echo "ALIGNMENT SUMMARY:" >> ${sample}_${strain}_star_alignment_stats.txt
    cat ${sample}_${strain}_Log.final.out >> ${sample}_${strain}_star_alignment_stats.txt
    """
}