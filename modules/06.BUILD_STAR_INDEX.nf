process BUILD_STAR_INDEX {
    tag "${strain}"
    publishDir "${params.star_index_dir}", mode: 'copy', pattern: "${strain}"
    cpus { Math.min(params.star_index_threads, params.max_cpus as int) }
    memory { params.star_index_mem }
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(strain), path(fasta), path(gtf)

    output:
    tuple val(strain), path("${strain}"), emit: star_index

    script:
    """
    set -euo pipefail
    
    echo "Building STAR index for strain: ${strain}"
    echo "This may take 30-60 minutes..."
    
    mkdir -p ${strain}
    
    if [[ ${fasta} == *.gz ]]; then
        zcat ${fasta} > ${strain}_genome.fa
        FASTA_FILE="${strain}_genome.fa"
    else
        FASTA_FILE="${fasta}"
    fi
    
    if [[ ${gtf} == *.gz ]]; then
        zcat ${gtf} > ${strain}_annotation.gtf
        GTF_FILE="${strain}_annotation.gtf"
    else
        GTF_FILE="${gtf}"
    fi
    
    ${params.star_bin} \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${strain} \\
        --genomeFastaFiles \${FASTA_FILE} \\
        --sjdbGTFfile \${GTF_FILE} \\
        --sjdbOverhang 149
    
    sed 's///g'

    echo "STAR index built successfully for ${strain}"
    """
}