process DECONTAMINER_STEP1_STAR_MAPPED {
    tag "${sample} ${reads_type}"
    publishDir "${params.outdir}/Temporary/decontaminer/${sample}/step1_decontaminer_mapped", mode: 'copy'
    cpus { Math.min(8, params.max_cpus as int) }
    memory { params.max_mem }

    input:
    tuple val(sample), path(mapped_bam), path(mapped_bai), val(reads_type)

    output:
    tuple val(sample), val(reads_type), path("decontaminer_output"), emit: output_dir
    path "decontaminer_output/RESULTS/**", optional: true

    when:
    params.run_decontaminer

    script:
    def organisms = ""
    if (params.decontaminer_organisms.contains('b')) organisms += " -b"
    if (params.decontaminer_organisms.contains('f')) organisms += " -f"
    if (params.decontaminer_organisms.contains('v')) organisms += " -v"
    
    def quality_filter = params.decontaminer_quality_filter ? "-Q ${params.decontaminer_quality_filter}" : ""
    def ribo_filter = params.decontaminer_ribo_filter ? "-R ${params.decontaminer_ribo_filter}" : ""
    """
    set -euo pipefail
    
    SAMPLE="${sample}"
    mkdir -p input_fastq

    ${params.samtools_bin} fastq -1 input_fastq/\${SAMPLE}_mapped_1.fq -2 input_fastq/\${SAMPLE}_mapped_2.fq -0 /dev/null -s /dev/null -N ${mapped_bam}
    
    echo "Running DecontaMiner on MAPPED reads...."
    bash ${params.decontaminer_dir}/shell_scripts/decontaMiner.sh \\
        -i \$(pwd)/input_fastq \\
        -o \$(pwd)/decontaminer_output \\
        -c ${params.decontaminer_config} \\
        -F fastq \\
        -s ${params.decontaminer_pairing} \\
        ${quality_filter} \\
        ${ribo_filter} \\
        ${organisms} || echo "DecontaMiner step 1 (mapped) completed with warnings"
    
    mkdir -p decontaminer_output
    """
}