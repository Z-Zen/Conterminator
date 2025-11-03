process DECONTAMINER_STEP1_STAR_UNMAPPED {
    tag "${sample} ${reads_type}"
    publishDir "${params.outdir}/Temporary/decontaminer/${sample}/step1_decontaminer_unmapped", mode: 'copy'
    cpus { Math.min(8, params.max_cpus as int) }
    memory { params.max_mem }

    input:
    tuple val(sample), path(unmapped_r1), path(unmapped_r2), val(reads_type)

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
    
    # Add /1 and /2 suffixes to read names (only to header lines)
    awk 'NR%4==1 {sub(/^@/, ""); print "@"\$0"/1"} NR%4!=1' ${unmapped_r1} > input_fastq/\${SAMPLE}_unmapped_1.fq
    awk 'NR%4==1 {sub(/^@/, ""); print "@"\$0"/2"} NR%4!=1' ${unmapped_r2} > input_fastq/\${SAMPLE}_unmapped_2.fq
    
    echo "Running DecontaMiner on UNMAPPED reads...."
    bash ${params.decontaminer_dir}/shell_scripts/decontaMiner.sh \\
        -i \$(pwd)/input_fastq \\
        -o \$(pwd)/decontaminer_output \\
        -c ${params.decontaminer_config} \\
        -F fastq \\
        -s ${params.decontaminer_pairing} \\
        ${quality_filter} \\
        ${ribo_filter} \\
        ${organisms} || echo "DecontaMiner step 1 (unmapped) completed with warnings"
    """
}