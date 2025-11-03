process DECONTAMINER_STEP2 {
    tag "${sample} ${reads_type}"
    publishDir "${params.outdir}/Temporary/decontaminer/${sample}/step2_filtering", mode: 'copy'
    cpus { Math.min(4, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'ignore'

    input:
    tuple val(sample), val(reads_type), path(output_dir)

    output:
    tuple val(sample), val(reads_type), path("${output_dir}"), emit: filtered_dir

    when:
    params.run_decontaminer

    script:
    """
    set -euo pipefail
    
    if [ -d "${output_dir}/RESULTS/BACTERIA" ]; then
        bash ${params.decontaminer_dir}/shell_scripts/filterBlastInfo.sh \\
            -i \$(pwd)/${output_dir}/RESULTS/BACTERIA/ \\
            -s ${params.decontaminer_pairing} \\
            -g ${params.decontaminer_gap} \\
            -m ${params.decontaminer_mismatch} \\
            -l ${params.decontaminer_match_len} \\
            -V O || echo "Bacteria/Fungi filtering completed with warnings"
    fi

    if [ -d "${output_dir}/RESULTS/FUNGI" ]; then
        bash ${params.decontaminer_dir}/shell_scripts/filterBlastInfo.sh \\
            -i \$(pwd)/${output_dir}/RESULTS/FUNGI/ \\
            -s ${params.decontaminer_pairing} \\
            -g ${params.decontaminer_gap} \\
            -m ${params.decontaminer_mismatch} \\
            -l ${params.decontaminer_match_len} \\
            -V O || echo "Bacteria/Fungi filtering completed with warnings"
    fi
    
    if [ -d "${output_dir}/RESULTS/VIRUSES" ]; then
        bash ${params.decontaminer_dir}/shell_scripts/filterBlastInfo.sh \\
            -i \$(pwd)/${output_dir}/RESULTS/VIRUSES \\
            -s ${params.decontaminer_pairing} \\
            -g ${params.decontaminer_gap} \\
            -m ${params.decontaminer_mismatch} \\
            -l ${params.decontaminer_match_len} \\
            -V V || echo "Virus filtering completed with warnings"
    fi
    
    echo "Filtering step completed for ${sample}"
    """
}