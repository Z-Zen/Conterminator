process DECONTAMINER_STEP3 {
    tag "${sample} ${reads_type}"
    publishDir "${params.outdir}/Output/decontaminer/${sample}/", mode: 'copy'
    cpus { Math.min(2, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'ignore'

    input:
    tuple val(sample), val(reads_type), path(filtered_dir)

    output:
    path "bacteria_reports/**", optional: true, emit: bacteria_reports
    path "fungi_reports/**", optional: true, emit: fungi_reports
    path "virus_reports/**", optional: true, emit: virus_reports
    path "HTML_REPORTS/**", optional: true, emit: html_reports

    when:
    params.run_decontaminer

    script:
    def plots_flag = params.decontaminer_generate_plots ? "-P y" : ""
    """
    set -euo pipefail

    export PATH="/opt/R/4.5.1/bin:\$PATH"

    if [ -d "${filtered_dir}/RESULTS/BACTERIA/COLLECTED_INFO" ]; then
        bash ${params.decontaminer_dir}/shell_scripts/collectInfo.sh \\
            -i \$(pwd)/${filtered_dir}/RESULTS/BACTERIA/COLLECTED_INFO \\
            -t ${params.decontaminer_match_threshold} \\
            -V O \\
            ${plots_flag} || echo "Bacteria collection completed with warnings"

        mkdir -p bacteria_reports
        # Copy COLLECTED_INFO outputs
        cp -r ${filtered_dir}/RESULTS/BACTERIA/COLLECTED_INFO/* bacteria_reports/ 2>/dev/null || true
    fi

    if [ -d "${filtered_dir}/RESULTS/FUNGI/COLLECTED_INFO" ]; then
        bash ${params.decontaminer_dir}/shell_scripts/collectInfo.sh \\
            -i \$(pwd)/${filtered_dir}/RESULTS/FUNGI/COLLECTED_INFO \\
            -t ${params.decontaminer_match_threshold} \\
            -V O \\
            ${plots_flag} || echo "Fungi collection completed with warnings"

        mkdir -p fungi_reports
        cp -r ${filtered_dir}/RESULTS/FUNGI/COLLECTED_INFO/* fungi_reports/ 2>/dev/null || true
    fi

    if [ -d "${filtered_dir}/RESULTS/VIRUSES/COLLECTED_INFO" ]; then
        bash ${params.decontaminer_dir}/shell_scripts/collectInfo.sh \\
            -i \$(pwd)/${filtered_dir}/RESULTS/VIRUSES/COLLECTED_INFO \\
            -t ${params.decontaminer_match_threshold} \\
            -V V \\
            ${plots_flag} || echo "Virus collection completed with warnings"

        mkdir -p virus_reports
        cp -r ${filtered_dir}/RESULTS/VIRUSES/COLLECTED_INFO/* virus_reports/ 2>/dev/null || true
    fi

    if [ -d "${filtered_dir}/RESULTS/HTML_REPORTS" ]; then
        echo "Copying shared HTML_REPORTS directory..."
        cp -r ${filtered_dir}/RESULTS/HTML_REPORTS . || echo "Warning: Could not copy HTML_REPORTS"
    else
        echo "No HTML_REPORTS directory found - this is normal if no valid contamination was detected"
    fi
    """
}