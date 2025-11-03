process MULTIQC {
    publishDir "${params.outdir}/Output/multiqc", mode: 'copy'
    cpus { Math.min(4, params.max_cpus as int) }
    memory '8 GB'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data, optional: true
    path "multiqc_plots", optional: true, emit: plots

    when:
    params.run_multiqc

    script:
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    """
    set -euo pipefail
    
    # Count input files
    FILE_COUNT=\$(ls -1 2>/dev/null | wc -l)
    echo "MultiQC found \$FILE_COUNT input files"
    
    if [ "\$FILE_COUNT" -eq 0 ]; then
        echo "ERROR: No input files for MultiQC"
        exit 1
    fi
    
    # Run MultiQC
    ${params.multiqc_bin} . \\
        --force \\
        --title "Conterminator QC Report" \\
        --comment "RNA-seq Quality Control & Contamination Detection Pipeline" \\
        --filename multiqc_report.html \\
        --dirs \\
        --dirs-depth 2 \\
        ${config} \\
        --verbose
    
    echo "MultiQC report generated successfully"
    echo "Found \$FILE_COUNT QC output files"
    """
}