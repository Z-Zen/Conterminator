process CONTAMINATION_FINAL_REPORT_MULTI {
    tag "${sample}"
    publishDir "${params.outdir}/Output/contamination_check/${sample}/final_report", mode: 'copy'
    cpus 2
    memory '4 GB'
    errorStrategy 'ignore'

    input:
    tuple val(sample), path(mapped_stats), path(unmapped_stats), path(blast_summaries), path(top_hits_files)

    output:
    path "${sample}.contamination_final_summary.txt", emit: summary

    when:
    params.run_contamination_check

    script:
    """
    set -euo pipefail
    
    cat > ${sample}.contamination_final_summary.txt << 'EOF'
COMPREHENSIVE CONTAMINATION REPORT: ${sample}
========================================
EOF
    
    echo "" >> ${sample}.contamination_final_summary.txt
    echo "MAPPED READS ANALYSIS:" >> ${sample}.contamination_final_summary.txt
    cat ${mapped_stats} >> ${sample}.contamination_final_summary.txt
    
    echo "" >> ${sample}.contamination_final_summary.txt
    echo "UNMAPPED READS ANALYSIS:" >> ${sample}.contamination_final_summary.txt
    cat ${unmapped_stats} >> ${sample}.contamination_final_summary.txt
    
    echo "" >> ${sample}.contamination_final_summary.txt
    echo "BLAST CONTAMINATION SCREENING (ALL DATABASES):" >> ${sample}.contamination_final_summary.txt
    echo "===========================================" >> ${sample}.contamination_final_summary.txt
    
    for summary in ${blast_summaries}; do
        echo "" >> ${sample}.contamination_final_summary.txt
        cat \$summary >> ${sample}.contamination_final_summary.txt
        echo "" >> ${sample}.contamination_final_summary.txt
    done
    """
}