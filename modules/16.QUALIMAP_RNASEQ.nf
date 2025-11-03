process QUALIMAP_RNASEQ {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/qualimap/${sample}/rnaseq", mode: 'copy'
    cpus { Math.min(params.qualimap_threads, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'retry'; maxRetries 2
    maxForks params.max_parallel_samples

    input:
    tuple val(sample), val(strain), path(bam), path(bai), val(bam_type)

    output:
    path "qualimap_rnaseq_unique_mapped_reads/**", emit: reports
    path "qualimap_rnaseq_unique_mapped_reads/qualimapReport.html", optional: true, emit: html
    path "qualimap_rnaseq_unique_mapped_reads/qualimapReport.pdf", optional: true, emit: pdf
    path "qualimap_rnaseq_unique_mapped_reads/rnaseq_results.txt", optional: true, emit: results
    path "qualimap_rnaseq_proportional/**", emit: reports_prop
    path "qualimap_rnaseq_proportional/qualimapReport.html", optional: true, emit: html_prop
    path "qualimap_rnaseq_proportional/qualimapReport.pdf", optional: true, emit: pdf_prop
    path "qualimap_rnaseq_proportional/rnaseq_results.txt", optional: true, emit: results_prop

    when:
    params.run_qualimap && (params.qualimap_mode == "rnaseq" || params.qualimap_mode == "both")

    script:
    def gtf = "/home/abadreddine/rcp_storage/common/Users/vonalven/HDP_pseudogenomes_construction/Data/HPC_results/HDP_pseudogenomes/${strain}/HDP_merge_splitnorm_v1__pseudogenome__strain_${strain}.gtf.gz"
    """
    set -euo pipefail

    if [[ ${gtf} == *.gz ]]; then
        zcat ${gtf} > ${strain}_annotation.gtf
        GTF_FILE="${strain}_annotation.gtf"
    else
        GTF_FILE="${gtf}"
    fi
    
    ${params.qualimap_bin} rnaseq \\
        --algorithm uniquely-mapped-reads \\
        -bam ${bam} \\
        -gtf \${GTF_FILE} \\
        -outdir qualimap_rnaseq_unique_mapped_reads \\
        -outformat PDF:HTML \\
        -p ${params.qualimap_protocol} \\
        -pe \\
        -s

    ${params.qualimap_bin} rnaseq \\
        --algorithm proportional \\
        -bam ${bam} \\
        -gtf \${GTF_FILE} \\
        -outdir qualimap_rnaseq_proportional \\
        -outformat PDF:HTML \\
        -p ${params.qualimap_protocol} \\
        -pe \\
        -s 
    """
}