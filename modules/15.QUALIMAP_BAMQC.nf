process QUALIMAP_BAMQC {
    tag "${sample}_${bam_type}"
    publishDir "${params.outdir}/Output/qualimap/${sample}/bamqc", mode: 'copy'
    cpus { Math.min(params.qualimap_threads, params.max_cpus as int) }
    memory { params.max_mem }
    errorStrategy 'retry'; maxRetries 2
    maxForks params.max_parallel_samples

    input:
    tuple val(sample), val(strain), path(bam), path(bai), val(bam_type)
    path ref

    output:
    path "qualimap_bamqc/**", emit: reports
    path "qualimap_bamqc/qualimapReport.html", optional: true, emit: html
    path "qualimap_bamqc/qualimapReport.pdf", optional: true, emit: pdf
    path "qualimap_bamqc/genome_results.txt", optional: true, emit: results

    when:
    params.run_qualimap && (params.qualimap_mode == "bamqc" || params.qualimap_mode == "both")

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
    
    ${params.qualimap_bin} bamqc \\
        -bam ${bam} \\
        -c \\
        -gd ${params.qualimap_genome} \\
        -gff \${GTF_FILE} \\
        -outdir qualimap_bamqc \\
        -outformat PDF:HTML \\
        -p ${params.qualimap_protocol} \\
        -nt ${task.cpus}
    """
}