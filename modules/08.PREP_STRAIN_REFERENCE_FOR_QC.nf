// Reference preparation for QC tools
process PREP_STRAIN_REFERENCE_FOR_QC {
    tag "${strain}"
    publishDir "${params.outdir}/Input/strain_references/${strain}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(strain), path(fasta_gz), path(gtf)

    output:
    tuple val(strain), path("${strain}.fa"), path("${strain}.fa.fai"), path("${strain}.dict"), path("${strain}.bed"), path("${strain}_effective_size.txt"), path("${strain}.2bit"), emit: qc_references

    script:
    """
    set -euo pipefail

    echo "Preparing QC references for strain: ${strain}"

    # Decompress fasta
    zcat ${fasta_gz} > ${strain}.fa

    # Create fasta index
    ${params.samtools_bin} faidx ${strain}.fa

    # Create sequence dictionary
    ${params.samtools_bin} dict ${strain}.fa > ${strain}.dict

    # Create BED file from fai
    awk 'BEGIN{OFS="\\t"} {print \$1, 0, \$2}' ${strain}.fa.fai > ${strain}.bed

    # Calculate effective genome size
    ${params.facount_bin} ${strain}.fa > ${strain}.facount.txt
    awk 'BEGIN{len=0; n=0} \$1=="total"{len=\$2; n=\$7} END{print len-n}' ${strain}.facount.txt > ${strain}_effective_size.txt
    echo "Effective genome size for ${strain}: \$(cat ${strain}_effective_size.txt)" >&2

    # Create 2bit file if deeptools is enabled
    if [[ "${params.run_deeptools}" == "true" ]]; then
        ${params.fa2bit_bin} ${strain}.fa ${strain}.2bit
    else
        # Create empty placeholder
        touch ${strain}.2bit
    fi

    echo "QC references prepared for ${strain}"
    """
}