process PREPARE_STRAIN_REFERENCE {
    tag "${strain}"
    publishDir "${params.star_index_dir}", mode: 'copy', pattern: "${strain}"
    
    input:
    val strain

    output:
    tuple val(strain), path("${strain}_genome.fa.gz"), path("${strain}_annotation.gtf.gz"), emit: references

    script:
    """
    set -euo pipefail
    
    STRAIN_DIR="${params.strains_base_dir}/${strain}"
    
    FASTA=\$(find "\${STRAIN_DIR}" -name "*pseudogenome__strain_${strain}.fa.gz" | head -1)
    if [ -z "\${FASTA}" ]; then
        echo "ERROR: Could not find pseudogenome fasta for ${strain}" >&2
        exit 1
    fi
    
    GTF=\$(find "\${STRAIN_DIR}" -name "*pseudogenome__strain_${strain}.gtf.gz" | head -1)
    if [ -z "\${GTF}" ]; then
        echo "ERROR: Could not find GTF for ${strain}" >&2
        exit 1
    fi
    
    ln -s "\${FASTA}" ${strain}_genome.fa.gz
    ln -s "\${GTF}" ${strain}_annotation.gtf.gz
    
    echo "Prepared references for ${strain}"
    echo "FASTA: \${FASTA}"
    echo "GTF: \${GTF}"
    """
}