process COPY_SAMPLE_SHEET {
    publishDir "${params.outdir}/Input/", mode: 'copy'

    input:
    path sample_sheet

    output:
    path sample_sheet

    script:
    """
    # Just pass through - publishDir will copy it
    """
}