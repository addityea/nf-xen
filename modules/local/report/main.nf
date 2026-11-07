process REPORT {
    label 'process_single'

    publishDir "${params.outdir}",
        mode     : 'copy',
        overwrite: true

    container params.container_quarto

    input:
    path report_md
    path figures

    output:
    path "nf-xen-report.html", emit: report
    path "figures/*",          emit: figures

    script:
    """
    mkdir -p figures
    cp ${figures} figures/
    quarto render ${report_md}
    mv report.html nf-xen-report.html
    """
}
