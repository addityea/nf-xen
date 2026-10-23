process CELLDEX_DOWNLOAD {
    tag "${ref}"
    label 'lowMemST'

    container params.container_celldex

    input:
    val(ref)

    output:
    path("celldex_${ref}_h5_se.tar.gz"), emit: tar

    script:
    template 'celldex_download.R'
}
