process SPATIAL {
    tag "${meta.sample}"
    label 'highMemMT'

    publishDir { "${params.outdir}/${meta.sample}/spatial" },
        mode     : 'copy',
        overwrite: true

    container params.container_spatial

    input:
    tuple val(meta), path(h5ad)
    val(spatial_knn)
    val(top_svg)

    output:
    tuple val(meta), path("${meta.sample}_morans_i.csv"),        emit: tables
    tuple val(meta), path("${meta.sample}_spatial.pdf"),         emit: report
    tuple val(meta), path("${meta.sample}_spatial_summary.json"), emit: summary

    script:
    template 'spatial.py'
}
