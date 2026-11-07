process CLUSTERING {
    tag "${meta.sample}"
    label 'highMemMT'

    publishDir { "${params.outdir}/${meta.sample}/clustering" },
        mode     : 'copy',
        overwrite: true

    container params.container_clustering

    input:
    tuple val(meta), path(h5ad)
    val(clustering_method)
    val(n_neighbors)

    output:
    tuple val(meta), path("${meta.sample}_${clustering_method}.h5ad"),        emit: clusters
    tuple val(meta), path("${meta.sample}_clustering.pdf"),                   emit: report
    tuple val(meta), path("${meta.sample}_clustering_summary.json"),          emit: summary

    script:
    template 'clustering.py'
}
