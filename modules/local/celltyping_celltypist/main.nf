process CELLTYPES_CELLTYPIST {
    tag "${meta.sample}"
    label 'highMemMT'

    publishDir { "${params.outdir}/${meta.sample}/celltypes" },
        mode     : 'copy',
        overwrite: true

    container params.container_celltypist

    input:
    tuple val(meta), path(h5ad)
    val(majority_voting)

    output:
    tuple val(meta), path("${meta.sample}_celltypist_predictions.csv"),   emit: predictions
    tuple val(meta), path("${meta.sample}_celltypist.pdf"),               emit: report
    tuple val(meta), path("${meta.sample}_celltypist_summary.json"),      emit: summary

    script:
    template 'celltypist.py'
}
