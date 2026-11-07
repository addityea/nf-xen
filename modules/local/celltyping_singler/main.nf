process CELLTYPES_SINGLER {
    tag "${meta.sample}"
    label 'highMemST'

    publishDir { "${params.outdir}/${meta.sample}/celltypes" },
        mode     : 'copy',
        overwrite: true

container params.container_singler

    input:
    tuple val(meta), path(h5ad)
    path(references)

    output:
    tuple val(meta), path("${meta.sample}_predictions.csv"),      emit: predictions
    tuple val(meta), path("${meta.sample}_*.png"),                emit: report
    tuple val(meta), path("${meta.sample}_singler_summary.json"), emit: summary

    script:
    template 'singler.R'
}
