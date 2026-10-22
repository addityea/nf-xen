process QC {
    tag "${meta.sample}"
    label 'lowMemMT'

    publishDir { "${params.outdir}/${meta.sample}/qc" },
        mode     : 'copy',
        overwrite: true

    container params.container_python

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${meta.sample}_qc.h5ad"),          emit: processed
    tuple val(meta), path("${meta.sample}_qc.pdf"),           emit: report
    tuple val(meta), path("${meta.sample}_qc_summary.json"),  emit: summary

    script:
    template 'qc.py'
}
