process CELLTYPES_SINGLER {
    publishDir "${params.outdir}/${sampleName}/celltypes",
    mode : 'copy',
    overwrite : true
    tag "$sampleName"
    label 'highMemST'

    container 'docker.io/saditya88/singler:0.0.1'


    input:
    val(sampleName)
    path(h5ad)
    path(reference)
    val(label)

    output:
    //tuple val(meta), path("*.h5ad"), emit: h5ad
    path "*.pdf"                   , emit: pdf
    val(sampleName)
    path("*.csv")                   ,emit: obs

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'singleR.R'
}
