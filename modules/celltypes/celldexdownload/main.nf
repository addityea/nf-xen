process CELLTYPES_CELLDEXDOWNLOAD {
    label 'lowMemST'
    tag "${ref}"

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-celldex_bioconductor-hdf5array_bioconductor-singlecellexperiment_r-yaml:c4e76f99d7b45118':
        'community.wave.seqera.io/library/bioconductor-celldex_bioconductor-hdf5array_bioconductor-singlecellexperiment_r-yaml:13bf33457e3e7490' }"

    input:
    val(ref)

    output:
    path("celldex_${ref}_h5_se.tar.gz"),        emit: tar

    script:
    template("celldexDownload.R")

}
