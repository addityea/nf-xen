process CLUSTERING {
    publishDir "${params.outdir}/${sampleName}/clustering",
    mode : 'copy',
    overwrite : true

    tag "${sampleName}"
    label "lowMemMT"

    container "docker.io/saditya88/nf-xen:clust"

    input:
    val(sampleName)
    path(h5ad)
    val(clusteringMethod)
    val(resolutions)

    output:
    path "${sampleName}_${clusteringMethod}.h5ad", emit: "clustering_h5ad"
    path "${sampleName}_clustering.pdf", emit: "clustering_pdf"

    script:
    template "clustering.py"
}