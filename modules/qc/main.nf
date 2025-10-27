process QC {
    publishDir "${params.outdir}/${sampleName}/qc",
    mode : 'copy',
    overwrite : true

    tag "${sampleName}"
    label "lowMemST"

    container "docker.io/saditya88/nf-xen:qc"

    input:
    tuple val(sampleName), path(h5ad), val(reg_col),val(min_counts), val(max_counts),
    val(min_genes), val(max_genes), val(min_cell_area), val(max_cell_area),
    val(min_area_ratio), val(max_area_ratio), val(max_nucleus_area), val(min_cells_per_gene),
    val(qc), val(clust), val(clust_res), val(celldex_ref), val(celldex_labs)

    output:
    path("${sampleName}_qc.h5ad"), emit: "h5ad"
    path "${sampleName}_qc.pdf", emit: "qc_pdf"
    val sampleName, emit: "sampleName"
    val clust_res, emit: "clust_res"
    val celldex_ref, emit: "celldex_ref"
    val celldex_labs, emit: "celldex_labs"

    script:
    template "qc.py"
}