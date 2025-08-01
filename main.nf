/*
========================================================================================
    ATACseq Nextflow Workflow
========================================================================================
    GitHub   : addityea/nf-atac
    Contact  : aditya.singh@nbis.se
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

println """\
         N F - X E N    P I P E L I N E
         ===================================
         GitHub: addityea/nf-xen
         ___________________________________
         SAMPLE SHEET   : ${params.sampleSheet}
         OUTPUT DIR     : ${params.outdir}
         ___________________________________
         """
         .stripIndent()


/*
========================================================================================
    Include Sub-Workflows
========================================================================================
*/

include { RUNALL as RUNALL } from './subworkflows/runall'
include { QC_ONLY as QC_ONLY } from './subworkflows/qc_only'
include { QC_CLUST as QC_CLUST } from './subworkflows/qc_clust'
include { CLUST_ONLY as CLUST_ONLY } from './subworkflows/clust_only'
include { QC_CELLTYPE as QC_CELLTYPE } from './subworkflows/qc_celltype'
include { CLUST_CELLTYPE as CLUST_CELLTYPE } from './subworkflows/clust_celltype'
include { CELLTYPE_ONLY as CELLTYPE_ONLY } from './subworkflows/celltype_only'

if (params.sampleSheet) {
    ch_input = Channel
        .fromPath(params.sampleSheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sampleName, file(row.h5ad), row.region,row.min_counts,row.max_counts,row.min_genes,row.max_genes,row.min_cell_area,row.max_cell_area,row.min_area_ratio,row.max_area_ratio,row.max_nucleus_area,row.min_cells_per_gene,row.qc,row.clust,row.clust_res,row.celldex_ref,row.celldex_labs) }
} else {
    error "Please provide --sampleSheet"
}

ch_input_celltype_only = ch_input.filter { sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs -> qc == 'NO' && clust == 'NO' && celldex_ref != 'NA' && celldex_labs != 'NA' }
    
ch_input_qc_only = ch_input.filter { sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs -> qc == 'YES' && clust == 'NO' && celldex_ref == 'NA' && celldex_labs == 'NA' }
    
ch_input_qc_clust = ch_input.filter { sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs -> qc == 'YES' && clust == 'YES' && celldex_ref == 'NA' && celldex_labs == 'NA' }

ch_input_clust_only = ch_input.filter { sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs -> qc == 'NO' && clust == 'YES' && celldex_ref == 'NA' && celldex_labs == 'NA' }

ch_input_qc_celltype = ch_input.filter { sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs -> qc == 'YES' && clust == 'NO' && celldex_ref != 'NA' && celldex_labs != 'NA' }

ch_input_clust_celltype = ch_input.filter { sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs -> qc == 'NO' && clust == 'YES' && celldex_ref != 'NA' && celldex_labs != 'NA' }

ch_input_run_all = ch_input.filter { sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs -> qc == 'YES' && clust == 'YES' && celldex_ref != 'NA' && celldex_labs != 'NA' }

// Create a channel for each workflow step
workflow {
    RUNALL(ch_input_run_all)
    QC_ONLY(ch_input_qc_only)
    QC_CLUST(ch_input_qc_clust)
    CLUST_ONLY(ch_input_clust_only)
    QC_CELLTYPE(ch_input_qc_celltype)
    CLUST_CELLTYPE(ch_input_clust_celltype)
    CELLTYPE_ONLY(ch_input_celltype_only)
}


workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}