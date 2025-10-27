include { CELLTYPE_ASSIGNMENT as CELLTYPE_ASSIGNMENT } from '../celltyping'
include { CELLTYPES_CELLDEXDOWNLOAD as CELLTYPES_CELLDEXDOWNLOAD } from '../../modules/celltypes/celldexdownload'

workflow CELLTYPE_ONLY {
    take:
    ch_input // channel: [ sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs ]

    main:
    CELLTYPES_CELLDEXDOWNLOAD(ch_input.map { row -> row[16] }) // celldex_ref
    CELLTYPE_ASSIGNMENT(
        ch_input.map { row -> row[0] }, // sampleName
        ch_input.map { row -> row[1] }, // h5ad
        CELLTYPES_CELLDEXDOWNLOAD.out.tar,
        ch_input.map { row -> row[17] } // celldex_labs
    )

}