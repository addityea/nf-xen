include { QC as QC } from '../../modules/qc'
include { CELLTYPE_ASSIGNMENT as CELLTYPE_ASSIGNMENT } from '../celltyping'
include { CELLTYPES_CELLDEXDOWNLOAD as CELLTYPES_CELLDEXDOWNLOAD } from '../../modules/celltypes/celldexdownload'

workflow QC_CELLTYPE {
    take:
    ch_input // channel: [ sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs ]

    main:
    QC(ch_input)
    CELLTYPES_CELLDEXDOWNLOAD(QC.out.celldex_ref)
    CELLTYPE_ASSIGNMENT(QC.out.sampleName,QC.out.h5ad, CELLTYPES_CELLDEXDOWNLOAD.out.tar, QC.out.celldex_labs)

}