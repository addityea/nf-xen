include { QC as QC } from '../../modules/qc'
include { CLUSTERING as CLUSTERING } from '../../modules/clustering'

workflow QC_CLUST {
    take:
    ch_input // channel: [ sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs ]

    main:
    QC(ch_input)
    CLUSTERING(
        QC.out.sampleName,
        QC.out.h5ad,
        params.clust_method,
        QC.out.clust_res.map { it == 'NA' ? params.clust_res : it }
    )

}