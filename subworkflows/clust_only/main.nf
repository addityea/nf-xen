include { CLUSTERING as CLUSTERING } from '../../modules/clustering'

workflow CLUST_ONLY {
    take:
    ch_input // channel: [ sampleName, h5ad, region, min_counts, max_counts, min_genes, max_genes, min_cell_area, max_cell_area, min_area_ratio, max_area_ratio, max_nucleus_area, min_cells_per_gene, qc, clust, clust_res, celldex_ref, celldex_labs ]

    main:
    CLUSTERING(
        ch_input.map { row -> row.sampleName },
        ch_input.map { row -> row.h5ad },
        params.clust_method,
        ch_input.map { row -> row.clust_res == 'NA' ? params.clust_res : row.clust_res }
    )

}