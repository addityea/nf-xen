/*
========================================================================================
    NF_XEN — main workflow

    Smart entry point: every sample declares its own `steps` in the sample sheet
    (qc, clustering, celltyping, spatial). Each module runs exactly once; samples are
    routed to it by filtering the channel — no per-combination sub-workflows.
========================================================================================
*/

include { QC }          from '../../modules/local/qc'
include { CLUSTERING }  from '../../modules/local/clustering'
include { CELLTYPING }  from '../celltyping'
include { SPATIAL }     from '../../modules/local/spatial'
include { REPORT_DATA } from '../../modules/local/report_data'
include { REPORT }      from '../../modules/local/report'

workflow NF_XEN {
    take:
    ch_samples   // tuple( val(meta), path(h5ad) )
    params_file  // path: JSON dump of run parameters

    main:

    // ---------------------------------------------------------------- QC
    QC(
        ch_samples.filter { meta, h5ad -> meta.steps.contains('qc') }
    )

    // Downstream steps operate on the QC'd object when QC ran for that
    // sample, and on the raw input otherwise.
    ch_processed = QC.out.processed
        .mix( ch_samples.filter { meta, h5ad -> !meta.steps.contains('qc') } )

    // -------------------------------------------------------- Clustering
    CLUSTERING(
        ch_processed.filter { meta, h5ad -> meta.steps.contains('clustering') },
        params.clustering_method,
        params.n_neighbors,
    )

    // -------------------------------------------------------- Cell typing
    CELLTYPING(
        ch_processed.filter { meta, h5ad -> meta.steps.contains('celltyping') }
    )

    // -------------------------------------------------- Spatial analysis
    // Spatial runs on the clustered object when clustering ran for that
    // sample (cluster labels live there), else on the QC'd object.
    ch_spatial_direct = ch_processed
        .filter { meta, h5ad -> meta.steps.contains('spatial') && !meta.steps.contains('clustering') }
    ch_spatial_clustered = CLUSTERING.out.clusters
        .filter { meta, h5ad -> meta.steps.contains('spatial') }

    SPATIAL(
        ch_spatial_direct.mix(ch_spatial_clustered),
        params.spatial_knn,
        params.spatial_top_svg,
    )

    // ------------------------------------------------------------- Report
    // Every analysis process emits:
    //   report  : tuple(meta, path(*.pdf))
    //   summary : tuple(meta, path(*_summary.json))
    // Fold everything into two lists and render one Quarto report per run.
    ch_figures = QC.out.report
        .map { meta, pdfs -> pdfs }
        .flatten()
        .mix( CLUSTERING.out.report.map { meta, pdfs -> pdfs }.flatten() )
        .mix( CELLTYPING.out.report.map { meta, pdfs -> pdfs }.flatten() )
        .mix( SPATIAL.out.report.map { meta, pdfs -> pdfs }.flatten() )
        .toList()

    ch_summaries = QC.out.summary
        .map { meta, json -> json }
        .mix( CLUSTERING.out.summary.map { meta, json -> json } )
        .mix( CELLTYPING.out.summary.map { meta, json -> json } )
        .mix( SPATIAL.out.summary.map { meta, json -> json } )
        .toList()

    REPORT_DATA(ch_figures, ch_summaries, params_file)

    REPORT(
        REPORT_DATA.out.markdown,
        REPORT_DATA.out.figures,
    )

    emit:
    processed   = QC.out.processed
    clustering  = CLUSTERING.out.clusters
    predictions = CELLTYPING.out.predictions
    spatial     = SPATIAL.out.tables
    report      = REPORT.out.report
}
