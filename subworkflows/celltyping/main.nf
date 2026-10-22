/*
========================================================================================
    CELLTYPING — modular cell type annotation

    Backends (params.celltyping_method):
      * singler   — SingleR + celldex references (auto-downloaded or local tarballs)
      * celltypist — CellTypist logistic-regression models (downloaded or local .pkl)
      * both
========================================================================================
*/

include { CELLTYPES_SINGLER }    from '../../modules/local/celltyping_singler'
include { CELLTYPES_CELLTYPIST } from '../../modules/local/celltyping_celltypist'
include { CELLDEX_DOWNLOAD }     from '../../modules/local/celldex_download'

workflow CELLTYPING {
    take:
    ch  // tuple( val(meta), path(h5ad) )

    main:

    // ------------------------------------------------------- SingleR
    ch_singler = ch.filter { meta, h5ad ->
        ( params.celltyping_method in ['singler', 'both'] ) &&
        meta.celldex_ref && meta.celldex_labels
    }

    // celldex_ref may hold several '--'-joined entries. Entries ending in
    // .tar.gz are used as-is; anything else is fetched via CELLDEX_DOWNLOAD.
    ch_refs = ch_singler
        .map { meta, h5ad -> meta.celldex_ref.split('--').collect { r -> r.trim() } }
        .flatten()
        .unique()

    ch_local_refs = ch_refs.filter { r -> r.endsWith('.tar.gz') }.map { r -> file(r) }
    ch_remote_refs = ch_refs.filter { r -> !r.endsWith('.tar.gz') }

    CELLDEX_DOWNLOAD(ch_remote_refs)

    ch_all_refs = ch_local_refs
        .mix( CELLDEX_DOWNLOAD.out.tar )
        .collect()

    CELLTYPES_SINGLER(ch_singler, ch_all_refs)

    // ------------------------------------------------------ CellTypist
    ch_celltypist = ch.filter { meta, h5ad ->
        ( params.celltyping_method in ['celltypist', 'both'] ) && meta.celltypist_model
    }

    CELLTYPES_CELLTYPIST(ch_celltypist, params.celltypist_majority_voting)

    emit:
    report      = CELLTYPES_SINGLER.out.report.mix( CELLTYPES_CELLTYPIST.out.report )
    summary     = CELLTYPES_SINGLER.out.summary.mix( CELLTYPES_CELLTYPIST.out.summary )
    predictions = CELLTYPES_SINGLER.out.predictions.mix( CELLTYPES_CELLTYPIST.out.predictions )
}
