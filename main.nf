/*
========================================================================================
    nf-xen — end-to-end Nextflow pipeline for 10x Xenium spatial transcriptomics
----------------------------------------------------------------------------------------
    GitHub   : addityea/nf-xen
    Contact  : aditya.singh@nbis.se
    Language : Nextflow DSL2 (strict parser, Nextflow >= 25.10)
----------------------------------------------------------------------------------------
*/

include { NF_XEN } from './subworkflows/nf_xen'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf --sampleSheet <samplesheet.csv> -profile <local|docker|arm|apptainer|offline|pdc_kth|uppmax>
    Smoke test:   nextflow run main.nf -profile test,docker       (bundled toy data)
    Full test:    nextflow run main.nf -profile test_full,docker  (small public Xenium dataset)

    Required:
      --sampleSheet      Sample sheet (CSV) validated against assets/samplesheet_schema.json

    Main options:
      --outdir           Output directory (default: results)
      --clustering_method   leiden | louvain (default: leiden)
      --clustering_resolutions  Comma-separated resolutions (default: 0.5,0.8,1.0,1.2)
      --celltyping_method  singler | celltypist | both (default: singler)
      --spatial_knn      k for k-NN spatial graph (default: 6)

    For the full parameter list use: nextflow run main.nf --help
    """.stripIndent()
}

workflow {
    main:
    if (params.help) {
        helpMessage()
        exit 0
    }
    validateParameters()
    log.info paramsSummaryLog(workflow)

    if (!params.sampleSheet) {
        error "Please provide --sampleSheet. See --help for usage."
    }

    log.info """
         N F - X E N    P I P E L I N E
         ===================================
         GitHub       : addityea/nf-xen
         Nextflow     : ${nextflow.version}
         Sample sheet : ${params.sampleSheet}
         Output dir   : ${params.outdir}
         ___________________________________
         """.stripIndent()

    // Parse + validate the sample sheet. nf-schema returns one list per row:
    // [meta, input, steps, region_col, <filter thresholds...>, <annotation refs...>]
    // in the exact column order declared in assets/samplesheet_schema.json.
    def rows = samplesheetToList(params.sampleSheet, "${projectDir}/assets/samplesheet_schema.json")

    ch_samples = channel
        .fromList(rows)
        .map { meta, input, steps, region_col,
              min_counts, max_counts, min_genes, max_genes,
              min_cell_area, max_cell_area, min_area_ratio, max_area_ratio,
              max_nucleus_area, min_cells_per_gene, max_control_fraction,
              exclude_doublets, clustering_resolutions,
              celldex_ref, celldex_labels, celltypist_model ->
            def m = meta + [
                region_col          : region_col,
                steps               : steps.tokenize(',').collect { s -> s.trim() },
                filters             : [
                    min_counts        : min_counts,
                    max_counts        : max_counts,
                    min_genes         : min_genes,
                    max_genes         : max_genes,
                    min_cell_area     : min_cell_area,
                    max_cell_area     : max_cell_area,
                    min_area_ratio    : min_area_ratio,
                    max_area_ratio    : max_area_ratio,
                    max_nucleus_area  : max_nucleus_area,
                    min_cells_per_gene: min_cells_per_gene,
                    max_control_fraction: max_control_fraction,
                    exclude_doublets  : exclude_doublets,
                ],
                clustering_resolutions: clustering_resolutions ?: params.clustering_resolutions,
                celldex_ref           : celldex_ref,
                celldex_labels        : celldex_labels,
                celltypist_model      : celltypist_model ?: params.celltypist_model,
            ]
            tuple(m, file(input))
        }

    // Dump run parameters for the final report
    def params_file = file("${workflow.launchDir}/nf_xen_params.json")
    params_file.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))

    NF_XEN(ch_samples, params_file)

    onComplete:
    if (workflow.success) {
        log.info """
        Pipeline execution summary
        ---------------------------
        Completed at : ${workflow.complete}
        Duration     : ${workflow.duration}
        Success      : ${workflow.success}
        workDir      : ${workflow.workDir}
        Report       : ${params.outdir}/nf-xen-report.html
        """.stripIndent()
    } else {
        log.error "Pipeline failed: ${workflow.errorReport}"
    }
}
