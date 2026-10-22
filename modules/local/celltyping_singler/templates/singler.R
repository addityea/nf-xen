#!/usr/bin/env Rscript
# nf-xen cell typing via SingleR + celldex HDF5SummarizedExperiment references.
#
# References arrive as tarballs (downloaded by CELLDEX_DOWNLOAD or provided as
# local paths in the sample sheet). celldex_ref and celldex_labels in the
# sample sheet are '--'-joined and positionally aligned.
#
# Reference choice matters: for Xenium (single-cell resolution imaging),
# reference-based annotation with SingleR is a solid baseline (Cheng et al.
# 2025, BMC Bioinformatics 26:22). Use a cell-type label granularity that
# matches the panel size (label.main for small panels).

library(SingleR)
library(yaml)
library(ggplot2)
library(anndataR)
library(HDF5Array)
library(SingleCellExperiment)

h5ad_file <- "${h5ad}"
prefix <- "${meta.sample}"
refs_str <- "${meta.celldex_ref}"
labels_str <- "${meta.celldex_labels}"

# Positional ref -> label map from the sample sheet ('--'-joined)
references <- trimws(strsplit(refs_str, "--", fixed = TRUE)[[1]])
ref_labels <- trimws(strsplit(labels_str, "--", fixed = TRUE)[[1]])
stopifnot(length(references) == length(ref_labels))

# Staged reference tarballs (downloaded by CELLDEX_DOWNLOAD or local paths).
ref_files <- trimws(strsplit("${references}", " ")[[1]])

sce <- read_h5ad(h5ad_file, as = "SingleCellExperiment")
counts_assay <- if ("counts" %in% assayNames(sce)) assay(sce, "counts") else assay(sce)

Sys.setenv(XDG_CACHE_HOME = file.path(getwd(), ".cache"))

prediction_results <- list()
summary_lines <- c("{", paste0('  "step": "singler", "sample": "', prefix, '", "references": ['))
ref_entries <- c()

for (i in seq_along(ref_files)) {
  tarball <- ref_files[i]

  # Derive the celldex ref key from the tarball name:
  #   celldex_<ref>_h5_se.tar.gz (downloaded) or <ref>_h5_se.tar.gz (local)
  tag <- strsplit(basename(tarball), ".tar.gz", fixed = TRUE)[[1]][1]
  if (endsWith(tag, "_h5_se")) tag <- substring(tag, 1, nchar(tag) - 6)
  if (startsWith(tag, "celldex_")) tag <- substring(tag, 9)

  idx <- match(tag, references)
  reflabel <- if (!is.na(idx)) ref_labels[idx] else ref_labels[1]

  exdir <- paste0("ref_", i)
  untar(tarball, exdir = exdir)
  ref_dir <- Sys.glob(file.path(exdir, "*_h5_se"))[1]
  stopifnot(!is.na(ref_dir), dir.exists(ref_dir))

  reference <- loadHDF5SummarizedExperiment(dir = ref_dir)
  stopifnot(reflabel %in% colnames(colData(reference)))

  message(sprintf("SingleR: %s (label: %s)", tag, reflabel))

  predictions <- SingleR(
    test = counts_assay,
    ref = reference,
    labels = colData(reference)[[reflabel]]
  )

  p <- plotScoreHeatmap(
    predictions,
    main = paste0("SingleR: ", prefix, " [", tag, "]"),
    show_rownames = TRUE,
    show_colnames = FALSE
  )
  ggsave(filename = paste0(prefix, "_", tag, "_heatmap.pdf"), plot = p, width = 10, height = 8)

  p2 <- plotDeltaDistribution(predictions, ncol = 3) +
    ggtitle(paste0("SingleR delta: ", prefix, " [", tag, "]"))
  ggsave(filename = paste0(prefix, "_", tag, "_distribution.pdf"), plot = p2, width = 14, height = 12)

  # Confidence summary for the report (before column renaming below)
  ref_entries <- c(ref_entries, sprintf(
    '    {"reference": "%s", "label": "%s", "n_labels": %d, "median_delta": %.4f}',
    tag, reflabel, ncol(predictions\$scores), median(as.matrix(predictions\$delta.next))
  ))

  colnames(predictions) <- paste0(colnames(predictions), "_", tag)
  prediction_results[[tag]] <- predictions
}

ref_json <- paste(ref_entries, collapse = ",\\n")
writeLines(
  c(summary_lines, ref_json, "  ]", "}"),
  paste0(prefix, "_singler_summary.json")
)

prediction_nrows <- vapply(prediction_results, nrow, integer(1L))
stopifnot(all(prediction_nrows == prediction_nrows[1]))

predictions <- do.call(cbind, prediction_results)
write.csv(predictions, file = paste0(prefix, "_predictions.csv"), row.names = TRUE)

if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
message("SingleR annotation complete")
