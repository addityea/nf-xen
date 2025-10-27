#!/usr/bin/env Rscript

library(SingleR)
library(celldex)
library(yaml)
library(Seurat)
library(ggplot2)
library(anndataR)
library(HDF5Array)

# Read .h5ad file using zellkonverter
h5ad_file <- "${h5ad}" # Get the filename from environment variable
sce <- read_h5ad(h5ad_file, as = "SingleCellExperiment") # Converts .h5ad to a SingleCellExperiment object

# Split the references by -- and loop over each
references <- strsplit("${reference.join('--')}", ",")[[1]]

reference_labels <- strsplit("${label}", "--")[[1]]
stopifnot(
    #"Lengths of references and reference_labels vectors must match",
    length(references) == length(reference_labels)
)
prefix <- "${sampleName}"
Sys.setenv(XDG_CACHE_HOME = file.path(getwd(), ".cache"))
prediction_results <- list()
for (ref_idx in seq_along(references)) {
  ref <- references[ref_idx]
  reflabel <- reference_labels[ref_idx]
  ref_name <- strsplit(ref, "__")[[1]][1]
  ref_ver <- strsplit(ref, "__")[[1]][2]
  # Untar the reference files into a directory named after the reference without the extension
  ref_dir <- gsub(".tar.gz", "", ref)
  untar(ref, exdir = "./")
  # Read the SummarizedExperiment object from the provided path
  reference <- loadHDF5SummarizedExperiment(dir = ref_dir)
  stopifnot(
    reflabel %in% colnames(colData(reference))
  )
  predictions <- SingleR(
    test = assay(sce, 'counts'),
    ref = reference,
    labels = colData(reference)[[reflabel]]
  )

  # Plot and save heatmap
  p <- plotScoreHeatmap(
    predictions,
    main = paste0(
      "SingleR Predictions: ",
      basename(h5ad_file),
      " [", ref, "]"
    ),
    show_rownames = TRUE,
    show_colnames = FALSE
  )
  ggsave(
    filename = paste0(prefix, "_", ref, "_heatmap.pdf"),
    plot = p,
    width = 10,
    height = 8
  )

  # Plot and save distribution
  p2 <- plotDeltaDistribution(predictions, ncol = 3)
  p2 <- p2 + ggtitle(
    paste0(
      "SingleR Predictions: ",
      basename(h5ad_file),
      " [", ref, "]"
    )
  )
  ggsave(
    filename = paste0(prefix, "_", ref, "_distribution.pdf"),
    plot = p2,
    width = 14,
    height = 12
  )

  colnames(predictions) <- paste0(
    colnames(predictions), "_", ref_name, "_", ref_ver
  )
  prediction_results[[ref]] <- predictions
}

prediction_nrows <- lapply(prediction_results, nrow)
prediction_rownames <- lapply(prediction_results, rownames)


stopifnot(
  all(sapply(prediction_nrows, function(x) x == prediction_nrows[[1]])) &
  all(sapply(prediction_rownames, function(x) all(x == prediction_rownames[[1]])))
)

# This is predicated in the assumption that all prediction data frames have exactly
# the same rows ... see the stopifnot clause above
predictions <- do.call(cbind, prediction_results)

write.csv(
  predictions,
  file = paste0(prefix, "_predictions.csv"),
  row.names = TRUE
)
# Delete the Rplots.pdf file if it exists
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

