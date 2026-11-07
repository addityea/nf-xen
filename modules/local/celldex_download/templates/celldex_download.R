#!/usr/bin/env Rscript
# Fetch a celldex reference (name__version, e.g. hpca__2024-02-26) and store it
# as an HDF5SummarizedExperiment tarball reusable by CELLTYPES_SINGLER.
library(celldex)
library(SingleCellExperiment)
library(HDF5Array)

r <- "${ref}"
message(paste("Fetching celldex reference:", r))

parts <- strsplit(r, "__")[[1]]
stopifnot(length(parts) == 2)
ref_name <- parts[1]
ref_version <- parts[2]

reference <- fetchReference(ref_name, ref_version, cache = "./")

dir <- paste0("celldex_", r, "_h5_se")
saveHDF5SummarizedExperiment(reference, dir = dir, replace = TRUE)
tar(tarfile = paste0(dir, ".tar.gz"), files = dir)
