#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
library(celldex)
library(SingleCellExperiment)
library(yaml)
library(HDF5Array)
r <- "${ref}"
print(paste("Attempting to fetch reference:", r))
# Split the reference into ref_name and ref_version based on __
ref_name <- strsplit(r, "__")[[1]][1]
ref_version <- strsplit(r, "__")[[1]][2]
reference <- fetchReference(ref_name, ref_version, cache = "./")
# Save SummarizedExperiment to HDF5 files
saveHDF5SummarizedExperiment(
  reference,
  dir = paste0("celldex_", r, "_h5_se"),
  replace = TRUE
)
# Compress the HDF5 files into a tar.gz archive
tar(tarfile = paste0("celldex_", r, "_h5_se.tar.gz"), files = paste0("celldex_", r, "_h5_se"))

