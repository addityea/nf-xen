#!/usr/bin/env python3
import scanpy as sc
from pyclustree import clustree
import glob
import os
from PyPDF2 import PdfMerger
sc.set_figure_params(dpi_save=150, facecolor='white', frameon=True, vector_friendly=True)
def run_umap(adata):
    print("Running PCA...")
    sc.pp.pca(adata)
    print("Running neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=16)
    print("Running UMAP...")
    sc.tl.umap(adata)
def run_clust(sample_name, adata, resolutions, method):
    clust_keys = []
    plots = []
    for res in resolutions:
        print(f"Running clustering with resolution {res} using method {method}")
        if method == "leiden":
            sc.tl.leiden(adata, resolution=res,
                         key_added=f"leiden_{res}",
                         flavor = "igraph",
                         n_iterations=2)
            clust_keys.append(f"leiden_{res}")
            sc.pl.umap(adata, color=f"leiden_{res}",
                   title=f"Leiden clustering at resolution {res}",
                   save=f"{sample_name}_leiden_{res}.pdf")
            plots.append(f"figures/umap{sample_name}_leiden_{res}.pdf")

            sc.pl.spatial(adata, color=f"leiden_{res}",
                   title=f"Leiden clustering at resolution {res}",
                   save=f"{sample_name}_spatial_leiden_{res}.pdf",
                   spot_size = 100)
            plots.append(f"figures/show{sample_name}_spatial_leiden_{res}.pdf")
        elif method == "louvain":
            sc.tl.louvain(adata, resolution=res,
                          key_added=f"louvain_{res}",
                          flavor="vtraag")
            clust_keys.append(f"louvain_{res}")
            sc.pl.umap(adata, color=f"louvain_{res}",
                   title=f"Louvain clustering at resolution {res}",
                   save=f"{sample_name}_louvain_{res}.pdf")
            plots.append(f"figures/umap{sample_name}_louvain_{res}.pdf")
            sc.pl.spatial(adata, color=f"louvain_{res}",
                   title=f"Louvain clustering at resolution {res}",
                   save=f"{sample_name}_spatial_louvain_{res}.pdf",
                   spot_size = 100)
            plots.append(f"figures/show{sample_name}_spatial_louvain_{res}.pdf")
        else:
            raise ValueError(f"Unknown clustering method: {method}")
    if len(clust_keys) > 1:
        print("Running clustree...")
        fig = clustree(adata,clust_keys, title = f"Clustree for {sample_name}",
                    edge_weight_threshold = 0.00,
                    show_fraction = True)
        fig.set_size_inches(10, 15)
        fig.set_dpi(100)
        # Save as PDF
        fig.savefig(f"{sample_name}_clustree.pdf", bbox_inches='tight')
        plots.append(f"{sample_name}_clustree.pdf")
    return plots

# Function to concat all PDFs into one and save as sampleName_qc.pdf
def concat_pdfs(sampleName, pdf_files):
    output_pdf = f"{sampleName}_clustering.pdf"
    merger = PdfMerger()
    for pdf_file in pdf_files:
        if pdf_file != output_pdf:  # Avoid including the output file if it already exists
            merger.append(pdf_file)
    merger.write(output_pdf)
    merger.close()
    print(f"Concatenated PDFs saved to {output_pdf}")

adata = sc.read_h5ad("${h5ad}")
sampleName = "${sampleName}"
clusteringMethod = "${clusteringMethod}"
resolutions = "${resolutions}"
# Split resolutions by comma and convert to float
if isinstance(resolutions, float):
    resolutions = [resolutions]
elif isinstance(resolutions, str):
    resolutions = resolutions.split(",")
    resolutions = [float(res) for res in resolutions]
else:
    raise ValueError("Resolutions should be a float or a comma-separated string of floats.")

adata = sc.read_h5ad("${h5ad}")
run_umap(adata)
plots = run_clust(sampleName, adata, resolutions, clusteringMethod)
concat_pdfs(sampleName, plots)
# Write the updated AnnData object back to an h5ad file
adata.write(f"{sampleName}_{clusteringMethod}.h5ad")


