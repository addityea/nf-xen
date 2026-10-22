#!/usr/bin/env python3
"""nf-xen clustering step.

Leiden (default) / Louvain on the scanpy k-NN graph, PCA without HVG
selection (Xenium panels are too small for meaningful highly-variable-gene
choice), UMAP, multi-resolution scan with clustree, and spatial + UMAP
visualisations. Leiden is the recommended method; Louvain is deprecated in
scanpy and kept only for backwards compatibility.
"""

import json
import warnings

import pandas as pd
import scanpy as sc
from pypdf import PdfWriter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.set_figure_params(dpi_save=150, facecolor="white", frameon=True, vector_friendly=True)

SAMPLE = "${meta.sample}"
METHOD = "${clustering_method}"
RESOLUTIONS = [float(r) for r in "${meta.clustering_resolutions}".split(",")]
N_NEIGHBORS = int("${n_neighbors}")

summary = {"step": "clustering", "sample": SAMPLE, "method": METHOD, "resolutions": RESOLUTIONS, "n_neighbors": N_NEIGHBORS}

adata = sc.read_h5ad("${h5ad}")

# Dimensionality reduction. Xenium panels are small (tens–hundreds of genes),
# so no highly-variable-gene selection is needed.
sc.pp.pca(adata, n_comps=min(50, adata.n_vars - 1, adata.n_obs - 1))
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=min(30, adata.obsm["X_pca"].shape[1]))
sc.tl.umap(adata)

pdfs = []

def save(fig, name):
    fig.savefig(name, bbox_inches="tight")
    plt.close(fig)
    pdfs.append(name)

keys = []
for res in RESOLUTIONS:
    key = f"{METHOD}_{res}"
    print(f"Clustering: {METHOD} @ {res}")
    if METHOD == "leiden":
        sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph", n_iterations=2, directed=False)
    elif METHOD == "louvain":
        warnings.warn("Louvain is deprecated in scanpy; prefer leiden.", DeprecationWarning)
        sc.tl.louvain(adata, resolution=res, key_added=key, flavor="vtraag")
    else:
        raise ValueError(f"Unknown clustering method: {METHOD}")
    keys.append(key)
    summary.setdefault("n_clusters", {})[str(res)] = int(adata.obs[key].nunique())

    fig = sc.pl.umap(adata, color=key, title=f"{SAMPLE}: {METHOD} @ {res}", show=False, return_fig=True)
    save(fig, f"umap_{key}.pdf")
    if "spatial" in adata.obsm:
        fig = sc.pl.spatial(adata, color=key, spot_size=100, title=f"{SAMPLE}: {METHOD} @ {res}", show=False, return_fig=True)
        save(fig, f"spatial_{key}.pdf")

if len(keys) > 1:
    try:
        from pyclustree import clustree
        fig = clustree(adata, keys, title=f"{SAMPLE}: clustree", show_fraction=True)
        fig.set_size_inches(10, 12)
        save(fig, "clustree.pdf")
    except Exception as exc:  # noqa: BLE001 - clustree is optional
        print(f"WARNING: clustree failed ({exc}); skipping.")

# Use the middle resolution as the default cluster key for downstream steps
default_key = keys[len(keys) // 2]
adata.uns["nfxen_cluster_key"] = default_key

writer = PdfWriter()
for pdf in pdfs:
    writer.append(pdf)
writer.write(f"{SAMPLE}_clustering.pdf")
writer.close()

# anndata cannot write Arrow-backed string columns (pandas >= 3): normalise
for _df in (adata.obs, adata.var):
    for _c in _df.columns:
        if pd.api.types.is_string_dtype(_df[_c]) and not pd.api.types.is_object_dtype(_df[_c]):
            _df[_c] = _df[_c].astype("object")
    if not pd.api.types.is_object_dtype(_df.index):
        _df.index = _df.index.astype("object")
adata.write_h5ad(f"{SAMPLE}_{METHOD}.h5ad")
with open(f"{SAMPLE}_clustering_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print(f"Clustering done: {len(keys)} resolution(s), default key: {default_key}")
