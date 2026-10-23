#!/usr/bin/env python3
"""nf-xen spatial analysis step (squidpy).

Squidpy (Palla et al. 2022, Nat. Methods; 1.8.x actively maintained) provides
the standard Xenium spatial toolkit — see the official squidpy Xenium tutorial:
  * spatial k-NN graph on cell centroids
  * neighbourhood enrichment of clusters (sq.gr.nhood_enrichment)
  * cluster co-occurrence (sq.gr.co_occurrence)
  * spatially variable genes via Moran's I (sq.gr.spatial_autocorr)
"""

import json

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
import itertools

_pages = itertools.count()

class _PngPages:
    """Collect figures as PNG pages, merge into one compact PDF at close."""
    def __init__(self, out):
        self.out = out

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        imgs = [Image.open(p).convert("RGB") for p in sorted(pdfs)]
        if imgs:
            imgs[0].save(self.out, "PDF", resolution=150,
                         save_all=True, append_images=imgs[1:])
        return False

    def savefig(self, fig):
        name = f"fig_{next(_pages):02d}.png"
        fig.savefig(name, dpi=150, bbox_inches="tight")
        pdfs.append(name)


SAMPLE = "${meta.sample}"
KNN = int("${spatial_knn}")
TOP_SVG = int("${top_svg}")

adata = sc.read_h5ad("${h5ad}")
assert "spatial" in adata.obsm, "obsm['spatial'] missing — cannot run spatial analysis"

summary = {"step": "spatial", "sample": SAMPLE, "knn": KNN}
pdfs = []

# ---------------------------------------------------------------- graph
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=KNN, copy=False)
summary["n_cells"] = int(adata.n_obs)

# Cluster key: prefer an explicit one, else the nf-xen clustering default,
# else any leiden*/louvain*/cluster* column with >=2 categories.
cluster_key = adata.uns.get("nfxen_cluster_key")
if cluster_key not in adata.obs.columns:
    cands = [c for c in adata.obs.columns if c.startswith(("leiden", "louvain", "cluster"))]
    cluster_key = next((c for c in cands if adata.obs[c].nunique() >= 2), None)
if cluster_key is not None and adata.obs[cluster_key].nunique() < 2:
    print(f"WARNING: cluster key '{cluster_key}' has <2 categories; skipping enrichment plots.")
    cluster_key = None

with _PngPages(f"{SAMPLE}_spatial.pdf") as pdf:
    # ------------------------------------------------ neighbourhood enrichment
    if cluster_key:
        summary["cluster_key"] = str(cluster_key)
        sq.gr.nhood_enrichment(adata, cluster_key=cluster_key, copy=False)
        en = adata.uns[f"{cluster_key}_nhood_enrichment"]["zscore"]
        summary["max_positive_zscore"] = float(np.max(en))

        fig, ax = plt.subplots(figsize=(max(8, 0.6 * en.shape[0]), max(6, 0.5 * en.shape[0])))
        sq.pl.nhood_enrichment(adata, cluster_key=cluster_key, title=f"{SAMPLE}: neighbourhood enrichment", ax=ax)
        pdf.savefig(fig)
        plt.close(fig)

        sq.gr.co_occurrence(adata, cluster_key=cluster_key, copy=False)
        plt.figure(figsize=(max(8, 0.6 * en.shape[0]), max(6, 0.5 * en.shape[0])))
        sq.pl.co_occurrence(adata, cluster_key=cluster_key)
        plt.title(f"{SAMPLE}: co-occurrence ({cluster_key})")
        pdf.savefig(plt.gcf())
        plt.close("all")
    else:
        print("WARNING: no cluster/region key found; skipping neighbourhood enrichment.")

    # ------------------------------------------------------ spatially variable genes
    # Moran's I permutation test on normalised expression
    sq.gr.spatial_autocorr(adata, mode="moran", n_perms=100, corr_method="fdr_bh", copy=False, show_progress_bar=False)
    morans = adata.uns["moranI"]
    morans.to_csv(f"{SAMPLE}_morans_i.csv")
    summary["n_genes_tested"] = int(morans.shape[0])
    sig = morans[morans["pval_norm_fdr_bh"] < 0.05]
    summary["n_svg_fdr05"] = int(sig.shape[0])

    top = morans.sort_values("I", ascending=False).head(TOP_SVG)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.25 * len(top))))
    top["I"].sort_values().plot.barh(ax=ax)
    ax.set_title(f"{SAMPLE}: top {TOP_SVG} Moran's I genes")
    pdf.savefig(fig)
    plt.close(fig)

    genes = list(top.index[: min(9, TOP_SVG)])
    if genes:
        fig = sc.pl.spatial(
            adata, color=genes, spot_size=100, ncols=3,
            title=[f"{SAMPLE}: {g}" for g in genes], show=False, return_fig=True,
        )
        pdf.savefig(fig)
        plt.close(fig)

    # Cluster map for context
    if cluster_key:
        fig = sc.pl.spatial(adata, color=cluster_key, spot_size=100,
                            title=f"{SAMPLE}: {cluster_key}", show=False, return_fig=True)
        pdf.savefig(fig)
        plt.close(fig)

with open(f"{SAMPLE}_spatial_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print(f"Spatial analysis done: {summary.get('n_svg_fdr05', 'NA')} SVGs at FDR<0.05")
