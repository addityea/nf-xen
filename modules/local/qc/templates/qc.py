#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""nf-xen QC step for Xenium data.

Implements the QC guidance from:
  * Salmelin et al. 2025, Nature Methods — "Optimizing Xenium In Situ data
    utility by quality assessment and best-practice analysis workflows"
  * 10x Xenium Onboard Analysis output metrics (control probes / codewords,
    empty-cells, Q>=20 transcripts)
  * SPICA / Galaxy Xenium guidance: low-count floor + segmentation-area
    filters; upper count filters are usually unnecessary for Xenium because
    segmentation-driven doublets are removed via cell_area.

Conventions:
  * Any threshold set to -1 disables that filter.
  * Output AnnData keeps raw counts in layers['counts']; X is CPM(1e4)+log1p.
"""

import json
import os
import glob

import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sc.set_figure_params(dpi_save=150, facecolor="white", frameon=True, vector_friendly=True)

SAMPLE = "${meta.sample}"
H5AD = "${h5ad}"
REG_COL = "${meta.region_col}"
F = dict(
    min_counts=${meta.filters.min_counts},
    max_counts=${meta.filters.max_counts},
    min_genes=${meta.filters.min_genes},
    max_genes=${meta.filters.max_genes},
    min_cell_area=${meta.filters.min_cell_area},
    max_cell_area=${meta.filters.max_cell_area},
    min_area_ratio=${meta.filters.min_area_ratio},
    max_area_ratio=${meta.filters.max_area_ratio},
    max_nucleus_area=${meta.filters.max_nucleus_area},
    min_cells_per_gene=${meta.filters.min_cells_per_gene},
    max_control_fraction=${meta.filters.max_control_fraction},
)
EXCLUDE_DOUBLETS = "${meta.filters.exclude_doublets}".lower() == "true"

summary = {"step": "qc", "sample": SAMPLE, "filters": F, "exclude_doublets": EXCLUDE_DOUBLETS}
figures = []


def fig_to_pdf(fig, name):
    fig.savefig(name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    figures.append(name)


def off(v):
    return v is None or v == -1 or v == -1.0


adata = sc.read_h5ad(H5AD)
summary["n_cells_input"] = int(adata.n_obs)
summary["n_genes_input"] = int(adata.n_vars)

# Normalise spatial coordinates: many Xenium exports keep centroids in obs
# (x_centroid/y_centroid, X_centroid/Y_centroid, x_location/y_location, x/y)
# rather than obsm['spatial'].
if "spatial" not in adata.obsm:
    for xcol, ycol in [("x_centroid", "y_centroid"), ("X_centroid", "Y_centroid"),
                       ("x_location", "y_location"), ("x", "y")]:
        if xcol in adata.obs.columns and ycol in adata.obs.columns:
            adata.obsm["spatial"] = adata.obs[[xcol, ycol]].to_numpy(dtype=float)
            print(f"Set obsm['spatial'] from obs['{xcol}'], obs['{ycol}']")
            break

# ------------------------------------------------------------------ metrics
# Matrix is expected to hold raw counts (Xenium cell_feature_matrix /
# Xenium Analyzer exports). Verify and warn if it looks pre-processed.
X = adata.X
is_counts = np.issubdtype(X.dtype, np.integer) or np.allclose(X.data if hasattr(X, "data") else X, np.round(X.data if hasattr(X, "data") else X))
if not is_counts:
    print("WARNING: X does not look like raw counts — QC thresholds may not apply as documented.")

sc.pp.calculate_qc_metrics(adata, qc_vars=[], percent_top=None, log1p=False, inplace=True)

# Xenium-specific control metrics (present in v1 cell_feature_matrix obs)
for col, frac in [
    ("control_probe_counts", "control_probe_fraction"),
    ("control_codeword_counts", "control_codeword_fraction"),
    ("unassigned_codeword_counts", "unassigned_codeword_fraction"),
    ("genomic_control_counts", "genomic_control_fraction"),
]:
    if col in adata.obs.columns:
        adata.obs[frac] = adata.obs[col] / (adata.obs["total_counts"] + 1)
if {"control_probe_counts", "control_codeword_counts"}.issubset(adata.obs.columns):
    adata.obs["control_fraction"] = (
        adata.obs["control_probe_counts"] + adata.obs["control_codeword_counts"]
    ) / (adata.obs["total_counts"] + 1)

# Fallback: derive control fraction from negative-control genes flagged in
# var (h5ad exports) or named NegProbe*/Negative* (10x convention).
if "control_fraction" not in adata.obs.columns:
    if "negative_control" in adata.var.columns:
        _neg = adata.var["negative_control"].values.astype(bool)
    else:
        _neg = adata.var_names.astype(str).str.startswith(("NegProbe", "Negative")).values
    if _neg.any():
        _ctrl = np.asarray(adata.X[:, _neg].sum(axis=1)).ravel()
        adata.obs["control_fraction"] = _ctrl / (adata.obs["total_counts"] + 1)

# Segmentation metrics
has_area = "cell_area" in adata.obs.columns and "nucleus_area" in adata.obs.columns
if has_area and "area_ratio" not in adata.obs.columns:
    adata.obs["area_ratio"] = adata.obs["nucleus_area"] / adata.obs["cell_area"]

has_spatial = "spatial" in adata.obsm

for col in ["total_counts", "n_genes_by_counts", "cell_area", "nucleus_area", "area_ratio", "control_fraction"]:
    if col in adata.obs.columns:
        summary[f"median_{col}"] = float(np.nanmedian(adata.obs[col]))

# ------------------------------------------------------------------ plots
def violin(col, fname):
    if col not in adata.obs.columns:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(adata.obs[col].clip(upper=adata.obs[col].quantile(0.995)), bins=60, kde=True, ax=ax)
    ax.set_title(f"{SAMPLE}: {col}")
    fig_to_pdf(fig, fname)

violin("total_counts", "01_total_counts.png")
violin("n_genes_by_counts", "02_n_genes.png")
violin("cell_area", "03_cell_area.png")
violin("area_ratio", "04_area_ratio.png")

if "control_fraction" in adata.obs.columns or "unassigned_codeword_fraction" in adata.obs.columns:
    cols = [c for c in ["control_fraction", "unassigned_codeword_fraction"] if c in adata.obs.columns]
    fig, axes = plt.subplots(1, len(cols), figsize=(8 * len(cols), 6))
    if len(cols) == 1:
        axes = [axes]
    for ax, c in zip(axes, cols):
        sns.histplot(adata.obs[c], bins=60, kde=True, ax=ax)
        ax.set_title(f"{SAMPLE}: {c}")
    fig_to_pdf(fig, "05_controls.png")

# Top expressed genes — drawn manually (sc.pl.highest_expr_genes has
# scanpy/seaborn kwarg-compat issues depending on versions).
_total = np.asarray(adata.X.sum(axis=0)).ravel()
_frac = _total / _total.sum()
_top = np.argsort(-_frac)[:20]
fig, ax = plt.subplots(figsize=(8, 6))
ax.barh(range(len(_top))[::-1], _frac[_top])
ax.set_yticks(range(len(_top))[::-1], list(adata.var_names[_top]))
ax.set_xlabel("Fraction of total counts")
ax.set_title(f"{SAMPLE}: top expressed genes")
fig_to_pdf(fig, "06_highest_expr_genes.png")

if has_spatial and has_area:
    fig, axarr = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.spatial(adata, color="cell_area", spot_size=100, title=f"{SAMPLE}: cell area", ax=axarr[0], show=False)
    sc.pl.spatial(adata, color="area_ratio", spot_size=100, title=f"{SAMPLE}: nucleus/cell area ratio", ax=axarr[1], show=False)
    fig_to_pdf(fig, "07_spatial_areas.png")

# ------------------------------------------------------------------ filtering
mask = np.ones(adata.n_obs, dtype=bool)
removed = {}

def apply_mask(name, new_mask):
    global mask
    n = int((mask & ~new_mask).sum())
    removed[name] = n
    mask &= new_mask

if EXCLUDE_DOUBLETS:
    if "is_doublet" in adata.obs.columns:
        apply_mask("is_doublet", ~adata.obs["is_doublet"].astype(bool).values)
    else:
        print("WARNING: exclude_doublets requested but obs['is_doublet'] not found (Xenium Onboard Analysis >= 2.0 / Xenium Analyzer provide it). Skipping.")

if not off(F["min_counts"]):
    apply_mask("min_counts", (adata.obs["total_counts"] >= F["min_counts"]).values)
if not off(F["max_counts"]):
    apply_mask("max_counts", (adata.obs["total_counts"] <= F["max_counts"]).values)
if not off(F["min_genes"]):
    apply_mask("min_genes", (adata.obs["n_genes_by_counts"] >= F["min_genes"]).values)
if not off(F["max_genes"]):
    apply_mask("max_genes", (adata.obs["n_genes_by_counts"] <= F["max_genes"]).values)
if has_area:
    if not off(F["min_cell_area"]):
        apply_mask("min_cell_area", (adata.obs["cell_area"] >= F["min_cell_area"]).values)
    if not off(F["max_cell_area"]):
        apply_mask("max_cell_area", (adata.obs["cell_area"] <= F["max_cell_area"]).values)
    if not off(F["max_nucleus_area"]):
        apply_mask("max_nucleus_area", (adata.obs["nucleus_area"] <= F["max_nucleus_area"]).values)
    if not off(F["min_area_ratio"]):
        apply_mask("min_area_ratio", (adata.obs["area_ratio"] >= F["min_area_ratio"]).values)
    if not off(F["max_area_ratio"]):
        apply_mask("max_area_ratio", (adata.obs["area_ratio"] <= F["max_area_ratio"]).values)
if not off(F["max_control_fraction"]):
    if "control_fraction" in adata.obs.columns:
        apply_mask("max_control_fraction", (adata.obs["control_fraction"] <= F["max_control_fraction"]).values)
    else:
        print("WARNING: max_control_fraction set but control columns absent. Skipping.")

adata.obs["filtered_out"] = ~mask
if has_spatial:
    fig = sc.pl.spatial(
        adata, color="filtered_out", spot_size=100, palette=["lightgray", "red"],
        title=f"{SAMPLE}: filtered out (red)", show=False, return_fig=True,
    )
    fig_to_pdf(fig, "08_filtered_spatial.png")

adata = adata[mask].copy()
adata.obs["filtered_out"] = False
summary["cells_removed"] = removed
summary["n_cells_after_filter"] = int(adata.n_obs)

if not off(F["min_cells_per_gene"]):
    sc.pp.filter_genes(adata, min_cells=int(F["min_cells_per_gene"]), inplace=True)
summary["n_genes_after_filter"] = int(adata.n_vars)

# Post-filter distributions
for col, fname in [("total_counts", "09_post_filter_counts.png"), ("n_genes_by_counts", "10_post_filter_genes.png")]:
    violin(col, fname)

# ------------------------------------------------------------------ normalize
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Region-stratified QC (if the requested region column exists)
if REG_COL in adata.obs.columns and adata.obs[REG_COL].nunique() > 1:
    feats = [c for c in ["total_counts", "n_genes_by_counts", "cell_area", "area_ratio"] if c in adata.obs.columns]
    fig, axes = plt.subplots(1, len(feats), figsize=(6 * len(feats), 5), squeeze=False)
    for ax, feat in zip(axes[0], feats):
        sns.boxplot(x=adata.obs[REG_COL], y=adata.obs[feat], ax=ax)
        ax.set_title(f"{SAMPLE}: {feat} by {REG_COL}")
        ax.tick_params(axis="x", rotation=45)
    fig_to_pdf(fig, "11_region_qc.png")
else:
    if REG_COL != "all_regions":
        print(f"WARNING: region column '{REG_COL}' not found in obs; skipping region QC plots.")

# ------------------------------------------------------------------ outputs
imgs = [Image.open(p).convert("RGB") for p in sorted(figures)]
if imgs:
    imgs[0].save(f"{SAMPLE}_qc.pdf", "PDF", resolution=150,
                  save_all=True, append_images=imgs[1:])

# anndata cannot write Arrow-backed string columns (pandas >= 3): normalise
for _df in (adata.obs, adata.var):
    for _c in _df.columns:
        if pd.api.types.is_string_dtype(_df[_c]) and not pd.api.types.is_object_dtype(_df[_c]):
            _df[_c] = _df[_c].astype("object")
    if not pd.api.types.is_object_dtype(_df.index):
        _df.index = _df.index.astype("object")
adata.write_h5ad(f"{SAMPLE}_qc.h5ad")

summary["n_cells_output"] = int(adata.n_obs)
summary["n_genes_output"] = int(adata.n_vars)
with open(f"{SAMPLE}_qc_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)

print(f"QC done: {summary['n_cells_input']} -> {summary['n_cells_output']} cells, "
      f"{summary['n_genes_input']} -> {summary['n_genes_output']} genes")
