#!/usr/bin/env python3
"""nf-xen cell typing via CellTypist (Teichlab).

CellTypist (Zhang et al. 2022, Nat. Methods) annotates cells with pre-trained
logistic-regression models; actively maintained (1.7.1, 2025). Models are
downloaded on demand, or local .pkl paths can be given in the sample sheet
(offline runs). Works best on log1p-normalised data — the QC step provides
that; if celltyping runs without QC we normalise in-memory.
"""

import json
import os

import numpy as np
import pandas as pd
import scanpy as sc
os.environ.setdefault("CELLTYPIST_FOLDER", "./celltypist_models")
import celltypist
from celltypist import models

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
import itertools

pdfs = []

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
MODELS = [m.strip() for m in "${meta.celltypist_model}".split("--") if m.strip()]
MAJORITY_VOTING = ${majority_voting ? 'True' : 'False'}

adata = sc.read_h5ad("${h5ad}")

# CellTypist expects log-normalised expression
X = adata.X
maxv = X.max() if hasattr(X, "max") else np.asarray(X).max()
looks_counts = np.issubdtype(X.dtype, np.integer) or maxv > 100
if looks_counts:
    print("Input looks like raw counts — normalising in-memory for CellTypist.")
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

# celltypist 1.7.x annotate() takes a single model; loop over models and
# combine prediction columns.
pred_frames = {}
for m in MODELS:
    if not os.path.exists(m):
        models.download_models(model=m)
    model_obj = models.Model.load(m)
    predictions = celltypist.annotate(
        adata,
        model=model_obj,
        majority_voting=MAJORITY_VOTING,
    )
    po = predictions.to_adata()
    cols = [c for c in ["predicted_labels", "majority_voting", "conf_majority_voting"] if c in po.obs.columns]
    df = po.obs[cols].copy()
    if len(MODELS) > 1:
        tag = os.path.splitext(os.path.basename(m))[0]
        df.columns = [f"{c}_{tag}" for c in cols]
    pred_frames[m] = df

pred_df = pd.concat(pred_frames.values(), axis=1)
pred_df.to_csv(f"{SAMPLE}_celltypist_predictions.csv", index=True)

out = adata.copy()
for c in pred_df.columns:
    out.obs[c] = pred_df[c].values

def _pick(kind, exclude=()):
    return next((c for c in out.obs.columns if kind in c and not any(e in c for e in exclude)), None)

LABEL_COL = _pick("majority_voting", exclude=("conf",)) or _pick("predicted_labels")
CONF_COL = _pick("conf_majority_voting") or _pick("conf_score")

# ---- figures
with _PngPages(f"{SAMPLE}_celltypist.pdf") as pdf:
    vc = out.obs[LABEL_COL].value_counts()
    fig, ax = plt.subplots(figsize=(max(8, 0.4 * len(vc)), 6))
    vc.plot.barh(ax=ax)
    ax.set_title(f"{SAMPLE}: CellTypist labels ({', '.join(MODELS)})")
    ax.set_xlabel("n cells")
    pdf.savefig(fig)
    plt.close(fig)

    if CONF_COL:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(out.obs[CONF_COL].astype(float), bins=50)
        ax.set_title(f"{SAMPLE}: majority-vote confidence")
        pdf.savefig(fig)
        plt.close(fig)

    if "spatial" in out.obsm:
        out.obs[LABEL_COL] = out.obs[LABEL_COL].astype(str)
        fig = sc.pl.spatial(
            out, color=[LABEL_COL], spot_size=100,
            title=[f"{SAMPLE}: {LABEL_COL}"], show=False, return_fig=True,
        )
        pdf.savefig(fig)
        plt.close(fig)

summary = {
    "step": "celltypist",
    "sample": SAMPLE,
    "models": MODELS,
    "majority_voting": bool(MAJORITY_VOTING),
    "n_labels": int(vc.shape[0]),
    "top_labels": {str(k): int(v) for k, v in vc.head(10).items()},
    "median_confidence": float(pd.to_numeric(out.obs[CONF_COL], errors="coerce").median()) if CONF_COL else None,
}
with open(f"{SAMPLE}_celltypist_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print(f"CellTypist done: {vc.shape[0]} labels")
