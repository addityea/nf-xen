#!/usr/bin/env python3
"""Generate tests/toy_xenium.h5ad — a tiny synthetic Xenium-like AnnData used
by the `test` profile. Deterministic (fixed seed), ~400 cells x ~80 genes,
raw counts in X, Xenium-style obs columns and spatial coordinates.

Run with the h5ad-inspect venv or any env with anndata/numpy/pandas:
    python3 bin/make_toy_h5ad.py
"""

import numpy as np
import pandas as pd
import anndata as ad

SEED = 42
rng = np.random.default_rng(SEED)

N_CELLS = 400
CELL_TYPES = {
    "T_cell": ["CD3D", "CD3E", "CD3G", "CD2", "IL7R", "CD27", "CD28", "TRAC",
               "TRBC2", "LCK", "ZAP70", "CD247", "THEMIS", "SKAP1"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "HLA-DRA", "HLA-DRB1",
               "PAX5", "EBF1", "IGLL1", "CD22", "TNFRSF13C", "BLK"],
    "Monocyte": ["CD14", "FCGR3A", "LYZ", "S100A8", "S100A9", "CSF3R",
                 "FCN1", "VCAN", "CD68", "ITGAM", "MNDA", "CEBPE"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1", "GZMB", "GZMK",
                "NCAM1", "FCGR3B", "IL2RB"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "CDH1", "CLDN4", "KRT19",
                   "TJP1", "DSP", "MAL", "GATA6", "SOX9", "CFTR"],
}
# Shared "housekeeping-ish" background genes across all cell types
BACKGROUND = [
    "ACTB", "GAPDH", "B2M", "TMSB4X", "FTL", "MT1G", "MT2A", "MALAT1",
    "EEF1A1", "RPL13A", "TPT1", "SDS", "FTH1", "PABPC1", "ODC1", "SERPINA1",
]
NEG_PROBES = ["NegProbeWF-BK001", "NegProbeWF-BK002", "NegProbeWB-BK003"]

genes = []
for g in CELL_TYPES.values():
    genes.extend(g)
genes += BACKGROUND + NEG_PROBES
genes = list(dict.fromkeys(genes))

# cell-type assignment with spatial clustering: each type occupies a blob
centers = {ct: (rng.uniform(0, 2000), rng.uniform(0, 2000)) for ct in CELL_TYPES}
labels = rng.choice(list(CELL_TYPES), size=N_CELLS)
coords = np.array([
    (lambda c: (c[0] + rng.normal(0, 120), c[1] + rng.normal(0, 120)))(centers[ct])
    for ct in labels
])

rows = []
for i, ct in enumerate(labels):
    marker_genes = CELL_TYPES[ct]
    for g in genes:
        if g in NEG_PROBES:
            n = rng.poisson(0.15)
        elif g in marker_genes:
            n = rng.poisson(8)
        elif g in BACKGROUND:
            n = rng.poisson(2)
        else:
            n = rng.poisson(0.05)
        rows.append(n)
X = np.array(rows, dtype=np.int64).reshape(N_CELLS, len(genes))

# inject some QC-relevant outliers
X[0] = 0                       # empty cell
X[1, :5] = 1                   # very low complexity
X[2] = rng.poisson(60, len(genes))  # doublet-ish high counts

cell_area = rng.normal(120, 40, N_CELLS).clip(10, 400)
nucleus_area = cell_area * rng.uniform(0.2, 0.8, N_CELLS)
cell_area[3] = 5               # too small
nucleus_area[4] = cell_area[4] * 0.99  # area ratio outlier

obs = pd.DataFrame(
    {
        "cell_area": cell_area,
        "nucleus_area": nucleus_area,
        "area_ratio": nucleus_area / cell_area,
        "is_doublet": np.array([i == 2 for i in range(N_CELLS)]),
        "predicted_celltype": labels,
    },
    index=pd.Index([f"cell_{i}" for i in range(N_CELLS)], name="cell_id"),
)

var = pd.DataFrame(index=pd.Index(genes, name="gene_id"))
var["negative_control"] = var.index.str.startswith("NegProbe")

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.layers["counts"] = X.copy()
adata.obsm["spatial"] = coords.astype(np.float64)
adata.uns["sample"] = "toy"

adata.write_h5ad("tests/toy_xenium.h5ad")
print(f"wrote tests/toy_xenium.h5ad: {adata.n_obs} cells x {adata.n_vars} genes")
