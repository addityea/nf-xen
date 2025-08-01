#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import scanpy as sc
from scipy.sparse import csr_matrix
import glob
from PyPDF2 import PdfMerger

sc.set_figure_params(dpi_save=150, facecolor='white', frameon=True, vector_friendly=True)

def plot_counts_per_cell(adata_list, sample_names, pdf_filename="01_counts_per_cell.pdf"):
    all_counts = []
    all_samples = []
    for adata, name in zip(adata_list, sample_names):
        if 'total_counts' not in adata.obs.columns:
            if isinstance(adata.X, csr_matrix):
                total_counts = np.array(adata.X.sum(axis=1)).flatten()
            else:
                total_counts = adata.X.sum(axis=1)
            adata.obs['total_counts'] = total_counts
        all_counts.append(adata.obs['total_counts'])
        all_samples.extend([name] * adata.shape[0])
    all_counts = np.concatenate(all_counts)
    df = pd.DataFrame({'total_counts': all_counts, 'sample': all_samples})
    with PdfPages(pdf_filename) as pdf:
        plt.figure(figsize=(8, 6))
        sns.violinplot(x='sample', y='total_counts', hue='sample', data=df, legend=False)
        plt.ylabel("Total Transcript Counts per Cell")
        plt.title("Total Transcript Counts per Cell by Sample")
        plt.tight_layout()
        pdf.savefig()
        plt.close()
    print(f"Saved counts per cell plot to {pdf_filename}")

def plot_ngenes_per_cell(adata_list, sample_names, pdf_filename="02_ngenes_per_cell.pdf"):
    all_n_genes = []
    all_samples = []
    for adata, name in zip(adata_list, sample_names):
        if 'n_genes' not in adata.obs.columns:
            if isinstance(adata.X, csr_matrix):
                n_genes = np.array((adata.X > 0).sum(axis=1)).flatten()
            else:
                n_genes = (adata.X > 0).sum(axis=1)
            adata.obs['n_genes'] = n_genes
        all_n_genes.append(adata.obs['n_genes'])
        all_samples.extend([name] * adata.shape[0])
    all_n_genes = np.concatenate(all_n_genes)
    df = pd.DataFrame({'n_genes': all_n_genes, 'sample': all_samples})
    with PdfPages(pdf_filename) as pdf:
        plt.figure(figsize=(8, 6))
        sns.violinplot(x='sample', y='n_genes', hue='sample', legend=False, data=df)
        plt.xlabel("Sample")
        plt.ylabel("Number of Genes Detected per Cell")
        plt.title("Number of Genes Detected per Cell by Sample")
        plt.tight_layout()
        pdf.savefig()
        plt.close()
    print(f"Saved n_genes per cell plot to {pdf_filename}")

def plot_highest_expr_genes(adata_list, sample_names, pdf_filename="03_highest_expr_genes.pdf"):
    n = len(adata_list)
    with PdfPages(pdf_filename) as pdf:
        if n == 1:
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.highest_expr_genes(adata_list[0], n_top=20, ax=ax, show=False)
            ax.set_title(f"Sample {sample_names[0]}")
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
        else:
            nrows = int(np.ceil(n / 3))
            fig, axes = plt.subplots(nrows, 3, figsize=(8, 6 * nrows))
            axes = axes.flatten()
            for i, adata in enumerate(adata_list):
                sc.pl.highest_expr_genes(adata, n_top=20, ax=axes[i], show=False)
                axes[i].set_title(f"Sample {sample_names[i]}")
            for j in range(i+1, nrows*3):
                fig.delaxes(axes[j])
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Saved highest expressed genes plots to {pdf_filename}")

def plot_area_distributions(adata_list, sample_names, pdf_filename="04_area_distributions.pdf"):
    all_cell_areas = []
    all_nucleus_areas = []
    all_samples = []
    for adata, name in zip(adata_list, sample_names):
        all_cell_areas.append(adata.obs['cell_area'])
        all_nucleus_areas.append(adata.obs['nucleus_area'])
        all_samples.extend([name] * adata.shape[0])
    all_cell_areas = np.concatenate(all_cell_areas)
    all_nucleus_areas = np.concatenate(all_nucleus_areas)
    df_cell = pd.DataFrame({'cell_area': all_cell_areas, 'sample': all_samples})
    df_nucleus = pd.DataFrame({'nucleus_area': all_nucleus_areas, 'sample': all_samples})
    with PdfPages(pdf_filename) as pdf:
        plt.figure(figsize=(8, 6))
        sns.violinplot(x='sample', y='cell_area', hue='sample', legend=False, data=df_cell)
        plt.title("Cell Area Distribution by Sample")
        plt.xlabel("Sample")
        plt.ylabel("Cell Area")
        plt.tight_layout()
        pdf.savefig()
        plt.close()
        plt.figure(figsize=(8, 6))
        sns.violinplot(x='sample', y='nucleus_area', hue='sample', legend=False, data=df_nucleus)
        plt.title("Nucleus Area Distribution by Sample")
        plt.xlabel("Sample")
        plt.ylabel("Nucleus Area")
        plt.tight_layout()
        pdf.savefig()
        plt.close()
    print(f"Saved area distributions plots to {pdf_filename}")

def plot_cell_vs_nucleus_area(adata_list, sample_names, pdf_filename="05_cell_vs_nucleus_area.pdf"):
    n = len(adata_list)
    with PdfPages(pdf_filename) as pdf:
        if n == 1:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(adata_list[0].obs['cell_area'], adata_list[0].obs['nucleus_area'], alpha=0.5)
            ax.set_xlabel("Cell Area")
            ax.set_ylabel("Nucleus Area")
            ax.set_title(f"{sample_names[0]}: Cell Area vs Nucleus Area")
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
        else:
            nrows = int(np.ceil(n / 3))
            fig, axes = plt.subplots(nrows, 3, figsize=(8, 6 * nrows))
            axes = axes.flatten()
            for i, adata in enumerate(adata_list):
                axes[i].scatter(adata.obs['cell_area'], adata.obs['nucleus_area'], alpha=0.5)
                axes[i].set_xlabel("Cell Area")
                axes[i].set_ylabel("Nucleus Area")
                axes[i].set_title(f"{sample_names[i]}: Cell Area vs Nucleus Area")
            for j in range(i+1, nrows*3):
                fig.delaxes(axes[j])
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Saved cell vs nucleus area plots to {pdf_filename}")

def plot_area_ratio_distributions(adata_list, sample_names, pdf_filename="06_area_ratio_distributions.pdf"):
    all_area_ratios = []
    all_samples = []
    for adata, name in zip(adata_list, sample_names):
        if 'area_ratio' not in adata.obs.columns:
            adata.obs['area_ratio'] = adata.obs['nucleus_area'] / adata.obs['cell_area']
        all_area_ratios.append(adata.obs['area_ratio'])
        all_samples.extend([name] * adata.shape[0])
    all_area_ratios = np.concatenate(all_area_ratios)
    df = pd.DataFrame({'area_ratio': all_area_ratios, 'sample': all_samples})
    with PdfPages(pdf_filename) as pdf:
        plt.figure(figsize=(8, 6))
        sns.violinplot(x='sample', y='area_ratio', hue='sample', data=df, legend=False)
        plt.xlabel("Sample")
        plt.ylabel("Nucleus-to-Cell Area Ratio")
        plt.title("Nucleus-to-Cell Area Ratio Distribution by Sample")
        plt.tight_layout()
        pdf.savefig()
        plt.close()
    print(f"Saved area ratio distributions plot to {pdf_filename}")

def plot_cell_area_vs_total_counts(adata_list, sample_names, pdf_filename="07_cell_area_vs_total_counts.pdf"):
    n = len(adata_list)
    with PdfPages(pdf_filename) as pdf:
        if n == 1:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(adata_list[0].obs['cell_area'], adata_list[0].obs['total_counts'], alpha=0.5)
            ax.set_xlabel("Cell Area")
            ax.set_ylabel("Total Transcript Counts")
            ax.set_title(f"{sample_names[0]}: Cell Area vs Total Transcript Counts")
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
        else:
            nrows = int(np.ceil(n / 3))
            fig, axes = plt.subplots(nrows, 3, figsize=(8, 6 * nrows))
            axes = axes.flatten()
            for i, adata in enumerate(adata_list):
                axes[i].scatter(adata.obs['cell_area'], adata.obs['total_counts'], alpha=0.5)
                axes[i].set_xlabel("Cell Area")
                axes[i].set_ylabel("Total Transcript Counts")
                axes[i].set_title(f"{sample_names[i]}: Cell Area vs Total Transcript Counts")
            for j in range(i+1, nrows*3):
                fig.delaxes(axes[j])
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Saved cell area vs total counts plots to {pdf_filename}")

def plot_spatial_cell_and_nucleus_area(adata_list, sample_names, spot_size=100, cmap='viridis_r', pdf_filename="08_spatial_cell_nucleus_area.pdf"):
    with PdfPages(pdf_filename) as pdf:
        for adata, name in zip(adata_list, sample_names):
            fig, axarr = plt.subplots(1, 2, figsize=(16, 6))
            sc.pl.spatial(adata, color='cell_area', spot_size=spot_size, title=f"{name}: Cell Area", cmap=cmap, ax=axarr[0], show=False)
            sc.pl.spatial(adata, color='nucleus_area', spot_size=spot_size, title=f"{name}: Nucleus Area", cmap=cmap, ax=axarr[1], show=False)
            plt.suptitle(f"{name}: Spatial Distribution of Cell and Nucleus Area", fontsize=14)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Saved spatial cell and nucleus area plots to {pdf_filename}")

def plot_spatial_gene_expression(adata_list, sample_names, gene, spot_size=100, cmap='viridis_r', pdf_filename="09_spatial_gene_expression.pdf"):
    with PdfPages(pdf_filename) as pdf:
        for adata, name in zip(adata_list, sample_names):
            fig, axes = plt.subplots(1, 2, figsize=(8, 6))
            sc.pl.spatial(adata, color='total_counts', spot_size=spot_size, title=f"{name}: Total Counts", cmap=cmap, ax=axes[0], show=False)
            sc.pl.spatial(adata, color=[gene], spot_size=spot_size, title=f"{name}: {gene} Expression", cmap=cmap, ax=axes[1], show=False)
            plt.suptitle(f"{name}: Total counts and {gene} Expression", fontsize=14)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Saved spatial gene expression plots to {pdf_filename}")

def plot_cp_uc(adata_list, sample_names, spot_size=100, cmap='viridis_r', pdf_filename="10_cp_uc.pdf"):
    n = len(adata_list)
    with PdfPages(pdf_filename) as pdf:
        fig, axes = plt.subplots(n, 3, figsize=(18, 6 * n), constrained_layout=True)
        if n == 1:
            axes = np.expand_dims(axes, 0)
        for i, (adata, name) in enumerate(zip(adata_list, sample_names)):
            sns.histplot(adata.obs['control_probe_counts'], kde=True, bins=50, ax=axes[i, 0])
            axes[i, 0].set_xlabel("Control Probe Counts")
            axes[i, 0].set_ylabel("Number of Cells")
            axes[i, 0].set_title(f"{name}: Control Probe Counts")
            sc.pl.spatial(adata, color='control_probe_counts', spot_size=spot_size,
                          title=None, cmap=cmap, ax=axes[i, 1], show=False)
            axes[i, 1].set_title(f"{name}: Spatial Control Probe Counts")
            sc.pl.spatial(adata, color='unassigned_codeword_counts', spot_size=spot_size,
                          title=None, cmap=cmap, ax=axes[i, 2], show=False)
            axes[i, 2].set_title(f"{name}: Spatial Unassigned Codeword Counts")
        pdf.savefig(fig)
        plt.close(fig)
    print(f"Saved control probe/unassigned codeword plots to {pdf_filename}")

def plot_extreme_counts(adata_list, sample_names, spot_size, pdf_filename="11_extreme_counts.pdf"):
    all_counts = []
    all_area_ratios = []
    all_samples = []
    with PdfPages(pdf_filename) as pdf:
        for adata, name in zip(adata_list, sample_names):
            adata.obs['min10'] = adata.obs['total_counts'] < 10
            adata.obs['min30'] = adata.obs['total_counts'] < 30
            adata.obs['max350'] = adata.obs['total_counts'] > 350
            adata.obs['extreme_area_ratio'] = (adata.obs['area_ratio'] >= 0.9) | (adata.obs['area_ratio'] <= 0.1)
            all_counts.append(adata.obs['total_counts'])
            all_area_ratios.append(adata.obs['area_ratio'])
            all_samples.extend([name] * adata.shape[0])
            fig = sc.pl.spatial(
                adata,
                color=['min10', 'min30', 'max350', 'extreme_area_ratio'],
                spot_size=spot_size,
                palette=["lightgray", "blue"],
                title=[f"{name}: min10", f"{name}: min30", f"{name}: max350", f"{name}: extreme_area_ratio"],
                show=False,
                return_fig=True
            )
            pdf.savefig(fig)
            plt.close(fig)
        all_counts = np.concatenate(all_counts)
        all_area_ratios = np.concatenate(all_area_ratios)
        df = pd.DataFrame({
            'total_counts': all_counts,
            'area_ratio': all_area_ratios,
            'sample': all_samples
        })
        fig, axes = plt.subplots(1, 2, figsize=(14, 8))
        sns.violinplot(x='sample', y='total_counts', hue='sample', legend=False, data=df, ax=axes[0])
        axes[0].set_title("Total Counts Distribution by Sample")
        axes[0].set_xlabel("Sample")
        axes[0].set_ylabel("Total Counts")
        sns.violinplot(x='sample', y='area_ratio', data=df, ax=axes[1])
        axes[1].set_title("Area Ratio Distribution by Sample")
        axes[1].set_xlabel("Sample")
        axes[1].set_ylabel("Area Ratio")
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    print(f"Saved extreme counts plots to {pdf_filename}")

def filter_adata(
    adata,
    min_counts=30,
    max_counts=350,
    min_genes=5,
    max_genes=110,
    min_cell_area=20,
    max_cell_area=900,
    max_area_ratio=0.9,
    min_area_ratio=0.1,
    max_nucleus_area=110,
    min_cells_per_gene=5,
    pdf_filename="12_filtered_cells_spatial.pdf"
):
    initial_cells_count = adata.n_obs
    initial_genes_count = adata.n_vars
    original_indices = adata.obs_names.copy()
    mask = (
        (adata.obs['total_counts'] >= min_counts) &
        (adata.obs['total_counts'] <= max_counts) &
        (adata.obs['n_genes'] >= min_genes) &
        (adata.obs['n_genes'] <= max_genes) &
        (adata.obs['cell_area'] > min_cell_area) &
        (adata.obs['cell_area'] < max_cell_area) &
        (adata.obs['area_ratio'] < max_area_ratio) &
        (adata.obs['area_ratio'] > min_area_ratio) &
        (adata.obs['nucleus_area'] < max_nucleus_area)
    )
    adata.obs['filtered_out'] = ~mask
    with PdfPages(pdf_filename) as pdf:
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.spatial(
            adata,
            color='filtered_out',
            spot_size=100,
            palette=["lightgray", "red"],
            title="Cells Filtered Out (red)",
            ax=ax,
            show=False
        )
        pdf.savefig(fig)
        plt.close(fig)
    sc.pp.filter_cells(adata, min_counts=min_counts, inplace=True)
    sc.pp.filter_cells(adata, max_counts=max_counts, inplace=True)
    sc.pp.filter_cells(adata, min_genes=min_genes, inplace=True)
    sc.pp.filter_cells(adata, max_genes=max_genes, inplace=True)
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene, inplace=True)
    adata = adata[adata.obs['cell_area'] > min_cell_area, :]
    adata = adata[adata.obs['cell_area'] < max_cell_area, :]
    adata = adata[adata.obs['area_ratio'] < max_area_ratio, :]
    adata = adata[adata.obs['area_ratio'] > min_area_ratio, :]
    adata = adata[adata.obs['nucleus_area'] < max_nucleus_area, :]
    filtered_cells_count = adata.n_obs
    filtered_genes_count = adata.n_vars
    print(f"Filtered {initial_cells_count - filtered_cells_count} (out of initial {initial_cells_count} cells)")
    print(f"Filtered {initial_genes_count - filtered_genes_count} (out of initial {initial_genes_count} genes)")
    print(f"Saved filtered cells spatial plot to {pdf_filename}")
    return adata

def plot_filtered_distributions(adata_list, sample_names, pdf_filename="13_filtered_distributions.pdf"):
    all_counts = []
    all_nucleus_areas = []
    all_area_ratios = []
    all_samples = []
    for adata, name in zip(adata_list, sample_names):
        all_counts.append(adata.obs['total_counts'])
        all_nucleus_areas.append(adata.obs['nucleus_area'])
        all_area_ratios.append(adata.obs['area_ratio'])
        all_samples.extend([name] * adata.shape[0])
    all_counts = np.concatenate(all_counts)
    all_nucleus_areas = np.concatenate(all_nucleus_areas)
    all_area_ratios = np.concatenate(all_area_ratios)
    df = pd.DataFrame({
        'total_counts': all_counts,
        'nucleus_area': all_nucleus_areas,
        'area_ratio': all_area_ratios,
        'sample': all_samples
    })
    with PdfPages(pdf_filename) as pdf:
        fig, axs = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
        sns.violinplot(x='sample', y='total_counts', hue='sample', legend=False, data=df, ax=axs[0])
        axs[0].set_title("Total Transcript Counts per Cell After Filtering")
        axs[0].set_xlabel("Sample")
        axs[0].set_ylabel("Total Counts")
        sns.violinplot(x='sample', y='nucleus_area', hue='sample', legend=False, data=df, ax=axs[1])
        axs[1].set_title("Nucleus Area After Filtering")
        axs[1].set_xlabel("Sample")
        axs[1].set_ylabel("Nucleus Area")
        sns.violinplot(x='sample', y='area_ratio', hue='sample', legend=False, data=df, ax=axs[2])
        axs[2].set_title("Nucleus-to-Cell Area Ratios After Filtering")
        axs[2].set_xlabel("Sample")
        axs[2].set_ylabel("Area Ratio")
        pdf.savefig(fig)
        plt.close(fig)
    print(f"Saved filtered distributions plots to {pdf_filename}")
# Violin plot of total_counts, n_genes, cell_area, nucleus_area, area_ratio grouped based on the regions
def plot_violin_by_region(
    adata_list, 
    sample_names, 
    region_column='all_regions', 
    features=['total_counts', 'n_genes', 'cell_area', 'nucleus_area', 'area_ratio'],
    pdf_filename="14_violin_by_region.pdf"
):

    n_samples = len(adata_list)
    n_features = len(features)
    with PdfPages(pdf_filename) as pdf:
        fig, axes = plt.subplots(n_samples, n_features, figsize=(8 * n_features, 6 * n_samples), squeeze=False)

        for i, (adata, sample_name) in enumerate(zip(adata_list, sample_names)):
            for j, feature in enumerate(features):
                ax = axes[i, j]
                sns.violinplot(
                    x=region_column, 
                    y=feature, 
                    hue=region_column, 
                    legend=False, 
                    data=adata.obs, 
                    ax=ax
                )
                ax.set_title(f"{sample_name}: {feature} by {region_column}")
                ax.set_xlabel(region_column)
                ax.set_ylabel(feature)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    print(f"Saved violin plots by region to {pdf_filename}")

def plot_normalization_effects(adata_list, sample_names, pdf_filename="15_normalization_effects.pdf"):
    with PdfPages(pdf_filename) as pdf:
        for adata, name in zip(adata_list, sample_names):
            original_counts = adata.layers['raw'].sum(axis=1)
            normalized_counts = adata.X.sum(axis=1)
            plt.figure(figsize=(8, 6))
            sns.histplot(original_counts.A1, color="blue", label="Before Normalization", kde=True)
            sns.histplot(normalized_counts.A1, color="orange", label="After Normalization", kde=True)
            plt.xlabel("Total Expression per Cell")
            plt.ylabel("Number of Cells")
            plt.title(f"Effect of Normalization on Expression Distribution: {name}")
            plt.legend()
            pdf.savefig()
            plt.close()
    print(f"Saved normalization effect plots to {pdf_filename}")

def plot_region_annotations_to_pdf(adata_list, sample_names, pdf_filename="16_region_annotations.pdf", subsample_fraction=0.1):
    with PdfPages(pdf_filename) as pdf:
        for adata, name in zip(adata_list, sample_names):
            fig, axs = plt.subplots(1, 3, figsize=(24, 5))
            adata_sub = sc.pp.subsample(adata, fraction=subsample_fraction, copy=True)
            sc.pl.spatial(adata_sub, color="all_regions", spot_size=100, ax=axs[0], show=False)
            sc.pl.spatial(adata_sub, color="region_invasive", spot_size=100, ax=axs[1], show=False)
            sc.pl.spatial(adata_sub, color="region_extended_invasive", spot_size=100, ax=axs[2], show=False)
            plt.suptitle(f"Sample: {name}", fontsize=20)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            print(f"Saved region annotations for {name} to {pdf_filename}")

# Function to concat all PDFs into one and save as sampleName_qc.pdf
def concat_pdfs(sampleName, folder):
    pdf_files = glob.glob(os.path.join(folder, f"*.pdf"))
    # Sort the PDF files to ensure consistent order, using numerical sorting
    pdf_files.sort(key=lambda x: int(os.path.basename(x).split('_')[0]) if '_' in os.path.basename(x) else x)
    if not pdf_files:
        print("No PDF files found to concatenate.")
        return
    output_pdf = os.path.join(folder, f"{sampleName}_qc.pdf")
    merger = PdfMerger()
    for pdf_file in pdf_files:
        if pdf_file != output_pdf:  # Avoid including the output file if it already exists
            merger.append(pdf_file)
    merger.write(output_pdf)
    merger.close()
    print(f"Concatenated PDFs saved to {output_pdf}")

sampleNames = ["${sampleName}"]
h5ad = "${h5ad}"
reg_col = "${reg_col}"
min_counts = int("${min_counts}")
max_counts = int("${max_counts}")
min_genes = int("${min_genes}")
max_genes = int("${max_genes}")
min_cell_area = int("${min_cell_area}")
max_cell_area = int("${max_cell_area}")
min_area_ratio = float("${min_area_ratio}")
max_area_ratio = float("${max_area_ratio}")
min_cells_per_gene = int("${min_cells_per_gene}")
max_nucleus_area = int("${max_nucleus_area}")
outFile = f"{sampleNames[0]}_qc.h5ad"
adataList = [sc.read_h5ad(h5ad)]

plot_counts_per_cell(adata_list=adataList,sample_names=sampleNames)
plot_ngenes_per_cell(adata_list=adataList,sample_names=sampleNames)
plot_highest_expr_genes(adata_list=adataList, sample_names=sampleNames)
plot_area_distributions(adata_list=adataList, sample_names=sampleNames)
plot_cell_vs_nucleus_area(adata_list=adataList, sample_names=sampleNames)
plot_area_ratio_distributions(adata_list=adataList, sample_names=sampleNames)
plot_cell_area_vs_total_counts(adata_list=adataList, sample_names=sampleNames)
plot_spatial_cell_and_nucleus_area(adata_list=adataList, sample_names=sampleNames)
plot_cp_uc(adata_list=adataList, sample_names=sampleNames, spot_size=100, cmap='viridis_r')
plot_extreme_counts(adata_list=adataList, sample_names=sampleNames, spot_size=100)
# If any of the filter parameters are set to -1, skip filtering
if (any(param == -1 for param in [
    min_counts, max_counts, min_genes, max_genes,
    min_cell_area, max_cell_area, min_area_ratio, max_area_ratio,
    max_nucleus_area, min_cells_per_gene
])):
    print("Skipping filtering due to -1 parameters.")
else:
    print("Filtering adata based on provided parameters...")
    # Filter the adata object based on the provided parameters
    adataList[0] = filter_adata(
        adata=adataList[0],
        min_counts=min_counts,
        max_counts=max_counts,
        min_genes=min_genes,
        max_genes=max_genes,
        min_cell_area=min_cell_area,
        max_cell_area=max_cell_area,
        min_area_ratio=min_area_ratio,
        max_area_ratio=max_area_ratio,
        max_nucleus_area=max_nucleus_area,
        min_cells_per_gene=min_cells_per_gene
    )
    plot_filtered_distributions(adata_list=adataList, sample_names=sampleNames)
    plot_highest_expr_genes(adata_list=adataList, sample_names=sampleNames)
adataList[0].layers['raw'] = adataList[0].X.copy()
sc.pp.normalize_total(adataList[0], target_sum=1e6, inplace=True)
sc.pp.log1p(adataList[0])
plot_normalization_effects(adata_list=adataList, sample_names=sampleNames)
plot_violin_by_region(
    adata_list=adataList, 
    sample_names=sampleNames, 
    region_column=reg_col
)

concat_pdfs(sampleName=sampleNames[0], folder="./")
# Save the final adata object
adataList[0].write_h5ad(outFile)
print("QC analysis completed successfully.")
