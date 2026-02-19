# Complete Workflow

This document describes the end-to-end workflow implemented by the notebooks in `scripts/`.

## What You Get

- A fully reproducible path from raw Visium sample downloads to a clustered + annotated AnnData object.
- QC, preprocessing, clustering, marker discovery, and publication-ready spatial/UMAP figures.

## Prerequisites

- Python 3.10+ recommended.
- `wget` and `tar` available (Notebook 01 uses them). On Windows, the easiest options are:
  - Git Bash (includes `wget`/`tar` if installed) or
  - WSL, or
  - conda-installed `wget`.

Core Python packages used:
- `scanpy`
- `squidpy`
- `anndata`, `numpy`, `pandas`, `matplotlib`

## Repository Conventions

- Raw data lives in `data/raw/Visium_Human_Breast_Cancer/` (ignored by git).
- Processed checkpoints live in `data/processed/` (ignored by git).
- Figures are written to `figures/` (tracked).
- Results tables can be written to `results/` (tracked).

## Notebook 01 - Data Setup

Path: `scripts/01_Data setup.ipynb`

Goal:
- Download the 10x Visium Human Breast Cancer v1.3.0 sample.
- Load Visium data using `scanpy.read_visium()`.
- Confirm the slide image/coordinates are available and filter to tissue spots.

What it downloads (to `data/raw/Visium_Human_Breast_Cancer/`):
- `Visium_Human_Breast_Cancer_spatial.tar.gz`
- `Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5`
- `Visium_Human_Breast_Cancer_raw_feature_bc_matrix.h5` (optional, kept for completeness)
- `Visium_Human_Breast_Cancer_molecule_info.h5` (optional)
- `Visium_Human_Breast_Cancer_analysis.tar.gz` (optional)

Key checks:
- `adata.obsm['spatial']` exists and has shape `(n_spots, 2)`.
- Filter to tissue spots:
  - `adata = adata[adata.obs['in_tissue'] == 1].copy()`

Outputs:
- No checkpoint is written here (Notebook 02 produces the first checkpoint).

## Notebook 02 - QC and Preprocessing

Path: `scripts/02_QC & Preprocessing.ipynb`

Goal:
- Compute QC metrics.
- Apply filtering thresholds.
- Normalize + log-transform.
- Compute HVGs, PCA, neighbors, UMAP.
- Save a clean preprocessed checkpoint.

What it does (high level):

1) Load Visium data from `data/raw/Visium_Human_Breast_Cancer/`.
2) Filter to tissue spots (`in_tissue == 1`).
3) QC metrics:
- mitochondrial (`MT-`), ribosomal (`RPS`/`RPL`), hemoglobin (`HB*`)
- `sc.pp.calculate_qc_metrics(... qc_vars=['mt','ribo','hb'])`
4) Filtering (values are intended for this dataset; adjust if you switch samples):
- `min_genes` (e.g. 200)
- `min_counts` (e.g. 500)
- `min_cells` (e.g. 3)
- `pct_counts_mt` cutoff (e.g. < 15)
- Recompute QC metrics after filtering (recommended).
5) Preserve counts and normalized data:
- Save raw counts to `adata.layers['counts']`.
- Run `sc.pp.normalize_total(... target_sum=1e4)` and `sc.pp.log1p(adata)`.
- Save log-normalized snapshot to `adata.raw`.
6) HVGs:
- `sc.pp.highly_variable_genes(... flavor='seurat_v3', layer='counts')`
7) PCA, neighbors, UMAP.

Outputs:
- Checkpoint: `data/processed/adata_preprocessed.h5ad`
- Figures:
  - `figures/qc/` (violins, scatter QC, HVGs, UMAP QC)
  - `figures/spatial/` (spatial QC overlays)

## Notebook 03 - Clustering and Annotation

Path: `scripts/03_Clustering.ipynb`

Goal:
- Cluster the dataset (Leiden), choose a sensible resolution.
- Find markers per cluster.
- Overlay marker genes on tissue/UMAP.
- Create human-readable cluster labels (`adata.obs['cell_type']`).
- Save an annotated checkpoint.

Main steps:

1) Load checkpoint:
- `data/processed/adata_preprocessed.h5ad`

2) Leiden clustering (multiple resolutions):
- Typically tries: 0.3, 0.5, 0.6, 1.0
- Decide best resolution by:
  - spatial coherence on tissue (contiguous regions),
  - marker gene interpretability,
  - avoiding over-fragmentation.
- For this dataset, `leiden_0.5` is a good default.

3) Marker genes per cluster:
- `sc.tl.rank_genes_groups(... method='wilcoxon', groupby='leiden', use_raw=True)`
- Dotplot for top markers.

4) Known marker overlays (examples):
- epithelial/tumor: `EPCAM`, `KRT*`, `ESR1`, `ERBB2`, `MKI67`
- stroma/CAF: `COL1A1`, `VIM`
- immune: `CD68`, `CD3D`
- endothelial: `PECAM1`

5) Annotation:
- Create a `cluster_labels` mapping and write `adata.obs['cell_type']`.

Outputs:
- Checkpoint: `data/processed/adata_annotated.h5ad`
- Figures:
  - `figures/clustering/` (UMAP overview, annotated UMAP panels)
  - `figures/spatial/` (Leiden overlays, marker overlays, annotated spatial panel)
  - `figures/markers/` (dotplots)

## Common Issues and Fixes

- NameError for paths like `fig_cl_dir` / `fig_sp_dir` / `fig_mk_dir`:
  - You ran a later cell without running the setup cell (or restarted the kernel). Re-run the first "IMPORTS & PATHS" cell.

- KeyError: 'leiden' when ranking genes:
  - Run the "CHOOSE BEST RESOLUTION" cell first (it creates `adata.obs['leiden']`).

- Black background in saved PNGs:
  - Transparent PNGs can look black in Windows Photos. The notebooks save non-transparent PNGs (white background) for spatial/UMAP panels.

## Next Notebooks

The `results/` and `figures/` folders include placeholders for downstream work (spatial statistics, deconvolution, imaging, interpretation). If you want, we can extend the workflow with:
- Spatially variable genes (Moran's I)
- Spatial domains (e.g. STAGATE)
- Deconvolution (Cell2Location / SPOTlight)
