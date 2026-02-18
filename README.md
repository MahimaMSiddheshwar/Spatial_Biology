# Visium Spatial Transcriptomics Workflow (Single Sample)

This project is organized for a complete Visium Human Breast Cancer analysis workflow, from raw data loading to final biological interpretation.

## Folder Structure

```text
spatial_biology_project/
|-- data/
|   |-- raw/                          # Input Visium files (already provided)
|   `-- processed/                    # Intermediate and final .h5ad outputs
|-- scripts/                          # Notebooks (01 to 06)
|-- results/
|   |-- qc/                           # QC summaries and filtering reports
|   |-- clustering/                   # Cluster labels and marker tables
|   |-- spatial/                      # Spatial variable genes and domains
|   |-- deconvolution/                # Cell type proportion outputs
|   |-- machine_learning/             # Supervised/unsupervised ML results
|   |-- imaging/                      # Image-feature and integration outputs
|   `-- interpretation/               # Final biological summary tables
`-- figures/
    |-- qc/                           # QC plots on tissue and distributions
    |-- clustering/                   # UMAP + tissue cluster overlays
    |-- spatial/                      # SVG and spatial domain plots
    |-- machine_learning/             # ML performance and feature plots
    |-- imaging/                      # Patch/image integration plots
    `-- interpretation/               # Publication-ready summary figures
```

## Workflow (Notebook Plan)

1. `01_data_loading_qc.ipynb`
   - Load data with `sc.read_visium()`
   - Compute QC metrics and apply filtering
   - Normalize, log-transform, set `adata.raw`, select HVGs once
   - Save: `data/processed/visium_qc_normalized.h5ad`

2. `02_clustering.ipynb`
   - PCA, neighbors, UMAP, Leiden clustering (multiple resolutions)
   - Plot clusters on UMAP and tissue image
   - Marker gene discovery and cluster annotation
   - Save: `data/processed/visium_clustered.h5ad`

3. `03_spatial_analysis.ipynb`
   - Spatial neighbors and Moran's I (SVG ranking)
   - Spatial domain detection (baseline + STAGATE optional)
   - Cell type deconvolution (Cell2Location/SPOTlight)
   - Save: `data/processed/visium_deconvolved.h5ad`

4. `04_machine_learning.ipynb`
   - Feature engineering (PCA, scVI, neighborhood/deconvolution features)
   - Unsupervised and supervised ML (RF/XGBoost)
   - Save ML metrics and feature importances

5. `05_image_integration.ipynb`
   - Tissue image preprocessing and patch extraction
   - Image feature extraction (CNN or Squidpy)
   - Multi-modal integration with expression features

6. `06_biological_interpretation.ipynb`
   - Integrate all outputs into tumor microenvironment maps
   - Produce publication figures and final exports
   - Save: `data/processed/visium_final.h5ad`

## Core Rules

- Always use `sc.read_visium()` (not `sc.read_10x_h5()`).
- Perform HVG selection once on the full normalized gene set.
- Keep `adata.raw = adata` before subsetting.
- Plot clusters on both UMAP and tissue image.
- Save an intermediate `.h5ad` after each major notebook.
