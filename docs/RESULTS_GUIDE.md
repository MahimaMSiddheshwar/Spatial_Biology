# Spatial Biology Results Guide

This document explains how to read the outputs of each step and what conclusions you can draw.

## Step 1 Results: Single-Cell Analysis

### Files to review

- `data/processed/human_breast_cancer_single_cell.h5ad`
- `results/single_cell/*_clusters.csv`
- QC and UMAP figures in `figures/`

### How to interpret

- **QC plot (`n_counts`, `n_genes`)**
  - Very low values can indicate poor-quality spots/cells
  - Extremely high counts may indicate technical artifacts
- **UMAP with Leiden clusters**
  - Clear separated groups suggest distinct cell populations
  - Heavy overlap suggests similar transcriptional states
- **Cluster counts CSV**
  - Shows relative abundance of discovered populations

### Expected insight

- Initial map of tumor microenvironment cell populations and heterogeneity.

## Step 2 Results: Machine Learning

### Files to review

- `results/machine_learning/gene_importance.csv`
- `results/machine_learning/cell_type_classifier.pkl`
- `figures/top_biomarkers.png`

### How to interpret

- **Classification report (console output)**
  - Focus on precision, recall, F1 for each cell type
  - Lower scores indicate classes that are hard to separate
- **Gene importance CSV/plot**
  - Top genes are candidate biomarkers driving class separation
  - Use as shortlist for downstream validation

### Expected insight

- Candidate biomarkers and a reusable classifier for cell-type prediction.

## Step 2B Results (Optional): scVI Analysis

### Files to review

- `data/processed/human_breast_cancer_scvi.h5ad`
- `results/machine_learning/scvi_model/`
- `results/machine_learning/scvi_latent_embeddings.csv`
- `results/machine_learning/scvi_cluster_counts.csv`
- `figures/human_breast_cancer_scvi_umap.png`

### How to interpret

- **scVI UMAP (`leiden_scvi`)**
  - Better separation can indicate cleaner latent structure than baseline PCA space
- **Latent embeddings CSV**
  - Can be used for downstream tasks (classification, trajectory, nearest-neighbor analysis)
- **scVI cluster counts**
  - Shows abundance distribution of latent-space clusters

### Expected insight

- A deep-learning-based view of cellular structure that complements Random Forest results.

## Step 3 Results: Spatial Analysis

### Files to review

- `data/processed/human_breast_cancer_spatial.h5ad`
- `results/spatial/*_spatial_clusters.csv`
- Spatial plots in `figures/`

### How to interpret

- **Spatial cluster plot**
  - Neighboring spots with same label indicate tissue domains
- **Spatial autocorrelation output**
  - Higher spatial autocorrelation suggests non-random local structure
- **Spatial cluster counts CSV**
  - Quantifies size of each tissue domain

### Expected insight

- Spatial architecture of tumor regions and local microenvironment organization.

## Step 4 Results: Imaging Integration

### Files to review

- `results/imaging/human_breast_cancer_integration_report.csv`
- `figures/Human Breast Cancer_multi_modal.png`

### How to interpret

- **Overlay figure**
  - Confirms whether molecular spatial clusters align with histology patterns
- **Integration report CSV**
  - Quick technical summary (spots, genes, image shape, completion status)

### Expected insight

- Visual confirmation of how transcriptional regions map to tissue morphology.

## Quality Checks Before Moving to the Next Step

- Step 1 -> Step 2: ensure `.h5ad` from Step 1 exists
- Step 2/2B -> Step 3: ML outputs generated without errors
- Step 3 -> Step 4: spatial `.h5ad` and tissue image are present
