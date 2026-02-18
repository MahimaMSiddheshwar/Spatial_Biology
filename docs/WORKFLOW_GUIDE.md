# Spatial Biology Workflow Guide

This guide explains what each pipeline step does, why it is needed, and the order to run the project.

## What You Need Before Running

### 1) Python environment

- Python 3.10+ recommended
- Create and activate a virtual environment

### 2) Required libraries

- Core: `scanpy`, `squidpy`, `numpy`, `pandas`, `matplotlib`, `scikit-learn`, `joblib`
- Optional for Step 2B: `scvi-tools`
- Install example:

```bash
pip install scanpy squidpy numpy pandas matplotlib scikit-learn joblib
```

Optional Step 2B install:

```bash
pip install scvi-tools
```

### 3) Required dataset files

Keep these in `data/raw/`:

- `Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5`
- `spatial/` folder containing files such as:
  - `tissue_hires_image.png`
  - `tissue_lowres_image.png`
  - `scalefactors_json.json`
  - `tissue_positions_list.csv`

## Run Order (One by One)

Run from project root:

```bash
jupyter lab
```

Do not run all steps together. Run one notebook, check outputs, then continue.
Open notebooks from `scripts/` in this order:

1. `01_single_cell_analysis.ipynb`
2. `02_machine_learning.ipynb`
3. `02b_scvi_analysis.ipynb` (optional)
4. `03_spatial_analysis.ipynb`
5. `04_imaging_integration.ipynb`

## What Each Step Does and Why

### Step 1: `01_single_cell_analysis.ipynb`

- **What it does:** loads expression matrix, performs QC, normalization, HVG selection, PCA, clustering, UMAP
- **Why it matters:** creates clean, biologically meaningful cell groups used by later steps
- **Main outputs:**
  - `data/processed/human_breast_cancer_single_cell.h5ad`
  - `results/single_cell/*`
  - `figures/*` (QC and UMAP plots)

### Step 2: `02_machine_learning.ipynb`

- **What it does:** trains a Random Forest to classify cell types from processed expression data
- **Why it matters:** provides automated labeling and biomarker ranking
- **Main outputs:**
  - `results/machine_learning/gene_importance.csv`
  - `results/machine_learning/cell_type_classifier.pkl`
  - `figures/top_biomarkers.png`

### Step 2B (Optional): `02b_scvi_analysis.ipynb`

- **What it does:** trains an scVI deep generative model and computes latent-space clusters/UMAP
- **Why it matters:** improves representation learning and can capture complex biological structure
- **Main outputs:**
  - `data/processed/human_breast_cancer_scvi.h5ad`
  - `results/machine_learning/scvi_model/`
  - `results/machine_learning/scvi_latent_embeddings.csv`
  - `results/machine_learning/scvi_cluster_counts.csv`
  - `figures/human_breast_cancer_scvi_umap.png`

### Step 3: `03_spatial_analysis.ipynb`

- **What it does:** computes spatial neighbors, spatial clusters, and spatial autocorrelation
- **Why it matters:** links gene patterns to tissue location and architecture
- **Main outputs:**
  - `data/processed/human_breast_cancer_spatial.h5ad`
  - `results/spatial/*`
  - `figures/*` (spatial clustering/domain plots)

### Step 4: `04_imaging_integration.ipynb`

- **What it does:** overlays spatial transcriptomics outputs on tissue image
- **Why it matters:** creates interpretable multi-modal visualization for biological context
- **Main outputs:**
  - `results/imaging/human_breast_cancer_integration_report.csv`
  - `figures/Human Breast Cancer_multi_modal.png`
