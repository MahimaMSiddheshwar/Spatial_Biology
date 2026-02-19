
METHODS
=======

Dataset
-------
Publicly available 10x Genomics Visium spatial transcriptomics data
from human breast cancer fresh-frozen tissue (v1.3.0) was obtained
from the 10x Genomics dataset portal. The dataset comprised
4,869 tissue-covered spots with a median of 9,720 UMIs
and 3,654 genes detected per spot.

Quality Control
---------------
Spots with fewer than 200 detected genes, fewer than 500 total UMI
counts, or greater than 15% mitochondrial gene expression were
excluded. Genes detected in fewer than 3 spots were removed,
yielding 4,869 spots and 21,349 genes for
downstream analysis. Hemoglobin gene contamination was negligible
(mean HB% = 0.01%).

Normalization and Preprocessing
--------------------------------
Raw counts were preserved in adata.layers['counts']. Library sizes
were normalized to 10,000 counts per spot followed by log1p
transformation. The top 3,000 highly variable genes were selected
using the Seurat v3 method. Principal component analysis (PCA) was
performed on scaled HVGs (50 components), and a k-nearest neighbor
graph (k=15, 30 PCs) was constructed for UMAP embedding.

Clustering and Annotation
--------------------------
Leiden community detection (resolution=0.5) identified 9
transcriptionally distinct spatial clusters. Clusters were annotated
using known breast cancer marker genes from CellMarker 2.0 and
PanglaoDB databases. Spatially variable genes were identified using
Moran's I statistic (n_perms=100) with Benjamini-Hochberg FDR
correction, yielding 2,486 significant SVGs (FDR < 0.05).

Cell Type Deconvolution
------------------------
Cell type composition per spot was estimated using gene signature
scoring (sc.tl.score_genes, Scanpy v1.12) with 24
curated signatures derived from Wu et al. (2021) Nature Genetics
breast cancer single-cell RNA-seq atlas. Signatures covered tumor
epithelial, stromal, immune, and special functional categories
including TLS, hypoxia, EMT, and T cell exhaustion.

Cell-Cell Communication
-----------------------
Ligand-receptor co-expression scores were computed as the geometric
mean of ligand and receptor expression per spot. A total of 39 LR
pairs from CellChatDB and NicheNet databases covering 8 signaling
pathways were analyzed. Statistical significance was assessed using
permutation testing (n=500) implemented in Squidpy sq.gr.ligrec.

PAM50 Subtyping
---------------
Molecular subtype scores were computed for all 5 PAM50 subtypes
(Luminal A, Luminal B, HER2-enriched, Basal-like, Normal-like)
using curated gene signatures and spatially mapped across tissue.

Software
--------
All analyses were performed in Python 3.13 using Scanpy v1.12,
Squidpy v1.8.1, NumPy, Pandas, and Matplotlib.
Clustering used the leidenalg package. Visualizations were generated
with Matplotlib and Squidpy spatial plotting functions.

References
----------
1. Wolf et al. (2018) Scanpy. Genome Biology.
2. Palla et al. (2022) Squidpy. Nature Methods.
3. Wu et al. (2021) Breast cancer scRNA atlas. Nature Genetics.
4. Traag et al. (2019) Leiden algorithm. Scientific Reports.
5. Jin et al. (2021) CellChat. Nature Communications.
