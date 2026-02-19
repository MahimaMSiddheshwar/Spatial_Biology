# Biological Interpretation of Your Leiden res=0.5 Map

| Cluster | Color | Likely Identity | Evidence |
|---:|---|---|---|
| 0 | 🔵 Blue | Normal breast / stromal tissue | Large homogeneous lower region |
| 1 | 🟠 Orange | Tumor epithelial | Scattered across upper tissue |
| 2 | 🟢 Green | Stromal / CAF-rich region | Upper right, diffuse |
| 3 | 🔴 Red | Tumor core | Dense, spatially confined — classic carcinoma nest |
| 4 | 🟣 Purple | Mixed stroma / border zone | Surrounds tumor regions |
| 5 | 🟤 Brown | Immune-rich / TIL hotspot | Small isolated regions |
| 6 | 🩷 Pink | Invasive margin / transition | Thin boundary regions |
| 7 | 🟡 Olive | Rare cell type / adipose | Very small population |
| 8 | 🩵 Cyan | Endothelial / vascular | Tiny — blood vessel spots |



# What the marker maps reveal

| Marker | Pattern | Biological Meaning |
|---|---|---|
| EPCAM | Hot spots upper-left + center | ✅ Tumor epithelial nests — confirms cluster 3 (red) is carcinoma |
| ESR1 | Strong upper region | ✅ Luminal A subtype — ER+ breast cancer dominates this section |
| ERBB2 | Similar to EPCAM | ✅ HER2 co-expressed in tumor zones — possible HER2+ component |
| KRT5 | Very low/sparse | ℹ️ Minimal basal-like component — this is likely ER+ luminal tumor |
| VIM | Diffuse, inverse of EPCAM | ✅ CAFs/stroma fills the spaces between tumor nests |
| COL1A1 | Strong lower region | ✅ Dense stromal collagen — cluster 0 (blue) = fibrous stroma |
| CD68 | Scattered spots | ✅ Macrophage infiltration — matches cluster 5 (brown) |
| CD3D | Diffuse low expression | ℹ️ T cells present but not highly infiltrated — immune-excluded? |
| MKI67 | Focal hot spots | ✅ Proliferating tumor cells — overlaps with EPCAM+ regions |
| PECAM1 | Very sparse | ✅ Rare endothelial spots — blood vessels (cluster 8 cyan) |
