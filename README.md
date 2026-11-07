# Human Dentate Gyrus Xenium Spatial Transcriptomics Analysis

## Overview

This pipeline analyzes human dentate gyrus (DG) granule cells using 10x Genomics Xenium spatial transcriptomics data (~366 gene panel). The dataset spans **10 DG FOVs** and **25 cortex FOVs** across **7 clinical conditions**: MAID, SEPSIS, HIBI, TBI, CARDIAC, METASTASIS, and EPILEPSY. The primary focus is on **granule cell (GC) subtypes** isolated from the combined dataset (~115K cells after QC).

**Lab:** Cembrowski Lab  
**Analyst:** Hope Luo  
**Platform:** 10x Genomics Xenium  
**Reference for annotation:** Thompson et al. 2025 human hippocampus snRNA-seq (humanHippocampus2024 package, ExperimentHub EH9606)

---

## File Structure

### Input Files

| File | Description |
|------|-------------|
| `DG_merged_raw_obj.rds` | Merged Seurat object with 10 DG FOVs |
| `DG_combined_raw_obj.rds` | Combined object: 10 DG + 25 cortex FOVs (~2.26M cells) |
| `gc_all_combined_updated.rds` | Isolated granule cells after QC (~115K cells) |
| `SingleR_output/DG_merged_clean_SingleR_with_all_labels.rds` | DG object with SingleR annotations at all levels |

### Output Directories

| Directory | Contents |
|-----------|----------|
| `SingleR_output/` | Per-FOV SingleR result .rds files |
| `GC_DEG_output/` | DEG result CSVs (pseudobulk and single-cell) |

### Script Files

| File | Description |
|------|-------------|
| `DG_xenium_analysis.R` | Main reorganized analysis script (this README documents it) |

---

## Objects and Their Roles

### `merged_obj` → `merged_obj_clean`
- **Scope:** 10 DG FOVs only
- **Used in:** Sections 1–3
- **Purpose:** Initial QC, normalization, clustering, and SingleR cell type annotation
- **Key operations:** SCTransform, PCA, UMAP, FindAllMarkers, SingleR at mid/fine/superfine levels

### `combined_obj`
- **Scope:** 10 DG FOVs + 25 cortex FOVs (~2.26M cells)
- **Used in:** Sections 4–5
- **Purpose:** Merge DG and cortex datasets, unify cell type labels into `combined_celltype`

### `gc_all` (primary working object)
- **Scope:** ~115K isolated granule cells from combined dataset
- **Used in:** Sections 6–16
- **Purpose:** All downstream GC analysis — clustering, excitability, maturity, POSTN biology, pseudotime, DEG
- **Key metadata columns:**

| Column | Description | Values |
|--------|-------------|--------|
| `condition` | Clinical condition | MAID, SEPSIS, HIBI, TBI, CARDIAC, METASTASIS, EPILEPSY |
| `sex` | Patient sex | M, F |
| `patient_age` | Patient age in years | 24–83 |
| `age_group` | Binned age | 20-30, 31-40, 41-50, 51-60, 61-70, 70+ |
| `ILAE_score` | Epilepsy severity (epilepsy only) | 0, 1, 2, 3; "PM" for postmortem |
| `tissue_type` | Tissue origin | surgical (epilepsy), postmortem (all others) |
| `group_label` | Condition + tissue type | e.g., "SEPSIS (postmortem)" |
| `seurat_clusters` | Unsupervised clusters (dims=1:25, res=0.2) | 0, 1, 2, ... |

---

## Script Sections

### Section 0: Environment Setup
**Object:** None (setup only)

Loads all required R packages and sets global options. Run this section every time you start a new R session. Sets file paths, creates output directories, and defines standard color palettes.

**Key settings:**
- `future.globals.maxSize = 250000 * 1024^2` — prevents memory errors in Seurat parallelization
- `set.seed(42)` — reproducibility

### Section 0B: Custom Helper Functions
**Object:** None (function definitions)

| Function | Purpose |
|----------|---------|
| `fix_sct_models()` | Fixes SCT model UMI assay references after merging (forces all to "Xenium") |
| `spatial_feature_plot()` | Plots spatial expression with high-expression cells rendered on top |
| `kw_within_group()` | Kruskal-Wallis test for within-condition sample heterogeneity |
| `plot_kw_gene()` | Violin + boxplot with Dunn's post-hoc significance brackets |
| `get_gene_attributes()` | Combines SCT feature attributes from merged object layers |
| `residualVarPlot()` | Residual variance vs. gene mean for SCTransform QC |
| `plot_pb_box()` | Pseudobulk boxplot for condition comparisons |

---

### PART 1: DG-Only Analysis

#### Section 1: Load DG Merged Object
**Object:** `merged_obj`

Loads the pre-merged object containing 10 DG FOVs (M38, M44, M56, M70, P28, P51, P60, P71, M25, M41).

#### Section 2: QC, Normalization, Dimensionality Reduction, Clustering
**Object:** `merged_obj` → `merged_obj_clean`

Pipeline: QC filter (nFeature > 5, nCount > 0) → SCTransform → PCA (50 PCs) → FindNeighbors (dims 1:20) → FindClusters (res 0.2) → UMAP → subset clusters 0–10 → re-process clean subset → FindAllMarkers.

**Known fix applied:** After merging, SCT models reference wrong assay names. `fix_sct_models()` forces all to "Xenium" before `PrepSCTFindMarkers`.

#### Section 3: SingleR Cell Type Annotation
**Object:** `merged_obj_clean`

Runs SingleR per FOV (to avoid memory issues) at two levels:
- **fine.cell.class:** GC, CA1, CA3, Sub, Astro, Oligo, etc.
- **superfine.cell.class:** GC.1–GC.5, CA1, CA3.1–CA3.2, etc.

Saves per-FOV results and adds four metadata columns: `SR_fine`, `SR_fine_pruned`, `SR_superfine`, `SR_superfine_pruned`.

---

### PART 2: Combine DG + Cortex

#### Section 4: Merge DG and Cortex Datasets
**Object:** `combined_obj`

Loads the pre-merged combined object (~2.26M cells). The cortex dataset (25 FOVs from LK's epilepsy study) has its own SingleR-based labels (`cellsubtype`, `cellsubtype_fine`).

#### Section 5: Unify Cell Type Labels
**Object:** `combined_obj`

Merges DG labels (`SR_superfine_pruned`) and cortex labels (`cellsubtype`) into a single `combined_celltype` column. Cortex labels take priority where both exist.

---

### PART 3: Granule Cell Analysis

#### Section 6: Isolate Granule Cells
**Object:** `gc_all` (created from `combined_obj`)

Extracts all cells labeled as GC.1–GC.5 (DG) or ExN_GC (cortex). Collapses all subtypes to "GC" for unbiased re-clustering.

#### Section 7: Add Metadata
**Object:** `gc_all`

Assigns condition, sex, patient age, age group, ILAE score, tissue type, and group label based on `orig.ident` (sample ID).

**Sample-to-condition mapping:**
- MAID: M70, P28
- SEPSIS: P51, P60, P71, M56
- HIBI: M44, M25, M41
- TBI: M38
- CARDIAC: L10
- METASTASIS: L5
- EPILEPSY: all remaining (25 cortex FOVs)

#### Section 8: Remove Low-Quality Samples
**Object:** `gc_all`

Removes 4 samples with unacceptably low median nCount:
- M44 (HIBI): nCount = 15
- M25 (HIBI): nCount = 28
- L5 (METASTASIS): nCount = 31
- P60 (SEPSIS): previously removed

M41 (HIBI, nCount = 106) is retained as borderline acceptable.

#### Section 9: SCTransform, PCA, Clustering, UMAP
**Object:** `gc_all`

Full re-processing pipeline on isolated GCs:
- SCTransform normalization
- PCA (50 PCs computed)
- **dims = 1:25** (captures 64.1% variance; PCs checked for nCount confound — only PC2 rho=-0.25 and PC6 rho=0.19 are moderate)
- **resolution = 0.2** for clustering
- UMAP for visualization
- JoinLayers + NormalizeData on Xenium assay (required before FetchData/GetAssayData)
- FindAllMarkers for cluster markers

**Critical technical notes:**
- After SCTransform on a merged object, layers are split by sample. `JoinLayers` is required before any expression data access.
- `fix_sct_models()` must be called before `PrepSCTFindMarkers`.
- `FindClusters` uses `graph.name = "SCT_snn"` (check with `names(gc_all@graphs)`).

#### Section 10: Basic UMAP Visualizations
**Object:** `gc_all`

Standard UMAPs colored by cluster, FOV, condition, sex, age group, ILAE score, tissue type. Includes `FeaturePlot` for PROX1 and `spatial_feature_plot` examples.

#### Section 11: Excitability Analysis
**Object:** `gc_all`

Computes module scores for hypo- and hyper-excitability gene sets:
- **Hypo:** SLC12A5, CALB1, GABRA5, GRIN2A, GABRB3, SEMA5A
- **Hyper:** HOMER1, SNAP25, KCNMA1, DLG1, DLG2, GABRA2

`AddModuleScore` requires `ctrl = 5, nbin = 10` (Xenium panel is too small for defaults).

**Major caveat:** Tissue type (surgical vs. postmortem) is a systematic confound for excitability scores. Within-postmortem comparisons are valid; cross-tissue comparisons are not.

#### Section 12: nCount Confound Analysis
**Object:** `gc_all`

Quantifies the correlation between gene expression and sequencing depth (nCount).

**Key findings:**
| Gene | Spearman rho with nCount | Interpretation |
|------|--------------------------|----------------|
| PROX1 | -0.92 | Heavily confounded — DO NOT trust raw expression |
| GRIA1 | -0.97 | Heavily confounded |
| SLC17A7 | -0.95 | Heavily confounded |
| SNAP25 | -0.79 | Confounded |
| GRIN3A | -0.77 | Confounded — needs correction for maturity score |
| CALB1 | +0.29 | NOT confounded — safe for cross-FOV comparisons |
| POSTN | +0.08 | NOT confounded — safe for cross-FOV comparisons |

#### Section 13: CALB1 Maturity Axis and HOMER1 Anti-correlation
**Object:** `gc_all`

The UMAP gradient represents a genuine biological maturity continuum (not a technical artifact):
- **CALB1** marks mature GCs (not confounded by nCount)
- **GRIN3A** marks immature GCs (confounded — needs correction)
- **Maturity score** = CALB1 − GRIN3A, corrected for nCount via linear regression residuals: `residuals(lm(maturity_raw ~ nCount_Xenium))`

HOMER1/CALB1 anti-correlation (rho = -0.20) reflects an activity-state axis independent of maturity.

#### Section 14: POSTN Biology (Cluster 2)
**Object:** `gc_all`, `gc_epi`

POSTN (Periostin) is exclusively expressed in cluster 2 and is NOT confounded by nCount. Marks the most structurally stable, ECM-anchored mature GCs.

**Key finding:** POSTN+ proportion changes with ILAE severity:
- ILAE 0: ~11.3%
- ILAE 2: ~12.3% (compensatory ECM reinforcement)
- ILAE 3: ~7.5% (decompensation at severe sclerosis)

#### Section 15: Slingshot Pseudotime Trajectory
**Object:** `gc_all` (downsampled to 20K cells)

Infers a pseudotime trajectory using Slingshot on PCA embeddings. Full dataset (115K cells) is too large for Slingshot, so cells are downsampled to 20K.

**Critical fix:** SCE conversion does not automatically carry PCA embeddings from Seurat. Must manually set: `reducedDim(sce, "PCA") <- Embeddings(gc_sub, "pca")`

Literature-defined maturity stages:
1. **Immature:** DCX, STMN1, DPYSL3
2. **Early:** CALB2, GRIN3A
3. **Mid:** GRIA1, SLC17A7
4. **Mature:** CALB1, PROX1

#### Section 16: Differential Gene Expression Analysis
**Object:** `gc_all`, `pseudo_gc`

**16A — Single-cell DEG:** FindAllMarkers per cluster (already computed in Section 9F). Visualized with DotPlot and DoHeatmap.

**16B — Pseudobulk DEG (DESeq2):** Aggregates expression by condition × sample × cluster using `AggregateExpression()`, then runs pairwise DESeq2 comparisons. Most reliable comparison: SEPSIS vs HIBI (sufficient samples per group).

**16C — Within-group Kruskal-Wallis:** Tests whether individual samples within a condition differ from each other per gene. Useful for detecting within-group heterogeneity before interpreting between-group differences.

---

## Known Issues and Fixes

| Error | Cause | Fix |
|-------|-------|-----|
| `PrepSCTFindMarkers: Multiple UMI assays` | Merged SCT models reference wrong assay | `fix_sct_models()` — loop through models, set umi.assay to "Xenium" |
| `quantile.default` error in FeaturePlot | Split layers → NAs for cells from other samples | `JoinLayers()` before FeaturePlot |
| `future.globals.maxSize exceeded` | Seurat parallelization memory limit | `options(future.globals.maxSize = 250000 * 1024^2)` |
| `AddModuleScore: sample larger than population` | Xenium panel too small for default ctrl/nbin | `ctrl = 5, nbin = 10` |
| `FindClusters: graph.name not present` | Wrong graph name after SCTransform | Check `names(obj@graphs)`, use `graph.name = "SCT_snn"` |
| `GetAssayData: doesn't work for multiple layers` | Split layers in merged object | `JoinLayers()` first |
| `slingPseudotime: unable to find inherited method for NULL` | SCE missing PCA embeddings | `reducedDim(sce, "PCA") <- Embeddings(obj, "pca")` |
| `lm(): object 'nCount_Xenium' not found` | Formula evaluated without data argument | Pass `data = gc_all@meta.data` |
| `stat_pvalue_manual: aesthetics length mismatch` | ggpubr inherits aesthetics from base plot | Add `inherit.aes = FALSE` |
| Slingshot runs forever | 115K cells too large | Downsample to 20K cells |
| FeaturePlot shows dark grey cells | Split layers → NA expression rendered as dark grey | `JoinLayers()` before plotting |

---

## Analytical Caveats

1. **nCount confound:** Most genes are strongly anti-correlated with sequencing depth. Only CALB1 and POSTN are reliable for raw cross-FOV comparisons. Use nCount-corrected residuals for confounded genes.

2. **Tissue type confound:** Surgical (epilepsy) vs. postmortem (all others) tissue systematically affects excitability scores. Only within-tissue-type comparisons are valid.

3. **Low-quality samples:** M44, M25, and L5 were removed due to very low nCount. M41 is borderline but retained.

4. **Sample size limitations:** Some conditions have only 1–2 samples (TBI: M38 only; CARDIAC: L10 only), limiting statistical power for between-condition comparisons.

5. **Pseudobulk DEG:** Requires ≥3 biological replicates per group per cluster. Many condition pairs lack sufficient samples.
