################################################################################
#                                                                              #
#   HUMAN DENTATE GYRUS (DG) XENIUM SPATIAL TRANSCRIPTOMICS ANALYSIS           #
#   Cembrowski Lab — Hope Luo                                                  #
#                                                                              #
#   Platform: 10x Genomics Xenium (~366 gene panel)                            #
#   Dataset:  10 DG FOVs + 25 cortex FOVs across 7 conditions                 #
#             (MAID, SEPSIS, HIBI, TBI, CARDIAC, METASTASIS, EPILEPSY)         #
#                                                                              #
#   Primary focus: Granule cell (GC) subtypes from the combined dataset        #
#                                                                              #
#   STRUCTURE:                                                                 #
#     Section 0:  Environment setup (packages, paths, options)                 #
#     Section 0B: Custom helper functions                                      #
#     PART 1 — DG-only analysis (Sections 1–3)                                 #
#       Section 1: Load DG merged object                                       #
#       Section 2: QC, SCTransform, PCA, clustering, UMAP                      #
#       Section 3: SingleR annotation (mid, fine, superfine)                   #
#     PART 2 — Combine DG + cortex (Sections 4–5)                              #
#       Section 4: Merge DG and cortex datasets                                #
#       Section 5: Unify cell type labels                                      #
#     PART 3 — Granule cell analysis (Sections 6–16)                           #
#       Section 6:  Isolate granule cells (gc_all)                             #
#       Section 7:  Add metadata (condition, sex, age, ILAE, tissue type)      #
#       Section 8:  Remove low-quality samples                                 #
#       Section 9:  SCTransform, PCA, clustering (dims=1:25, res=0.2), UMAP   #
#       Section 10: Basic UMAP visualizations                                  #
#       Section 11: Excitability analysis (hypo/hyper module scores)           #
#       Section 12: nCount confound analysis                                   #
#       Section 13: CALB1 maturity axis and HOMER1 anti-correlation            #
#       Section 14: POSTN biology (cluster 2)                                  #
#       Section 15: Slingshot pseudotime trajectory                            #
#       Section 16: DEG analysis (single-cell and pseudobulk)                  #
#                                                                              #
################################################################################


###############################################################################
# SECTION 0: ENVIRONMENT SETUP                                                #
#                                                                             #
# Load all required packages and set global options.                          #
# Run this section every time you start a new R session.                      #
###############################################################################

suppressPackageStartupMessages({
  # Core Seurat ecosystem
  library(Seurat)              # Seurat v5 — main scRNA-seq framework
  library(SeuratObject)        # Object classes

  # Data wrangling
  library(dplyr)               # Data manipulation (filter, mutate, summarise)
  library(tidyr)               # Reshape data (pivot_longer, pivot_wider)
  library(tibble)              # Enhanced data frames
  library(stringr)             # String manipulation
  library(readxl)              # Read Excel files

  # Plotting
  library(ggplot2)             # Grammar of graphics
  library(patchwork)           # Combine multiple ggplots (+ and / operators)
  library(scattermore)         # Fast scatter plots for large datasets
  library(scales)              # Axis scaling helpers
  library(RColorBrewer)        # Color palettes
  library(ggrepel)             # Non-overlapping text labels
  library(ggpubr)              # Publication-ready plots, stat_pvalue_manual()

  # Cell type annotation
  library(SingleR)             # Reference-based cell type annotation
  library(SingleCellExperiment) # SCE class for SingleR
  library(SummarizedExperiment) # SE class
  library(ExperimentHub)       # Access reference datasets
  library(humanHippocampus2024) # Thompson et al. 2025 hippocampus reference

  # Differential expression
  library(DESeq2)              # Pseudobulk DEG analysis

  # Statistical testing
  library(FSA)                 # Dunn's test (dunnTest)
  library(rstatix)             # Pipe-friendly stats
  library(broom)               # Tidy model outputs

  # Trajectory analysis
  library(slingshot)           # Pseudotime trajectory inference
})

# --- Global options ---
options(future.globals.maxSize = 250000 * 1024^2)  # Prevent memory errors in Seurat
set.seed(42)                                        # Reproducibility

# --- File paths ---
# Base directory for all DG analysis files
base_dir     <- "X:/CembrowskiLab/Hope Luo/DG analysis"
singleR_dir  <- file.path(base_dir, "SingleR_output")
deg_dir      <- file.path(base_dir, "GC_DEG_output")
xenium_dir   <- file.path(base_dir, "XeniumOutput")

# Key RDS file paths
merged_raw_path   <- file.path(base_dir, "DG_merged_raw_obj.rds")
combined_raw_path <- file.path(base_dir, "DG_combined_raw_obj.rds")
gc_all_path       <- file.path(base_dir, "gc_all_combined_updated.rds")
singleR_obj_path  <- file.path(singleR_dir, "DG_merged_clean_SingleR_with_all_labels.rds")

# Create output directories if they don't exist
for (d in c(singleR_dir, deg_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Standard color palettes ---
# GC subtype colors (used consistently throughout)
gc_colors <- c(
  "GC.1" = "#E64B35",
  "GC.2" = "#4DBBD5",
  "GC.3" = "#00A087",
  "GC.4" = "#3C5488",
  "GC.5" = "#F39B7F"
)

# Condition colors
condition_colors <- c(
  "MAID"       = "#E64B35",
  "SEPSIS"     = "#4DBBD5",
  "HIBI"       = "#00A087",
  "TBI"        = "#3C5488",
  "CARDIAC"    = "#F39B7F",
  "METASTASIS" = "#91D1C2",
  "EPILEPSY"   = "#8491B4"
)


###############################################################################
# SECTION 0B: CUSTOM HELPER FUNCTIONS                                         #
#                                                                             #
# Reusable functions for common operations throughout the analysis.           #
###############################################################################

#' Fix SCT model UMI assay references
#'
#' After merging Seurat objects, SCT models may reference wrong assay names.
#' This forces all models to reference "Xenium", fixing PrepSCTFindMarkers errors.
#'
#' @param obj A Seurat object with SCT assay
#' @return The same object with fixed SCT models
fix_sct_models <- function(obj) {
  for (i in seq_along(obj[["SCT"]]@SCTModel.list)) {
    slot(obj[["SCT"]]@SCTModel.list[[i]], "umi.assay") <- "Xenium"
  }
  return(obj)
}


#' Spatial feature plot with high-expression cells on top
#'
#' Standard Seurat spatial plots render cells in random order, so high-expressing
#' cells can be hidden behind low-expressing ones. This function sorts cells by
#' expression (ascending) so that high-expression cells are plotted last and
#' appear on top of the plot.
#'
#' @param obj        Seurat object
#' @param feature    Gene name to plot
#' @param fov        FOV name (must match an existing FOV in the object)
#' @param pt_size    Point size (default 0.5)
#' @param low_color  Color for low expression (default light grey)
#' @param high_color Color for high expression (default dark red)
#' @param title      Optional plot title (defaults to "feature — fov")
#' @return A ggplot object
spatial_feature_plot <- function(obj, feature, fov,
                                 pt_size = 0.5,
                                 low_color = "#F0F0F0",
                                 high_color = "#B2182B",
                                 title = NULL) {

  # Get tissue coordinates for the specified FOV
  coords <- GetTissueCoordinates(obj, image = fov)

  # Get expression data for the feature
  # Ensure layers are joined and data is normalized first
  expr_values <- FetchData(obj, vars = feature, cells = rownames(coords))

  # Combine coordinates and expression
  plot_df <- cbind(coords[, c("x", "y")], expr = expr_values[[feature]])

  # Sort by expression ASCENDING so high-expression cells plot on top
  plot_df <- plot_df[order(plot_df$expr), ]

  # Build the plot
  if (is.null(title)) title <- paste0(feature, " — ", fov)

  ggplot(plot_df, aes(x = x, y = y, color = expr)) +
    geom_point(size = pt_size, shape = 16) +
    scale_color_gradient(low = low_color, high = high_color, name = feature) +
    coord_fixed() +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "right"
    ) +
    ggtitle(title)
}


#' Kruskal-Wallis test within a condition group
#'
#' Tests whether individual samples within a condition differ from each other
#' for each gene. Useful for detecting within-group heterogeneity.
#'
#' @param obj            Seurat object
#' @param condition_name The condition to subset to
#' @return Data frame with per-gene KW test results, or NULL if <3 samples
kw_within_group <- function(obj, condition_name) {
  sub_obj <- subset(obj, subset = condition == condition_name)

  n_samples <- length(unique(sub_obj$orig.ident))
  if (n_samples < 3) {
    message("Skipping ", condition_name, " — only ", n_samples, " sample(s)")
    return(NULL)
  }

  expr_mat <- GetAssayData(sub_obj, layer = "data", assay = "Xenium")
  sample_ids <- sub_obj$orig.ident

  results <- lapply(rownames(expr_mat), function(gene) {
    values <- expr_mat[gene, ]
    if (sd(values) == 0) return(NULL)

    kw <- kruskal.test(values ~ sample_ids)
    data.frame(
      gene = gene,
      condition = condition_name,
      kw_statistic = kw$statistic,
      p_value = kw$p.value,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, results)
  df$p_adj <- p.adjust(df$p_value, method = "BH")
  df <- df %>% arrange(p_adj)
  return(df)
}


#' Plot Kruskal-Wallis result with Dunn's post-hoc brackets
#'
#' Creates a violin + boxplot for one gene within one condition group,
#' with significance brackets from Dunn's pairwise comparisons.
#'
#' @param obj            Seurat object
#' @param gene           Gene name
#' @param condition_name Condition to subset to
#' @param kw_results     Data frame of KW results (for subtitle p-value)
#' @return A ggplot object
plot_kw_gene <- function(obj, gene, condition_name, kw_results = NULL) {
  sub_obj <- subset(obj, subset = condition == condition_name)

  df <- data.frame(
    expr = GetAssayData(sub_obj, layer = "data", assay = "Xenium")[gene, ],
    sample = sub_obj$orig.ident
  )

  # Run Dunn's post-hoc test
  dunn_res <- dunnTest(expr ~ sample, data = df, method = "bh")
  dunn_df <- dunn_res$res
  sig_pairs <- dunn_df %>% filter(P.adj < 0.05)

  # Base plot
  p <- ggplot(df, aes(x = sample, y = expr, fill = sample)) +
    geom_violin(scale = "width", alpha = 0.7, width = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
    theme_classic() +
    labs(title = paste0(gene, " — ", condition_name), x = NULL, y = "Normalized Expression") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  # Add significance brackets
  if (nrow(sig_pairs) > 0) {
    sig_pairs <- sig_pairs %>%
      mutate(
        group1 = trimws(sapply(strsplit(Comparison, " - "), `[`, 1)),
        group2 = trimws(sapply(strsplit(Comparison, " - "), `[`, 2)),
        p.adj.signif = case_when(
          P.adj < 0.001 ~ "***",
          P.adj < 0.01  ~ "**",
          P.adj < 0.05  ~ "*",
          TRUE          ~ "ns"
        )
      )

    y_max <- max(df$expr, na.rm = TRUE)
    y_step <- if (y_max == 0) 0.2 else y_max * 0.12
    sig_pairs$y.position <- y_max + (seq_len(nrow(sig_pairs)) * y_step)

    p <- p + stat_pvalue_manual(
      sig_pairs, label = "p.adj.signif", y.position = "y.position",
      step.increase = 0.08, tip.length = 0.01, size = 3.5, inherit.aes = FALSE
    ) +
      ylim(min(df$expr, na.rm = TRUE), max(sig_pairs$y.position) + y_step)
  }

  # Add KW p-value as subtitle
  if (!is.null(kw_results)) {
    kw_p <- kw_results %>% filter(gene == !!gene, condition == condition_name) %>% pull(p_adj)
    if (length(kw_p) > 0) {
      p <- p + labs(subtitle = paste0("KW adj. p = ", formatC(kw_p[1], format = "e", digits = 2)))
    }
  }

  return(p)
}


#' Extract and combine SCT feature attributes from a merged object
#'
#' For merged Seurat objects, SCTResults returns a list of data frames
#' (one per layer/sample). This function combines them into a single data frame.
#'
#' @param obj   Seurat object
#' @param assay Assay name (default "SCT")
#' @return A combined data frame of gene attributes across all models
get_gene_attributes <- function(obj, assay = "SCT") {
  res <- SCTResults(obj, slot = "feature.attributes", assay = assay)

  if (is.data.frame(res)) {
    res$gene <- rownames(res)
    res$model <- "Model_1"
    return(res)
  }

  combined_df <- do.call(rbind, lapply(names(res), function(model_name) {
    df <- res[[model_name]]
    df$gene <- rownames(df)
    df$model <- model_name
    return(df)
  }))

  return(combined_df)
}


#' Residual variance plot for SCTransform QC
#'
#' Plots residual variance vs gene mean expression to assess normalization quality.
#' Genes above the dashed line at 1 have more variance than expected.
#'
#' @param input      Seurat object or data frame of gene attributes
#' @param xaxis      Column name for x-axis (default "gmean")
#' @param max_resvar Cap residual variance at this value (default 100)
#' @param ntop       Number of top-variance genes to highlight (default 20)
#' @param annotate   Add gene name labels to top genes (default FALSE)
#' @param pt_size    Point size (default 1.1)
#' @return A ggplot object
residualVarPlot <- function(input, xaxis = "gmean", max_resvar = 100,
                            ntop = 20, annotate = FALSE, pt_size = 1.1) {

  if (inherits(input, "Seurat")) {
    model_name <- names(input[["SCT"]]@SCTModel.list)[1]
    gene_var <- input[["SCT"]]@SCTModel.list[[model_name]]@feature.attributes
  } else {
    gene_var <- input
  }

  gene_var$gene <- rownames(gene_var)
  topn <- subset(gene_var, rank(-gene_var[, "residual_variance"]) <= ntop)$gene
  gene_var[gene_var$residual_variance > max_resvar, "residual_variance"] <- max_resvar

  p <- ggplot(gene_var, aes(.data[[xaxis]], .data[["residual_variance"]])) +
    geom_scattermore(pointsize = pt_size, shape = 16, alpha = 0.5, color = "#43a2ca") +
    geom_scattermore(data = subset(gene_var, gene %in% topn),
                     pointsize = pt_size, shape = 16, alpha = 1.0, color = "deeppink") +
    geom_hline(yintercept = 1, color = "#4daf4a", linewidth = 0.9, linetype = "dashed") +
    geom_smooth(method = "loess", span = 0.1, linewidth = 0.9, formula = y ~ x, color = "#e41a1c") +
    scale_y_continuous(trans = "sqrt", breaks = c(0, 1, 10, 25, 50, 100, 150),
                       limits = c(0, max_resvar + 1)) +
    scale_x_continuous(trans = "log10", breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                       labels = MASS::rational) +
    xlab("Gene mean expression") + ylab("Residual variance")

  if (annotate) {
    p <- p + geom_text_repel(
      data = subset(gene_var, gene %in% topn), aes(label = gene),
      color = "gray25", size = 1.8, direction = "x",
      vjust = 0.5, hjust = 0.5, segment.size = 0.2, segment.alpha = 0.2
    )
  }

  cat("Top", ntop, "genes by residual variance:\n")
  print(topn)

  return(p)
}


#' Pseudobulk boxplot for condition comparisons
#'
#' Boxplot of aggregated expression for one gene in one cell type,
#' comparing two condition groups.
#'
#' @param pseudo    Pseudobulk Seurat object (from AggregateExpression)
#' @param gene      Gene name
#' @param celltype  Cell type identity label
#' @param group1    First condition group name
#' @param group2    Second condition group name
#' @return A ggplot object
plot_pb_box <- function(pseudo, gene, celltype, group1, group2) {
  ids <- c(paste(celltype, group1, sep = "_"),
           paste(celltype, group2, sep = "_"))

  df <- FetchData(pseudo, vars = c("subtype_condition", "condition", gene))
  colnames(df)[ncol(df)] <- "expr"
  df <- df %>% filter(subtype_condition %in% ids)

  ggplot(df, aes(x = subtype_condition, y = expr, color = condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.35) +
    geom_point(size = 2.8, position = position_jitter(width = 0.08, height = 0)) +
    theme_classic() +
    labs(title = paste0(gene, " in ", celltype), x = NULL, y = "Aggregated expression")
}


################################################################################
#                                                                              #
#  ██████   █████  ██████  ████████     ██                                     #
#  ██   ██ ██   ██ ██   ██    ██       ███                                     #
#  ██████  ███████ ██████     ██        ██                                     #
#  ██      ██   ██ ██   ██    ██        ██                                     #
#  ██      ██   ██ ██   ██    ██        ██                                     #
#                                                                              #
#  DG-ONLY ANALYSIS (Sections 1–3)                                             #
#  Object: merged_obj → merged_obj_clean                                       #
################################################################################


###############################################################################
# SECTION 1: LOAD DG MERGED OBJECT                                            #
#                                                                             #
# The merged object contains 10 DG FOVs merged from individual Xenium runs.  #
# FOVs: M38, M44, M56, M70, P28, P51, P60, P71, M25, M41                    #
#                                                                             #
# NOTE: Individual Xenium runs were loaded using LoadXeniumV3() and saved     #
# as individual .rds files, then merged with Seurat::merge().                 #
###############################################################################

merged_obj <- readRDS(merged_raw_path)


###############################################################################
# SECTION 2: QC, NORMALIZATION, DIMENSIONALITY REDUCTION, CLUSTERING          #
#                                                                             #
# Pipeline: QC filter → SCTransform → PCA → FindNeighbors → FindClusters →   #
#           RunUMAP → subset to clusters 0-10 → re-process clean subset       #
###############################################################################

# --- 2A: Quality control ---
# Visualize QC metrics across FOVs
# nFeature = number of unique genes detected per cell
# nCount = total transcript counts per cell
VlnPlot(merged_obj, features = "nFeature_Xenium", group.by = "orig.ident", pt.size = 0)
VlnPlot(merged_obj, features = "nCount_Xenium", group.by = "orig.ident", pt.size = 0)

# Filter: keep cells with >5 genes detected and >0 total counts
merged_obj <- subset(merged_obj, subset = nFeature_Xenium > 5 & nCount_Xenium > 0)

# --- 2B: SCTransform normalization ---
# SCTransform performs normalization, variance stabilization, and feature selection
# in one step. Preferred over log-normalization for spatial data.
merged_obj <- SCTransform(merged_obj, assay = "Xenium")

# Visualize normalization quality
gene_attr_merged <- get_gene_attributes(merged_obj, assay = "SCT")
residualVarPlot(gene_attr_merged, max_resvar = 15, pt_size = 3, annotate = TRUE)

# --- 2C: PCA ---
# Compute 50 PCs from variable features identified by SCTransform
merged_obj <- RunPCA(merged_obj, features = VariableFeatures(merged_obj),
                     npcs = 50, verbose = FALSE)

# Elbow plot to choose number of PCs for downstream analysis
ElbowPlot(merged_obj, ndims = 50)

# PCA heatmap for top features driving first 8 PCs
DimHeatmap(merged_obj, dims = 1:8, cells = 50, balanced = TRUE)

# --- 2D: Clustering ---
# Using dims 1:20 and resolution 0.2 for initial broad clustering
merged_obj <- FindNeighbors(merged_obj, dims = 1:20)
merged_obj <- FindClusters(merged_obj, resolution = 0.2)

# --- 2E: UMAP ---
merged_obj <- RunUMAP(merged_obj, dims = 1:20)

# --- 2F: Remove noise clusters ---
# Clusters beyond 10 typically represent low-quality cells or doublets
# Inspect QC metrics per cluster before subsetting
VlnPlot(merged_obj, features = "nCount_Xenium", group.by = "seurat_clusters", pt.size = 0.1)
VlnPlot(merged_obj, features = "nFeature_Xenium", group.by = "seurat_clusters", pt.size = 0.1)

# Keep only clusters 0–10
merged_obj_clean <- subset(merged_obj, idents = as.character(0:10))
table(Idents(merged_obj_clean))

# Re-run PCA and UMAP on clean subset
merged_obj_clean <- RunPCA(merged_obj_clean, verbose = FALSE)
merged_obj_clean <- RunUMAP(merged_obj_clean, dims = 1:20)

# --- 2G: Visualize clean UMAP ---
DimPlot(merged_obj_clean, label = TRUE) + coord_fixed()

# Spatial plot for one FOV
ImageDimPlot(merged_obj_clean, fov = "M38", size = 0.5, border.size = NA)

# --- 2H: Find cluster markers ---
# Fix SCT models (required after merge) and find markers
merged_obj_clean <- fix_sct_models(merged_obj_clean)
merged_obj_clean <- PrepSCTFindMarkers(merged_obj_clean)

all_markers <- FindAllMarkers(merged_obj_clean, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25)

top25 <- all_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 25) %>%
  ungroup()

top25_list <- split(top25$gene, top25$cluster)


###############################################################################
# SECTION 3: SINGLER CELL TYPE ANNOTATION                                     #
#                                                                             #
# Uses the Thompson et al. 2025 human hippocampus snRNA-seq reference         #
# (humanHippocampus2024 package, ExperimentHub EH9606).                       #
#                                                                             #
# Three annotation levels:                                                    #
#   mid.cell.class   → broad types (ExcN, InhN, Astro, Oligo, etc.)          #
#   fine.cell.class  → finer types (GC, CA1, CA3, Sub, etc.)                  #
#   superfine.cell.class → subtypes (GC.1–GC.5, CA1, CA3.1–CA3.2, etc.)     #
#                                                                             #
# SingleR runs per FOV to avoid memory issues with large merged objects.      #
###############################################################################

# --- 3A: Load reference ---
eh <- ExperimentHub()
sce_hpc <- eh[["EH9606"]]

ref_labels_fine      <- sce_hpc$fine.cell.class
ref_labels_superfine <- sce_hpc$superfine.cell.class
ref_expr             <- logcounts(sce_hpc)

message("Reference loaded: ", ncol(sce_hpc), " cells, ",
        length(unique(ref_labels_superfine)), " superfine cell types")

# --- 3B: Normalize query data ---
merged_obj_clean <- NormalizeData(merged_obj_clean, assay = "Xenium", verbose = FALSE)

# --- 3C: Run SingleR per FOV at fine and superfine levels ---
sample_names <- unique(merged_obj_clean$orig.ident)

# Pre-allocate label vectors
all_fine              <- setNames(rep(NA, ncol(merged_obj_clean)), colnames(merged_obj_clean))
all_fine_pruned       <- all_fine
all_superfine         <- all_fine
all_superfine_pruned  <- all_fine

for (sample_name in sample_names) {
  message("=== Running SingleR for FOV: ", sample_name, " ===")

  sample_obj <- subset(merged_obj_clean, subset = orig.ident == sample_name)
  query_expr <- GetAssayData(sample_obj, layer = "data", assay = "Xenium")

  common_genes <- intersect(rownames(ref_expr), rownames(query_expr))
  message("  Overlapping genes: ", length(common_genes))

  if (length(common_genes) < 10) {
    warning("  Very few overlapping genes for ", sample_name, " — skipping")
    next
  }

  ref_sub   <- ref_expr[common_genes, , drop = FALSE]
  query_sub <- query_expr[common_genes, , drop = FALSE]

  # Fine level (GC, CA1, Astro, Oligo, etc.)
  sr_fine <- SingleR(test = as.matrix(query_sub), ref = as.matrix(ref_sub),
                     labels = ref_labels_fine)

  # Superfine level (GC.1–GC.5, etc.)
  sr_superfine <- SingleR(test = as.matrix(query_sub), ref = as.matrix(ref_sub),
                          labels = ref_labels_superfine)

  sample_cells <- colnames(sample_obj)
  all_fine[sample_cells]             <- sr_fine$labels
  all_fine_pruned[sample_cells]      <- sr_fine$pruned.labels
  all_superfine[sample_cells]        <- sr_superfine$labels
  all_superfine_pruned[sample_cells] <- sr_superfine$pruned.labels

  # Save per-sample results
  saveRDS(sr_fine, file.path(singleR_dir, paste0("SingleR_fine_", sample_name, "_results.rds")))
  saveRDS(sr_superfine, file.path(singleR_dir, paste0("SingleR_superfine_", sample_name, "_results.rds")))

  message("  Done.")
}

# Add all labels to metadata
merged_obj_clean$SR_fine              <- all_fine
merged_obj_clean$SR_fine_pruned       <- all_fine_pruned
merged_obj_clean$SR_superfine         <- all_superfine
merged_obj_clean$SR_superfine_pruned  <- all_superfine_pruned

# Save annotated object
saveRDS(merged_obj_clean, file.path(singleR_dir, "DG_merged_clean_SingleR_with_all_labels.rds"))

# --- 3D: Visualize annotations ---
# UMAP by superfine cell type
DimPlot(merged_obj_clean, reduction = "umap", group.by = "SR_superfine_pruned",
        label = TRUE, repel = TRUE) +
  coord_fixed() + ggtitle("SingleR Annotations (superfine.cell.class)")

# UMAP highlighting only GC subtypes (all others greyed out)
gc_subtypes <- c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5")
all_types <- unique(merged_obj_clean$SR_superfine_pruned)
highlight_cols <- setNames(rep("lightgrey", length(all_types)), all_types)
highlight_cols[gc_subtypes] <- gc_colors[gc_subtypes]

DimPlot(merged_obj_clean, reduction = "umap", group.by = "SR_superfine_pruned",
        cols = highlight_cols, label = FALSE, order = gc_subtypes) +
  coord_fixed() + ggtitle("GC Subtypes Highlighted")

# Proportional bar chart of cell types per FOV
ggplot(merged_obj_clean@meta.data, aes(x = orig.ident, fill = SR_superfine_pruned)) +
  geom_bar(position = "fill") +
  theme_classic() +
  labs(x = "FOV", y = "Proportion of Cells", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cell Type Composition by FOV")


################################################################################
#                                                                              #
#  ██████   █████  ██████  ████████    ██████                                  #
#  ██   ██ ██   ██ ██   ██    ██            ██                                 #
#  ██████  ███████ ██████     ██        █████                                  #
#  ██      ██   ██ ██   ██    ██       ██                                      #
#  ██      ██   ██ ██   ██    ██       ███████                                 #
#                                                                              #
#  COMBINE DG + CORTEX (Sections 4–5)                                          #
#  Object: combined_obj                                                        #
################################################################################


###############################################################################
# SECTION 4: MERGE DG AND CORTEX DATASETS                                     #
#                                                                             #
# The cortex dataset (25 FOVs from LK's epilepsy study) has its own           #
# SingleR-based labels (cellsubtype, cellsubtype_fine).                        #
# We merge it with the DG data to create a combined dataset of ~2.26M cells.  #
###############################################################################

# NOTE: The merging step was done previously and saved. Load directly:
combined_obj <- readRDS(combined_raw_path)

table(combined_obj$orig.ident)
message("Total cells in combined object: ", ncol(combined_obj))


###############################################################################
# SECTION 5: UNIFY CELL TYPE LABELS                                           #
#                                                                             #
# DG cells have SR_superfine_pruned labels; cortex cells have cellsubtype     #
# labels. Merge them into a single combined_celltype column.                  #
###############################################################################

has_cortex <- !is.na(combined_obj$cellsubtype)
has_dg     <- !is.na(combined_obj$SR_superfine_pruned)

message("Cortex-labeled cells: ", sum(has_cortex))
message("DG-labeled cells: ", sum(has_dg))
message("Overlap: ", sum(has_cortex & has_dg))

# Create unified label column
combined_obj$combined_celltype <- NA
combined_obj$combined_celltype[has_cortex] <- combined_obj$cellsubtype[has_cortex]

dg_cells <- !is.na(combined_obj$SR_superfine_pruned) & is.na(combined_obj$combined_celltype)
combined_obj$combined_celltype[dg_cells] <- combined_obj$SR_superfine_pruned[dg_cells]

message("Total labeled: ", sum(!is.na(combined_obj$combined_celltype)),
        " / ", ncol(combined_obj))


################################################################################
#                                                                              #
#  ██████   █████  ██████  ████████    ██████                                  #
#  ██   ██ ██   ██ ██   ██    ██            ██                                 #
#  ██████  ███████ ██████     ██        █████                                  #
#  ██      ██   ██ ██   ██    ██            ██                                 #
#  ██      ██   ██ ██   ██    ██       ██████                                  #
#                                                                              #
#  GRANULE CELL ANALYSIS (Sections 6–16)                                       #
#  Object: gc_all (~115K cells)                                                #
#                                                                              #
#  *** THIS IS THE PRIMARY WORKING SECTION ***                                 #
################################################################################


###############################################################################
# SECTION 6: ISOLATE GRANULE CELLS                                             #
#                                                                             #
# Extract all granule cells from the combined dataset. This includes:         #
#   - GC.1–GC.5 (DG samples, labeled by SingleR superfine)                    #
#   - ExN_GC (cortex samples, labeled by cellsubtype)                         #
# All are collapsed into a single "GC" identity for unified analysis.         #
###############################################################################

gc_labels <- c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5", "ExN_GC")

gc_all <- subset(combined_obj, subset = combined_celltype %in% gc_labels)

# Collapse all GC subtypes into a single label for unsupervised clustering
gc_all$combined_celltype <- "GC"

table(gc_all$orig.ident)
message("Total GC cells: ", ncol(gc_all))


###############################################################################
# SECTION 7: ADD METADATA                                                      #
#                                                                             #
# Add experimental metadata: condition, sex, age, ILAE score, tissue type,    #
# and group labels for downstream comparisons.                                #
###############################################################################

# --- 7A: Condition ---
gc_all$condition <- NA
gc_all$condition[gc_all$orig.ident %in% c("M70", "P28")]              <- "MAID"
gc_all$condition[gc_all$orig.ident %in% c("P51", "P60", "P71", "M56")] <- "SEPSIS"
gc_all$condition[gc_all$orig.ident %in% c("M44", "M25", "M41")]       <- "HIBI"
gc_all$condition[gc_all$orig.ident %in% c("M38")]                     <- "TBI"
gc_all$condition[gc_all$orig.ident %in% c("L10")]                     <- "CARDIAC"
gc_all$condition[gc_all$orig.ident %in% c("L5")]                      <- "METASTASIS"
gc_all$condition[is.na(gc_all$condition)]                              <- "EPILEPSY"

# --- 7B: Sex ---
gc_all$sex <- NA
gc_all$sex[gc_all$orig.ident %in% c("B1","D11","D8","I1","C1","C2","E008",
                                     "D12","D13","D14","D9","E015",
                                     "M56","M70","P28","P51","P71")] <- "M"
gc_all$sex[gc_all$orig.ident %in% c("F1","H1","J1","L20","L14","L15","L18",
                                     "L1","L21","L3","L2","L10","L5")]  <- "F"

# --- 7C: Patient age ---
age_map <- c(
  B1 = 45, D11 = 46, D8 = 41, F1 = 70, H1 = 74, I1 = 53, J1 = 83, L20 = 52,
  C1 = 51, C2 = 62, E008 = 44, L14 = 46, L15 = 48, L18 = 50,
  D12 = 24, D13 = 31, L1 = 27, L21 = 32, L3 = 35,
  D14 = 32, D9 = 40, E015 = 46, L2 = 44,
  L10 = 35, L5 = 70,
  M56 = 68, M70 = 37, P28 = 55, P51 = 66, P71 = 67
)
for (sample_id in names(age_map)) {
  gc_all$patient_age[gc_all$orig.ident == sample_id] <- age_map[sample_id]
}

# Age groups
gc_all$age_group <- cut(gc_all$patient_age,
                        breaks = c(0, 30, 40, 50, 60, 70, Inf),
                        labels = c("20-30", "31-40", "41-50", "51-60", "61-70", "70+"),
                        right = TRUE)

# --- 7D: ILAE score (epilepsy severity) ---
gc_all$ILAE_score <- NA
gc_all$ILAE_score[gc_all$orig.ident %in% c("B1","D11","D8","F1","H1","I1","J1","L20")] <- "0"
gc_all$ILAE_score[gc_all$orig.ident %in% c("C1","C2","E008","L14","L15","L18")]        <- "1"
gc_all$ILAE_score[gc_all$orig.ident %in% c("D12","D13","L1","L21","L3")]                <- "2"
gc_all$ILAE_score[gc_all$orig.ident %in% c("D14","D9","E015","L2")]                     <- "3"
gc_all$ILAE_score[gc_all$condition != "EPILEPSY"]                                        <- "PM"

# --- 7E: Tissue type ---
# Epilepsy samples are surgical tissue; all others are postmortem
gc_all$tissue_type <- ifelse(gc_all$condition == "EPILEPSY", "surgical", "postmortem")

# --- 7F: Group label (condition + tissue type) ---
gc_all$group_label <- paste0(gc_all$condition, " (", gc_all$tissue_type, ")")

# Verify metadata
table(gc_all$condition, gc_all$tissue_type)
table(gc_all$sex, gc_all$condition, useNA = "ifany")


###############################################################################
# SECTION 8: REMOVE LOW-QUALITY SAMPLES                                       #
#                                                                             #
# Three samples have unacceptably low median nCount (transcript counts):      #
#   M44 (HIBI): nCount = 15                                                   #
#   M25 (HIBI): nCount = 28                                                   #
#   L5 (METASTASIS): nCount = 31                                              #
#                                                                             #
# These are removed. M41 (HIBI, nCount=106) is kept as borderline acceptable  #
# (comparable to CARDIAC L10, nCount=120).                                     #
#                                                                             #
# P60 was also previously removed from gc_all.                                #
###############################################################################

samples_to_remove <- c("M44", "M25", "L5", "P60")
cells_to_keep <- colnames(gc_all)[!(gc_all$orig.ident %in% samples_to_remove)]
gc_all <- subset(gc_all, cells = cells_to_keep)

message("Samples removed: ", paste(samples_to_remove, collapse = ", "))
message("Remaining cells: ", ncol(gc_all))
table(gc_all$orig.ident, gc_all$condition)


###############################################################################
# SECTION 9: SCTRANSFORM, PCA, CLUSTERING, UMAP (gc_all)                      #
#                                                                             #
# Process the isolated GC cells with:                                         #
#   - SCTransform normalization                                                #
#   - PCA with 50 components                                                   #
#   - dims = 1:25 captures 64.1% variance                                     #
#     (PCs checked for nCount confound: only PC2 rho=-0.25 and PC6 rho=0.19  #
#      show moderate correlation — acceptable)                                 #
#   - resolution = 0.2 for clustering                                          #
#   - UMAP for visualization                                                   #
#                                                                             #
# IMPORTANT: After SCTransform, layers are split by sample. JoinLayers is     #
# required before accessing expression data (FetchData, GetAssayData).        #
###############################################################################

# --- 9A: SCTransform ---
gc_all <- SCTransform(gc_all, assay = "Xenium")

# --- 9B: PCA ---
gc_all <- RunPCA(gc_all, features = VariableFeatures(gc_all), npcs = 50, verbose = FALSE)
ElbowPlot(gc_all, ndims = 50)

# --- 9C: Clustering with dims=1:25, resolution=0.2 ---
gc_all <- FindNeighbors(gc_all, dims = 1:25)
gc_all <- FindClusters(gc_all, resolution = 0.2)
table(Idents(gc_all))

# --- 9D: UMAP ---
set.seed(42)
gc_all <- RunUMAP(gc_all, dims = 1:25)

# --- 9E: Join layers and normalize Xenium assay ---
# Required for FetchData, FeaturePlot, GetAssayData on the Xenium assay
gc_all[["Xenium"]] <- JoinLayers(gc_all[["Xenium"]])
gc_all <- NormalizeData(gc_all, assay = "Xenium")

# --- 9F: Find cluster markers ---
gc_all <- fix_sct_models(gc_all)
gc_all <- PrepSCTFindMarkers(gc_all)

gc_markers <- FindAllMarkers(gc_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gc_top25 <- gc_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 25) %>%
  ungroup()

gc_top25_list <- split(gc_top25$gene, gc_top25$cluster)

# Save
saveRDS(gc_all, gc_all_path)


###############################################################################
# SECTION 10: BASIC UMAP VISUALIZATIONS                                        #
#                                                                             #
# Standard UMAP plots colored by cluster, FOV, condition, sex, age, ILAE.    #
###############################################################################

DefaultAssay(gc_all) <- "Xenium"

# By cluster
DimPlot(gc_all, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE) +
  coord_fixed() + ggtitle("GC Clusters (dims=1:25, res=0.2)")

# By FOV
DimPlot(gc_all, reduction = "umap", group.by = "orig.ident",
        shuffle = TRUE, pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by FOV")

# By condition
DimPlot(gc_all, reduction = "umap", group.by = "condition",
        shuffle = TRUE, pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by Condition")

# By sex
DimPlot(gc_all, reduction = "umap", group.by = "sex",
        shuffle = TRUE, pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by Sex")

# By age group
DimPlot(gc_all, reduction = "umap", group.by = "age_group",
        shuffle = TRUE, pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by Age Group")

# By ILAE score
DimPlot(gc_all, reduction = "umap", group.by = "ILAE_score", shuffle = TRUE) +
  coord_fixed() + ggtitle("GC Cells by ILAE Score")

# By tissue type (surgical vs postmortem)
DimPlot(gc_all, reduction = "umap", group.by = "tissue_type",
        shuffle = TRUE, pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by Tissue Type")

# Feature plot for key GC marker
FeaturePlot(gc_all, features = "PROX1", reduction = "umap",
            min.cutoff = "q10", max.cutoff = "q90")

# Spatial feature plot using custom function (high expression on top)
spatial_feature_plot(gc_all, feature = "PROX1", fov = "M38")
spatial_feature_plot(gc_all, feature = "CALB1", fov = "M38")


###############################################################################
# SECTION 11: EXCITABILITY ANALYSIS                                            #
#                                                                             #
# Compute hypo-excitability and hyper-excitability module scores to assess    #
# the excitability state of each granule cell.                                #
#                                                                             #
# NOTE: AddModuleScore requires ctrl=5, nbin=10 for Xenium's small panel     #
# (~366 genes). Default parameters fail with "sample larger than population". #
#                                                                             #
# CAVEAT: Tissue type (surgical vs postmortem) is a major confound.           #
# Within-postmortem comparisons are valid; cross-tissue comparisons are not.  #
###############################################################################

DefaultAssay(gc_all) <- "Xenium"

# Define gene sets
hypo_genes  <- list(c("SLC12A5", "CALB1", "GABRA5", "GRIN2A", "GABRB3", "SEMA5A"))
hyper_genes <- list(c("HOMER1", "SNAP25", "KCNMA1", "DLG1", "DLG2", "GABRA2"))

# Compute module scores (ctrl=5, nbin=10 for small Xenium panel)
gc_all <- AddModuleScore(gc_all, features = hypo_genes, name = "hypo_score",
                         assay = "Xenium", ctrl = 5, nbin = 10)
gc_all <- AddModuleScore(gc_all, features = hyper_genes, name = "hyper_score",
                         assay = "Xenium", ctrl = 5, nbin = 10)

# Excitability balance: positive = hyperexcitable, negative = hypoexcitable
gc_all$excitability_balance <- gc_all$hyper_score1 - gc_all$hypo_score1

# --- Visualize excitability ---
# UMAP with hypo, hyper, and balance scores
FeaturePlot(gc_all, features = "hypo_score1", cols = c("#F0F0F0", "blue")) +
  coord_fixed() + ggtitle("Hypoexcitability Score")

FeaturePlot(gc_all, features = "hyper_score1", cols = c("#F0F0F0", "red")) +
  coord_fixed() + ggtitle("Hyperexcitability Score")

FeaturePlot(gc_all, features = "excitability_balance",
            cols = c("blue", "#F0F0F0", "red")) +
  coord_fixed() + ggtitle("Excitability Balance (red=hyper, blue=hypo)")

# Excitability by condition
ggplot(gc_all@meta.data, aes(x = condition, y = excitability_balance, fill = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic() +
  labs(x = NULL, y = "Excitability Balance\n(+ hyper / - hypo)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Excitability Balance by Condition")

# Per-sample excitability (shows within-group variability)
ggplot(gc_all@meta.data, aes(x = orig.ident, y = excitability_balance, fill = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic() +
  labs(x = NULL, y = "Excitability Balance\n(+ hyper / - hypo)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Excitability Balance by Sample")

# Excitability by ILAE score (epilepsy only)
gc_epi <- subset(gc_all, subset = condition == "EPILEPSY")
gc_epi$ILAE_numeric <- as.numeric(gc_epi$ILAE_score)

cor.test(gc_epi$excitability_balance, gc_epi$ILAE_numeric, method = "spearman")

ggplot(gc_epi@meta.data, aes(x = ILAE_score, y = excitability_balance, fill = ILAE_score)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ggtitle("Excitability Balance by ILAE Score (Epilepsy Only)")


###############################################################################
# SECTION 12: nCOUNT CONFOUND ANALYSIS                                         #
#                                                                             #
# Many genes are strongly anti-correlated with nCount (sequencing depth),     #
# meaning their apparent expression gradient is driven by technical, not      #
# biological, variation.                                                       #
#                                                                             #
# KEY FINDINGS:                                                                #
#   Confounded (DO NOT trust raw expression):                                  #
#     PROX1 (rho=-0.92), GRIA1 (rho=-0.97), SLC17A7 (rho=-0.95),            #
#     SNAP25 (rho=-0.79), GRIN3A (rho=-0.77)                                  #
#   NOT confounded (safe for cross-FOV comparisons):                           #
#     CALB1 (rho=+0.29), POSTN (rho=+0.08)                                   #
#                                                                             #
# For confounded genes, use nCount-corrected residuals (Section 13).          #
###############################################################################

# Check nCount distribution on UMAP
FeaturePlot(gc_all, features = "nCount_Xenium", cols = c("#F0F0F0", "darkblue")) +
  coord_fixed() + ggtitle("Total Transcripts (nCount)")

# Check correlation of key genes with nCount
key_genes <- c("PROX1", "CALB1", "GRIA1", "SLC17A7", "SNAP25", "GRIN3A",
               "HOMER1", "POSTN", "GRIN2A", "GABRA5")

for (gene in key_genes) {
  rho <- cor(gc_all@meta.data$nCount_Xenium,
             FetchData(gc_all, vars = gene)[[1]],
             method = "spearman")
  message(gene, ": Spearman rho with nCount = ", round(rho, 3))
}

# Check PCA components for nCount confound
pca_embeddings <- Embeddings(gc_all, "pca")
for (pc in 1:10) {
  rho <- cor(gc_all$nCount_Xenium, pca_embeddings[, pc], method = "spearman")
  if (abs(rho) > 0.15) {
    message("PC", pc, ": rho with nCount = ", round(rho, 3), " (moderate)")
  }
}


###############################################################################
# SECTION 13: CALB1 MATURITY AXIS AND HOMER1 ANTI-CORRELATION                 #
#                                                                             #
# The UMAP gradient represents a genuine biological maturity continuum:       #
#   - CALB1 marks mature GCs (NOT confounded by nCount, rho=+0.29)           #
#   - GRIN3A marks immature GCs (but IS confounded — needs correction)        #
#   - HOMER1 anti-correlates with CALB1 (rho=-0.20), reflecting an           #
#     activity-state axis independent of maturity                              #
#                                                                             #
# Maturity score = CALB1 - GRIN3A (using nCount-corrected GRIN3A residuals)  #
###############################################################################

# --- 13A: Visualize maturity markers ---
FeaturePlot(gc_all, features = c("CALB1", "GRIN3A"), ncol = 2,
            cols = c("#F0F0F0", "darkblue"))

# --- 13B: Compute nCount-corrected maturity score ---
# Raw maturity score
gc_all$maturity_raw <- FetchData(gc_all, "CALB1")[[1]] - FetchData(gc_all, "GRIN3A")[[1]]

# nCount-corrected maturity score (using linear regression residuals)
# This removes the effect of sequencing depth on the maturity estimate
gc_all$maturity_corrected <- residuals(
  lm(maturity_raw ~ nCount_Xenium, data = gc_all@meta.data)
)

# Visualize both
FeaturePlot(gc_all, features = "maturity_corrected",
            cols = c("purple", "#F0F0F0", "darkgreen")) +
  coord_fixed() + ggtitle("Maturity Score (nCount-corrected)")

# --- 13C: HOMER1/CALB1 anti-correlation (activity-state axis) ---
cor_result <- cor.test(
  FetchData(gc_all, "HOMER1")[[1]],
  FetchData(gc_all, "CALB1")[[1]],
  method = "spearman"
)
message("HOMER1 vs CALB1: rho = ", round(cor_result$estimate, 3),
        ", p = ", formatC(cor_result$p.value, format = "e", digits = 2))

# Scatter plot
ggplot(data.frame(
  HOMER1 = FetchData(gc_all, "HOMER1")[[1]],
  CALB1 = FetchData(gc_all, "CALB1")[[1]]
), aes(x = CALB1, y = HOMER1)) +
  geom_point(alpha = 0.1, size = 0.3) +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  ggtitle(paste0("HOMER1 vs CALB1 (rho = ", round(cor_result$estimate, 3), ")"))

# Maturity by condition (within postmortem only — valid comparison)
gc_pm <- subset(gc_all, subset = tissue_type == "postmortem")

ggplot(gc_pm@meta.data, aes(x = condition, y = maturity_corrected, fill = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  theme_classic() +
  scale_fill_manual(values = condition_colors) +
  labs(x = NULL, y = "Maturity Score (corrected)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Maturity Score by Condition (Postmortem Only)")


###############################################################################
# SECTION 14: POSTN BIOLOGY (CLUSTER 2)                                        #
#                                                                             #
# POSTN (Periostin) is exclusively expressed in cluster 2 and is NOT          #
# confounded by nCount (rho=+0.08). POSTN marks the most structurally        #
# stable, ECM-anchored mature granule cells.                                  #
#                                                                             #
# KEY FINDING: POSTN+ proportion decreases with ILAE severity:               #
#   ILAE 0: ~11.3% → ILAE 2: ~12.3% → ILAE 3: ~7.5%                        #
#   Interpretation: compensatory ECM reinforcement at moderate severity,     #
#   followed by decompensation at severe sclerosis.                           #
###############################################################################

# POSTN expression on UMAP
FeaturePlot(gc_all, features = "POSTN", cols = c("#F0F0F0", "#B2182B")) +
  coord_fixed() + ggtitle("POSTN Expression")

# POSTN spatial plot (high expression on top)
spatial_feature_plot(gc_all, feature = "POSTN", fov = "M38")

# Confirm POSTN is not nCount-confounded
rho_postn <- cor(gc_all$nCount_Xenium, FetchData(gc_all, "POSTN")[[1]], method = "spearman")
message("POSTN vs nCount: rho = ", round(rho_postn, 3))

# POSTN+ proportion by ILAE score (epilepsy only)
postn_threshold <- 0  # or use a quantile-based threshold
gc_epi$POSTN_positive <- FetchData(gc_epi, "POSTN")[[1]] > postn_threshold

postn_by_ilae <- gc_epi@meta.data %>%
  group_by(ILAE_score) %>%
  summarise(
    n_total = n(),
    n_postn_pos = sum(POSTN_positive, na.rm = TRUE),
    pct_postn = round(100 * n_postn_pos / n_total, 1)
  )

print(postn_by_ilae)

# DotPlot of POSTN and related ECM genes by cluster
DotPlot(gc_all, features = c("POSTN", "COL1A1", "COL1A2", "FN1", "THBS1",
                              "MMP2", "TNC", "VCAN"),
        group.by = "seurat_clusters", cols = c("blue", "red")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("ECM Genes by Cluster")


###############################################################################
# SECTION 15: SLINGSHOT PSEUDOTIME TRAJECTORY                                  #
#                                                                             #
# Slingshot infers a pseudotime trajectory from PCA embeddings.               #
# Because gc_all has ~115K cells, we downsample to 20K for Slingshot          #
# (otherwise it takes forever). The trajectory follows maturity:              #
#   immature GCs → mature GCs                                                 #
#                                                                             #
# Literature-defined maturity stages (4 stages):                              #
#   Stage 1 (immature): DCX, STMN1, DPYSL3                                   #
#   Stage 2 (early):    CALB2, GRIN3A                                          #
#   Stage 3 (mid):      GRIA1, SLC17A7                                        #
#   Stage 4 (mature):   CALB1, PROX1                                          #
#                                                                             #
# NOTE: SCE conversion must manually carry PCA embeddings (Seurat doesn't    #
# automatically transfer them).                                                #
###############################################################################

# Downsample to 20K cells for tractable computation
set.seed(42)
cells_subset <- sample(colnames(gc_all), min(20000, ncol(gc_all)))
gc_sub <- subset(gc_all, cells = cells_subset)

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(gc_sub)

# Manually transfer PCA embeddings (Seurat doesn't do this automatically)
reducedDim(sce, "PCA") <- Embeddings(gc_sub, "pca")

# Run Slingshot
sce <- slingshot(sce, clusterLabels = gc_sub$seurat_clusters,
                 reducedDim = "PCA")

# Extract pseudotime
pseudotime <- slingPseudotime(sce)

# Map back to gc_all
gc_all$pseudotime <- NA
gc_all$pseudotime[cells_subset] <- pseudotime[, 1]

# Visualize pseudotime on UMAP
FeaturePlot(gc_all, features = "pseudotime",
            cols = c("purple", "yellow")) +
  coord_fixed() + ggtitle("Slingshot Pseudotime")

# Pseudotime by cluster
ggplot(gc_all@meta.data[!is.na(gc_all$pseudotime), ],
       aes(x = seurat_clusters, y = pseudotime, fill = seurat_clusters)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ggtitle("Pseudotime by Cluster")


###############################################################################
# SECTION 16: DIFFERENTIAL GENE EXPRESSION ANALYSIS                            #
#                                                                             #
# Two approaches:                                                             #
#   16A: Single-cell DEG (FindAllMarkers) — what genes distinguish clusters   #
#   16B: Pseudobulk DEG (DESeq2) — compare conditions within each cluster    #
#                                                                             #
# IMPORTANT: Pseudobulk DEG requires ≥3 samples per group per cluster.       #
# With the current dataset, the most reliable comparison is SEPSIS vs HIBI.   #
###############################################################################

# --- 16A: Single-cell DEG by cluster ---
# Markers already computed in Section 9F (gc_markers, gc_top25_list)
# Additional visualization:

# DotPlot of top markers
top_per_cluster <- gc_markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

DotPlot(gc_all, features = unique(top_per_cluster$gene),
        group.by = "seurat_clusters", cols = c("blue", "red")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Cluster Markers")

# Heatmap of top 5 markers per cluster
gc_top5 <- gc_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5)

DoHeatmap(gc_all, features = gc_top5$gene, group.by = "seurat_clusters") +
  ggtitle("Top 5 Markers per Cluster")


# --- 16B: Pseudobulk DEG (condition comparisons per cluster) ---
# Aggregate expression by condition × sample × cluster
pseudo_gc <- AggregateExpression(
  gc_all,
  assays = "Xenium",
  return.seurat = TRUE,
  group.by = c("condition", "orig.ident", "seurat_clusters")
)

pseudo_gc$subtype_condition <- paste(pseudo_gc$seurat_clusters, pseudo_gc$condition, sep = "_")
Idents(pseudo_gc) <- "subtype_condition"

available_groups <- levels(factor(Idents(pseudo_gc)))
message("Available pseudobulk groups:\n", paste(available_groups, collapse = "\n"))

# Run pairwise condition comparisons for each cluster
clusters <- unique(gc_all$seurat_clusters)
conditions_to_compare <- c("SEPSIS", "HIBI", "MAID", "TBI")
pairs <- combn(conditions_to_compare, 2, simplify = FALSE)

all_pb_results <- list()

for (cl in clusters) {
  for (pair in pairs) {
    id1 <- paste0(cl, "_", pair[1])
    id2 <- paste0(cl, "_", pair[2])

    if (id1 %in% available_groups & id2 %in% available_groups) {
      comparison <- paste0(pair[1], "_vs_", pair[2])
      message("Running DESeq2 for cluster ", cl, ": ", comparison)
      tryCatch({
        res <- FindMarkers(
          pseudo_gc, ident.1 = id1, ident.2 = id2,
          assay = "Xenium", test.use = "DESeq2", min.cells.group = 1
        )
        res$cluster <- cl
        res$comparison <- comparison
        res$gene <- rownames(res)
        all_pb_results[[paste0(cl, "_", comparison)]] <- res
      }, error = function(e) {
        message("  Skipped: ", e$message)
      })
    }
  }
}

# Combine and save
if (length(all_pb_results) > 0) {
  all_pb_combined <- do.call(rbind, all_pb_results)
  write.csv(all_pb_combined,
            file.path(deg_dir, "GC_pseudobulk_DEG_results.csv"),
            row.names = FALSE)
  message("Pseudobulk DEG results saved.")

  # Summary of significant DEGs
  all_pb_combined %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster, comparison) %>%
    summarise(n_up = sum(avg_log2FC > 0), n_down = sum(avg_log2FC < 0), .groups = "drop") %>%
    print(n = 60)
}

# --- 16C: Within-group Kruskal-Wallis (sample heterogeneity) ---
# Tests whether samples within a condition differ from each other
conditions <- unique(gc_all$condition)
all_kw_results <- list()

for (cond in conditions) {
  message("Running KW for: ", cond)
  res <- kw_within_group(gc_all, cond)
  if (!is.null(res)) {
    all_kw_results[[cond]] <- res
    sig <- sum(res$p_adj < 0.05)
    message("  Significant genes (FDR < 0.05): ", sig)
  }
}

if (length(all_kw_results) > 0) {
  all_kw_combined <- do.call(rbind, all_kw_results)

  # View top hits per group
  all_kw_combined %>%
    filter(p_adj < 0.05) %>%
    group_by(condition) %>%
    slice_head(n = 10) %>%
    print(n = 40)
}


###############################################################################
# FINAL SAVE                                                                   #
###############################################################################

saveRDS(gc_all, gc_all_path)
message("Final gc_all object saved to: ", gc_all_path)
message("Analysis complete.")
