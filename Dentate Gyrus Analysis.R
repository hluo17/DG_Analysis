################################################################################
# Below is Mehwish's training script 2 after merging the rds
# To do some simple plotting and examine the data
# ------------------------------------------------------------------
#QC, SCTransform, PCA, UMAP, Clustering - same as Single-Sample Analysis Workflow but done on the merged file
# ------------------------------------------------------------------
#

suppressPackageStartupMessages({
  library(Seurat)        # main package for Seurat objects and plotting
  library(dplyr)         # data wrangling
  library(tidyr)         # pivot_longer / pivot_wider
  library(ggplot2)       # plotting backbone
  library(patchwork)     # combine ggplots with / and wrap_plots()
  library(scCustomize)   # custom Seurat plotting helpers, e.g. VlnPlot_scCustom
  library(colorspace)    # color manipulation if needed later
  library(grid)          # helps allow the text box for cluster selection to download properly
  library(gridExtra)     # arranging multiple plots
  library(rlang)         # helper infrastructure used by some functions
  library(stringr)       # string handling
  library(scales)        # scaling helpers for plots
  library(dittoSeq)      # used indirectly for palettes in this project
  library(RColorBrewer)  # color palette utilities
  library(tibble)  # helps with deframe()
  library(tidyverse)
  library(FSA)     # For dunnTest()
  library(ggpubr)  # For stat_pvalue_manual()
  library(broom)
  library(readxl) # For reading excel
})

# Recall merged file 
merged_obj <- readRDS("X:\\CembrowskiLab\\Hope Luo\\DG analysis\\DG_combined_raw_obj.rds")

# 1. Quality Control 
# Visualize QC metrics
# NOTE: nFeature = number of unique genes expressed in a single cell
p1 <- VlnPlot(merged_obj, features = "nFeature_Xenium", group.by = "orig.ident", pt.size = 0)
p1
# NOTE: nCount = number of total genes for each cell
p2 <- VlnPlot(merged_obj, features = "nCount_Xenium", group.by = "orig.ident", pt.size = 0)
p2

# NOTE: this is for subsetting the sample!
# NOTE: DON'T NORMALLY RUN THIS
# Filter cells based on QC thresholds (only if subsetting is needed e.g. for a region or transcripts of a specific size)
merged_obj <- subset(merged_obj, subset = nFeature_Xenium > 5 & nCount_Xenium > 0) 


# 2. Normalization and Feature Selection with SCTransform
# SCTransform performs normalization, variance stabilization, and identifies variable features
merged_obj <- SCTransform(merged_obj, assay = "Xenium")

# function to plot residual variance vs gene expression:  to visualize variance stabilization
# updated for merged obj
residualVarPlot <- function(input, xaxis = "gmean", max_resvar = 100, ntop = 20, annotate = FALSE, pt_size = 1.1) {
  
  # If a Seurat object was passed, extract gene_var from it
  # If a data frame was passed, use it directly
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
    geom_scattermore(data = subset(gene_var, gene %in% topn), pointsize = pt_size, shape = 16, alpha = 1.0, color = "deeppink") +
    geom_hline(yintercept = 1, color = "#4daf4a", linewidth = 0.9, linetype = "dashed") +
    geom_smooth(method = "loess", span = 0.1, linewidth = 0.9, formula = y ~ x, color = "#e41a1c") +
    scale_y_continuous(trans = "sqrt", breaks = c(0, 1, 10, 25, 50, 100, 150), limits = c(0, max_resvar + 1)) +
    scale_x_continuous(trans = "log10", breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = MASS::rational) +
    xlab("Gene mean expression") +
    ylab("Residual variance")
  if (annotate) {
    p <- p + geom_text_repel(
      data = subset(gene_var, gene %in% topn), aes(label = gene), color = "gray25",
      size = 1.8,
      direction = "x",
      vjust = 0.5,
      hjust = 0.5,
      segment.size = 0.2,
      segment.alpha = 0.2
    )
  }
  
  cat("Top", ntop, "genes by residual variance:\n")
  print(topn)
  
  return(p)
}

# Extracts and combines SCT attributes for merged objects
get_gene_attributes <- function(obj, assay = "SCT") {
  # Get the results (returns a list for merged objects)
  res <- SCTResults(obj, slot = "feature.attributes", assay = assay)
  
  # If it's already a single data frame, return it
  if (is.data.frame(res)) {
    res$gene <- rownames(res)
    res$model <- "Model_1"
    return(res)
  }
  
  # If it's a list, bind all data frames together
  combined_df <- do.call(rbind, lapply(names(res), function(model_name) {
    df <- res[[model_name]]
    df$gene <- rownames(df)
    df$model <- model_name # Keep track of which FOV/Sample the data comes from
    return(df)
  }))
  
  return(combined_df)
}

# test normalization by plotting "residual variance" vs "gene expression" 
# SCTResults is a way to pull data from SCTAssay object:
# documentation: https://satijalab.org/seurat/reference/sctresults#:~:text=Arguments,just%20returns%20the%20slot%20directly).
gene_attr_merged <- get_gene_attributes(merged_obj, assay = "SCT")

residualVarPlot(gene_attr_merged, max_resvar = 15, pt_size = 3, annotate = T)


# 3. Dimensionality Reduction: PCA
# PCA reduces high-dimensional data to principal components
# - features: use VariableFeatures from SCTransform
# - npcs: number of PCs to compute (e.g., 50)
# - Balance is Key: 
# - Including more PCs captures more variation in the data, but may introduce noise.
# - Including fewer PCs reduces noise but may risk missing important biological signals.
# In the funtion below we are simply setting how many pcs to calculate, not how many to use (which is what dims is for)
merged_obj <- RunPCA(merged_obj, features = VariableFeatures(merged_obj), npcs = 50, verbose = FALSE)

Reductions(merged_obj)

# Option 1: Examine ElbowPlot to choose how many PCs to use downstream
ElbowPlot(merged_obj, ndims = 50) 

# Alternatively, we can also Visualize this as Percentage Variance Explained, 
# instead of change in Standard deviation this way:

# Extract standard deviations and calculate variance (stdev^2)
pca_var <- merged_obj[["pca"]]@stdev^2

# Calculate total variance from the scaled data layer
total_var <- sum(rowVars(LayerData(merged_obj, assay = "SCT", layer = "scale.data")))

# Create a dataframe for the first 50 PCs
plot_data <- data.frame(PC = 1:50, Pct_Variance = (pca_var[1:50] / total_var) * 100)

# Generate the Elbow Plot for Variance Explained
ggplot(plot_data, aes(x = PC, y = Pct_Variance)) +
  geom_point(size = 2) +
  labs(
    title = "Elbow Plot: Variance Explained",
    x = "Principal Component", 
    y = "% Variance explained"
  ) +
  theme_classic()

# Option 2: PCA Heatmap (can also create these heatmaps after selecting the PCA above)
# - Use DimHeatmap() to visualize the top features (genes) driving each PC.
# - Adjust 'dims' to specify which PCs to examine; 'cells' controls how many cells to display per PC.
# - 'balanced = TRUE' scales positive and negative loadings equally.
DimHeatmap(merged_obj, dims = 1:8, cells = 50, balanced = TRUE)

###############
# Recommended follow Larissa's loop from previous script.
# We skip that in interest of time, but you can reference 
# our run under resources on the github.
###############



# 4. Clustering
# Identify clusters
# - resolution: higher => more clusters
merged_obj <- FindNeighbors(merged_obj, dims = 1:20)
merged_obj <- FindClusters(merged_obj, resolution = 0.2)



# 5. Non-linear Embedding: UMAP
# UMAP parameters:
# - dims: PCs to include (e.g., 1:20)
merged_obj <- RunUMAP(merged_obj, dims = 1:20)



# 6.Visualizations

# Standardising colors
# Get all cluster IDs
global_clusters <- levels(Idents(merged_obj))
n <- length(global_clusters)

# After subsetting to clusters 0-10:
global_clusters <- levels(factor(Idents(merged_obj_clean)))
n <- length(global_clusters)

cluster_colors <- Polychrome::createPalette(
  N = n,
  seedcolors = c("#FF0000", "#00FF00", "#0000FF"),
  range = c(30, 80)
)
names(cluster_palette) <- global_clusters

merged_obj_clean@misc$cluster_colors <- cluster_colors

# Store mapping for reproducibility
merged_obj@misc$cluster_colors <- cluster_colors

# Ensure Idents is a factor with consistent levels
Idents(merged_obj) <- factor(Idents(merged_obj), levels = global_clusters)

# Plot UMAPs
# NOTE: Plotting this UMAP with the merged obj will result in a large amount of 
# Clusters due to the nature of merging separate samples
# The section below this will run through harmony and integration
# So that the umap from the merged object will actually make sense
DimPlot(
  merged_obj,
  reduction = "umap",
  cols = merged_obj@misc$cluster_colors
) + DarkTheme() + coord_fixed() + ggtitle("merged_obj")

DimPlot(
  merged_obj,
  reduction = "umap",
  cols = merged_obj@misc$cluster_colors,
  group.by = "Group"
) + DarkTheme() + coord_fixed() + ggtitle("merged_obj")

# Plot spatial FOV view
ImageDimPlot(
  merged_obj,
  fov = "X3fov",
  cols = merged_obj@misc$cluster_colors,
  size = 0.5,
  border.size = NA,
  axes = TRUE,
  dark.background = TRUE
) + DarkTheme()


# You can also perform subseting as discussed in previous section!




############# SECTION: HARMONY INTEGRATION ####################################################

# --- Merge vs. Integration: When to Use Which ---
# Merge (Seurat::merge):
#   • Combines multiple Seurat objects into one.
#   • Ideal when samples share similar technical conditions and batch effects are minor.
# Integration (Harmony):
#   • Actively corrects technical or batch differences between datasets.
#   • Essential when combining data from different experimental dates, sample preps, timepoints.
#   • Aligns shared biological signals across batches, improving clustering and visualization.
#   • For some datasets, it can wash away TRUE biological differences (e.g. integration healthy vs disease conditions).
# Decide what is best for your dataset!!

# For more detailed coverage, paper on Harmony: https://www.nature.com/articles/s41592-019-0619-0

# Note that for this dataset, Harmony is not quite applicable for our data, as it has
# very similar samples, and low batch effects. The purpose of the following script
# is to showcase how Harmony integration would be performed

# First lets clear previously loaded items
rm(list = ls())
gc()


# And using original merged object from the previous section 
# ensure the merged obj did not go through any manipulation
merged_obj <- readRDS("DG_merged_raw_obj.rds")

#Normalize using SCTransform 
int_obj<- SCTransform(merged_obj,assay = "Xenium", verbose = FALSE)

#Run PCA 
int_obj<- RunPCA(int_obj,features = rownames(merged_obj), npcs = 100,
                 verbose = FALSE)

#Select PCs via ElbowPlot (or other methods discussed above)
ElbowPlot(int_obj,ndims = 100) 

# Integrate samples with Harmony
#    - orig.reduction: PCA embeddings to correct
#    - new.reduction: stores the batch-corrected embeddings
#    - normalization.method: use SCT-normalized data
int_obj<- IntegrateLayers(
  object = int_obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = TRUE)

#Run UMAP, find neighbours and clustering
int_obj<- RunUMAP(int_obj, reduction = "harmony", dims = 1:20, verbose = FALSE)
int_obj<- FindNeighbors(int_obj, reduction = "harmony", dims = 1:20)
int_obj<- FindClusters(int_obj, resolution = 0.2)



#Perform visualization and subseting as discussed in previous sections!

# Standardising colors
# Get all cluster IDs
global_clusters <- levels(Idents(merged_obj))
n <- length(global_clusters)

# Generate distinct color palette (Polychrome handles up to 40–50 unique colors)
if (n <= 36) {
  palette <- Polychrome::palette36.colors(n)
} else {
  # For >36 clusters, extend smoothly with hue palette
  palette <- scales::hue_pal()(n)
}

# Create a named color vector for consistent mapping
cluster_colors <- setNames(palette, global_clusters)

# Store mapping for reproducibility
int_obj@misc$cluster_colors <- cluster_colors

# Ensure Idents is a factor with consistent levels
Idents(int_obj) <- factor(Idents(int_obj), levels = global_clusters)

# Check cell counts per cluster
table(int_obj$seurat_clusters)

# Groups cells by the metadata "group"
# NOTE: change "group" based on what is labelled in the metadata
DimPlot(
  int_obj,
  reduction = "umap",
  cols = int_obj@misc$cluster_colors,
  group.by = "Group"
) + DarkTheme() + coord_fixed() + ggtitle("int_obj")

### NOTE: the above clustering gave me too many clusters with very smaller cell counts
# Will need to subset those clusters out:

# First inspect QC metrics (transcripts per cell, gene count) across clusters
# To determine which clusters to subset out
VlnPlot(int_obj, features = "nCount_Xenium", group.by = "seurat_clusters", pt.size = 0.1)

VlnPlot(int_obj, features = "nFeature_Xenium", group.by = "seurat_clusters", pt.size = 0.1)


# 1. Run FindClusters at 0.1 resolution
# (Uses existing SNN graph from FindNeighbors)
merged_obj <- FindClusters(merged_obj, resolution = 0.2)

# 2. Verify that SCT_snn_res.0.1 was created and check cluster IDs
table(Idents(merged_obj))

# 3. Subset to keep ONLY clusters 0 through 10
merged_obj_clean <- subset(merged_obj, idents = as.character(0:10))

# 4. Verify remaining clusters in the clean object
table(Idents(merged_obj_clean))

# 5. Re-run PCA and UMAP so the layout recalculates without noise clusters
merged_obj_clean <- RunPCA(merged_obj_clean, verbose = FALSE)
merged_obj_clean <- RunUMAP(merged_obj_clean, dims = 1:20)

# 6. Plot the clean UMAP
DimPlot(merged_obj_clean, label = TRUE)

# Plot spatial FOV view
# NOTE: will need to change the fov name accordingly
# NOTE: Try out a few different fovs
# Extract the default colors DimPlot uses
n <- length(levels(factor(Idents(merged_obj_clean))))
default_colors <- scales::hue_pal()(n)
names(default_colors) <- levels(factor(Idents(merged_obj_clean)))

# Spatial plot using the same default palette
ImageDimPlot(merged_obj_clean, 
             fov = "M38",
             cols = default_colors,
             size = 0.5, 
             border.size = NA)


#When subseting integrated data, make sure switch assay to 'Xenium' and 
#follow block of code above: run SCT, PCA, integratelayers, dimensionality reduction, clustering etc. 

# Force all SCT models to reference Xenium
for (model_name in names(merged_obj_clean[["SCT"]]@SCTModel.list)) {
  slot(merged_obj_clean[["SCT"]]@SCTModel.list[[model_name]], "umi.assay") <- "Xenium"
}

# Now prep and find markers
merged_obj_clean <- PrepSCTFindMarkers(merged_obj_clean)

all_markers <- FindAllMarkers(merged_obj_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top25 <- all_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 25) %>%
  ungroup()

top25_list <- split(top25$gene, top25$cluster)
top25_list

#####################################################################################
# Below is the modified trainign script 3 to run SingleR analysis

###############################################################################
# LOAD THOMPSON HUMAN HIPPOCAMPUS REFERENCE                                   #
###############################################################################

library(humanHippocampus2024)
library(ExperimentHub)
library(SingleR)
library(SummarizedExperiment)
library(SingleCellExperiment)

# Load the snRNA-seq reference (already downloaded to cache)
eh <- ExperimentHub()
sce_hpc <- eh[["EH9606"]]

# Use mid.cell.class for the primary annotation
# (ExcN, InhN, Astro, Oligo, OPC, Micro/Macro/T, CSF, Vascular)
ref_labels <- sce_hpc$mid.cell.class

# NOTE: this level contains granule cell markings
# Try out the line above first to get a general labelling, 
# and then use this line
#ref_labels <- sce_hpc$fine.cell.class

# Extract the log-normalized expression matrix
ref_expr <- logcounts(sce_hpc)

message("Reference loaded: ", ncol(sce_hpc), " cells, ",
        length(unique(ref_labels)), " cell types")
message("Cell types: ", paste(sort(unique(ref_labels)), collapse = ", "))


###############################################################################
# NORMALIZE QUERY DATA                                                        #
###############################################################################

merged_obj_clean <- NormalizeData(merged_obj_clean, assay = "Xenium", verbose = FALSE)


###############################################################################
# RUN SINGLER PER FOV                                                         #
###############################################################################

output_dir <- "X:/CembrowskiLab/Hope Luo/DG analysis/SingleR_output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

sample_names <- unique(merged_obj_clean$orig.ident)

all_labels <- rep(NA, ncol(merged_obj_clean))
all_pruned <- rep(NA, ncol(merged_obj_clean))
names(all_labels) <- colnames(merged_obj_clean)
names(all_pruned) <- colnames(merged_obj_clean)

# NOTE: shouldnt have to install BioNeighbours again but keeping the line below for now,
# just in case
# BiocManager::install("BiocNeighbors", type = "win.binary")
library(SingleR)

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
  
  ref_sub <- ref_expr[common_genes, , drop = FALSE]
  query_sub <- query_expr[common_genes, , drop = FALSE]
  
  singleR_results <- SingleR(
    test = as.matrix(query_sub),
    ref = as.matrix(ref_sub),
    labels = ref_labels
  )
  
  saveRDS(singleR_results,
          file = file.path(output_dir, paste0("SingleR_", sample_name, "_results.rds")))
  
  sample_cells <- colnames(sample_obj)
  all_labels[sample_cells] <- singleR_results$labels
  all_pruned[sample_cells] <- singleR_results$pruned.labels
  
  message("  Done. Unique labels: ", length(unique(singleR_results$labels)))
}

merged_obj_clean$SingleR_label <- all_labels
merged_obj_clean$SingleR_pruned <- all_pruned

saveRDS(merged_obj_clean,
        file = file.path(output_dir, "DG_merged_clean_SingleR.rds"))
message("Annotated object saved.")
gc()


###############################################################################
# VISUALIZATION                                                                #
###############################################################################

table(merged_obj_clean$SingleR_pruned)

# Compare SingleR labels against your Seurat clusters
table(merged_obj_clean$seurat_clusters, merged_obj_clean$SingleR_pruned)

# UMAP colored by SingleR labels
DimPlot(merged_obj_clean, reduction = "umap", group.by = "SingleR_pruned",
        label = TRUE, repel = TRUE) +
  coord_fixed() + ggtitle("SingleR Cell Type Annotations (mid.cell.class)")

# Spatial plot for one FOV
ImageDimPlot(merged_obj_clean, fov = "M38",
             group.by = "SingleR_pruned",
             size = 0.5, border.size = NA)

# Proportional bar chart per FOV
ggplot(merged_obj_clean@meta.data, aes(x = orig.ident, fill = SingleR_pruned)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(x = "FOV", y = "Proportion of Cells", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cell Type Composition by FOV")


#############################################
# Above was the code to add labels of the general cell types
# After testing that it functions correctly, below is the code
# to add both the general and fine cell labelling onto the sample

# Use superfine labels for the reference — this gives you GC.1 through GC.5
ref_labels_superfine <- sce_hpc$superfine.cell.class

# Also keep fine labels for the broad cell type identification
ref_labels_fine <- sce_hpc$fine.cell.class

# Pre-allocate for both annotation levels
all_fine   <- rep(NA, ncol(merged_obj_clean))
all_superfine <- rep(NA, ncol(merged_obj_clean))
all_fine_pruned   <- rep(NA, ncol(merged_obj_clean))
all_superfine_pruned <- rep(NA, ncol(merged_obj_clean))
names(all_fine) <- names(all_superfine) <- colnames(merged_obj_clean)
names(all_fine_pruned) <- names(all_superfine_pruned) <- colnames(merged_obj_clean)

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
  
  ref_sub <- ref_expr[common_genes, , drop = FALSE]
  query_sub <- query_expr[common_genes, , drop = FALSE]
  
  # Run SingleR at FINE level (GC, CA1, Astro, Oligo, etc.)
  sr_fine <- SingleR(
    test = as.matrix(query_sub),
    ref = as.matrix(ref_sub),
    labels = ref_labels_fine
  )
  
  # Run SingleR at SUPERFINE level (GC.1, GC.2, GC.3, GC.4, GC.5, etc.)
  sr_superfine <- SingleR(
    test = as.matrix(query_sub),
    ref = as.matrix(ref_sub),
    labels = ref_labels_superfine
  )
  
  sample_cells <- colnames(sample_obj)
  all_fine[sample_cells]   <- sr_fine$labels
  all_fine_pruned[sample_cells] <- sr_fine$pruned.labels
  all_superfine[sample_cells]   <- sr_superfine$labels
  all_superfine_pruned[sample_cells] <- sr_superfine$pruned.labels
  
  # Save per-sample results for both annotation levels
  saveRDS(sr_fine,
          file = file.path(output_dir, paste0("SingleR_fine_", sample_name, "_results.rds")))
  saveRDS(sr_superfine,
          file = file.path(output_dir, paste0("SingleR_superfine_", sample_name, "_results.rds")))
  
  message("  Done.")
}

# Add all labels to metadata
merged_obj_clean$SR_fine   <- all_fine
merged_obj_clean$SR_fine_pruned <- all_fine_pruned
merged_obj_clean$SR_superfine   <- all_superfine
merged_obj_clean$SR_superfine_pruned <- all_superfine_pruned

# Save
saveRDS(merged_obj_clean,
        file = file.path(output_dir, "DG_merged_clean_SingleR_with_all_labels.rds"))


#################################
# VISUALIZATION

# UMAP colored by SingleR labels
DimPlot(merged_obj_clean, reduction = "umap", group.by = "SR_superfine_pruned",
        label = TRUE, repel = TRUE) +
  coord_fixed() + ggtitle("SingleR Cell Type Annotations (super.fine.cell.class)")

# Spatial plot for one FOV
ImageDimPlot(merged_obj_clean, fov = "M38",
             group.by = "SR_superfine_pruned",
             size = 0.5, border.size = NA)

# Proportional bar chart per FOV
ggplot(merged_obj_clean@meta.data, aes(x = orig.ident, fill = SR_superfine_pruned)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(x = "FOV", y = "Proportion of Cells", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cell Type Composition by FOV")

#################################
# Below are the same 3 plots as above, except they are all
# Adjusted into feature plots to only highlight the GC cell types

# Define GC colors, everything else is grey
gc_subtypes <- c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5")

all_types <- unique(merged_obj_clean$SR_superfine_pruned)
highlight_cols <- setNames(rep("lightgrey", length(all_types)), all_types)
highlight_cols["GC.1"] <- "#E64B35"
highlight_cols["GC.2"] <- "#4DBBD5"
highlight_cols["GC.3"] <- "#00A087"
highlight_cols["GC.4"] <- "#3C5488"
highlight_cols["GC.5"] <- "#F39B7F"

# 1. UMAP — GC subtypes highlighted, rest greyed out
DimPlot(merged_obj_clean, reduction = "umap", group.by = "SR_superfine_pruned",
        cols = highlight_cols, label = FALSE, 
        order = gc_subtypes) +
  coord_fixed() + ggtitle("GC Subtypes Highlighted")

# 2. Spatial plot — GC subtypes highlighted
ImageDimPlot(merged_obj_clean, fov = "P60",
             group.by = "SR_superfine_pruned",
             cols = highlight_cols,
             size = 0.5, border.size = NA)

# 3. Bar chart — collapse non-GC types into "Other"
bar_data <- merged_obj_clean@meta.data %>%
  mutate(cell_group = ifelse(SR_superfine_pruned %in% gc_subtypes, 
                             SR_superfine_pruned, "Other"))

bar_colors <- c(
  "GC.1" = "#E64B35", "GC.2" = "#4DBBD5", "GC.3" = "#00A087",
  "GC.4" = "#3C5488", "GC.5" = "#F39B7F", "Other" = "lightgrey"
)

ggplot(bar_data, aes(x = orig.ident, fill = cell_group)) +
  geom_bar(position = "fill") +
  theme_classic() +
  scale_fill_manual(values = bar_colors) +
  labs(x = "FOV", y = "Proportion of Cells", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("GC Subtype Proportions by FOV")

################################
### Then subset to just granule cells and visualize through plots

# Keep only the actual GC subtypes
gc_obj <- subset(gc_obj, subset = SR_superfine_pruned %in% c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5"))

# Verify it's clean now
Idents(gc_obj) <- "SR_superfine_pruned"
table(Idents(gc_obj))

# Check how many GC cells and what subtypes
table(gc_obj$SR_superfine_pruned)

# Re-process just the GC cells
gc_obj <- SCTransform(gc_obj, assay = "Xenium")
gc_obj <- RunPCA(gc_obj, verbose = FALSE)
set.seed(42)
gc_obj <- RunUMAP(gc_obj, dims = 1:15)


### Calculate only the upregulated genes
# Fix SCT models before finding markers
for (i in seq_along(gc_obj[["SCT"]]@SCTModel.list)) {
  slot(gc_obj[["SCT"]]@SCTModel.list[[i]], "umi.assay") <- "Xenium"
}
gc_obj <- PrepSCTFindMarkers(gc_obj)

# Find markers distinguishing GC subtypes
DefaultAssay(gc_obj) <- "SCT"
gc_markers <- FindAllMarkers(gc_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

# Top 25 per subtype
gc_top25 <- gc_markers %>% 
  group_by(cluster) %>% 
  arrange(p_val_adj) %>% 
  slice_head(n = 25)

gc_top25_list <- split(gc_top25$gene, gc_top25$cluster)
gc_top25_list

### Calculate only the downregulated genes
# Fix SCT models before finding markers
for (i in seq_along(gc_obj[["SCT"]]@SCTModel.list)) {
  slot(gc_obj[["SCT"]]@SCTModel.list[[i]], "umi.assay") <- "Xenium"
}
gc_obj <- PrepSCTFindMarkers(gc_obj)

# Find ALL markers (both positive and negative fold changes)
DefaultAssay(gc_obj) <- "SCT"
gc_markers <- FindAllMarkers(
  gc_obj, 
  only.pos = FALSE,     # Set to FALSE to capture downregulated genes
  min.pct = 0.1, 
  logfc.threshold = 0.1
)

# Top 25 most significantly downregulated genes per subtype
gc_top25_down <- gc_markers %>% 
  filter(avg_log2FC < 0) %>%     # Keep only downregulated genes
  group_by(cluster) %>% 
  arrange(avg_log2FC) %>%        # Sort by strongest negative fold change (or p_val_adj)
  slice_head(n = 25)

gc_top25_down_list <- split(gc_top25_down$gene, gc_top25_down$cluster)
gc_top25_down_list

# UMAP colored by GC subtype
DimPlot(gc_obj, reduction = "umap", group.by = "SR_superfine_pruned",
        label = TRUE, repel = TRUE) +
  coord_fixed() + ggtitle("DG Granule Cell Subtypes")

# UMAP colored by FOV
DimPlot(
  gc_obj, 
  reduction = "umap", 
  group.by = "orig.ident"  
) + ggtitle("UMAP Grouped by FOV")

# UMAP colored by condition
DimPlot(
  gc_obj, 
  reduction = "umap", 
  group.by = "condition"  
) + ggtitle("UMAP Grouped by condition")

# Spatial plot of GC subtypes in one FOV
ImageDimPlot(gc_obj, fov = "M38", group.by = "SR_superfine_pruned",
             size = 0.5, border.size = NA)

#########################################################
# Adjusted training script 4
###############################################################################
# DIFFERENTIAL GENE EXPRESSION (DEG) ANALYSIS — DG GRANULE CELLS             #
#                                                                             #
# Run this AFTER the SingleR + GC subsetting steps so that gc_obj exists      #
# with SR_superfine_pruned labels (GC.1–GC.5) and has been re-processed      #
# (SCTransform, PCA, UMAP, PrepSCTFindMarkers).                              #
#                                                                             #
# Two approaches:                                                             #
#   1. Single-cell DEG (FindAllMarkers) — what genes distinguish GC subtypes  #
#   2. Pseudobulk DEG (DESeq2) — compare GC subtypes across M vs P groups    #
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scattermore)
  library(DESeq2)
})

options(future.globals.maxSize = 250000 * 1024^2)


###############################################################################
# PREREQUISITES                                                                #
#                                                                             #
# gc_obj should already be in your environment with:                          #
#   - Only GC.1–GC.5 cells (subset from merged_obj_clean)                    #
#   - SCTransform, PCA, UMAP done                                            #
#   - SCT models fixed (umi.assay = "Xenium") and PrepSCTFindMarkers run     #
#   - Idents set to SR_superfine_pruned                                       #
#                                                                             #
# If not, reload and re-subset:                                               #
# merged_obj_clean <- readRDS("X:/CembrowskiLab/Hope Luo/DG analysis/SingleR_output/DG_merged_clean_SingleR_with_all_labels.rds")
# gc_obj <- subset(merged_obj_clean, subset = SR_superfine_pruned %in% c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5"))
# Idents(gc_obj) <- "SR_superfine_pruned"
# gc_obj <- SCTransform(gc_obj, assay = "Xenium")
# gc_obj <- RunPCA(gc_obj, verbose = FALSE)
# set.seed(42)
# gc_obj <- RunUMAP(gc_obj, dims = 1:15)
# for (i in seq_along(gc_obj[["SCT"]]@SCTModel.list)) {
#   slot(gc_obj[["SCT"]]@SCTModel.list[[i]], "umi.assay") <- "Xenium"
# }
# gc_obj <- PrepSCTFindMarkers(gc_obj)
###############################################################################

# Verify setup
stopifnot("SR_superfine_pruned" %in% colnames(gc_obj@meta.data))
Idents(gc_obj) <- "SR_superfine_pruned"
table(Idents(gc_obj))

output_dir <- "X:/CembrowskiLab/Hope Luo/DG analysis/GC_DEG_output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


###############################################################################
#                                                                             #
#   SECTION 1: SINGLE-CELL LEVEL DEG — WHAT DISTINGUISHES EACH GC SUBTYPE?   #
#                                                                             #
#   FindAllMarkers compares each GC subtype against all other GC subtypes     #
#   to find the genes that make GC.1 different from GC.2, etc.                #
#                                                                             #
###############################################################################

gc_obj <- NormalizeData(gc_obj, assay = "Xenium")

# Find markers for each GC subtype vs all others
# Lower thresholds because GC subtypes are closely related
gc_markers <- FindAllMarkers(
  gc_obj,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.1
)

# Check how many markers were found per subtype
table(gc_markers$cluster)

# Top 25 markers per GC subtype by adjusted p-value
gc_top25 <- gc_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 25)

gc_top25_list <- split(gc_top25$gene, gc_top25$cluster)
gc_top25_list

# Save full marker table
# Errored - will skip for now (09/28/26)
write.csv(gc_markers,
          file.path(output_dir, "GC_subtype_all_markers.csv"),
          row.names = FALSE)

# Filter for biologically meaningful DEGs
gc_markers_filtered <- subset(
  gc_markers,
  pct.1 > 0.1 & (pct.2 < 0.3 | abs(avg_log2FC) > 0.5)
)

# Also skipped this for now (09/28/26)
write.csv(gc_markers_filtered,
          file.path(output_dir, "GC_subtype_markers_filtered.csv"),
          row.names = FALSE)

# Top marker per subtype by fold change
top1_gc <- gc_markers_filtered %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC)
top_gc_genes <- unique(top1_gc$gene)

top_gc_genes

###############################################################################
# VISUALIZE SINGLE-CELL GC DEGs                                               #
###############################################################################

# Dot plot: expression of top markers across GC subtypes
DotPlot(
  gc_obj,
  features = top_gc_genes,
  group.by = "SR_superfine_pruned",
  cols = c("red", "blue")
) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Marker per GC Subtype")

# Violin plots of top 3 markers
# Not used for now - not sure how useful this is in this case
if (length(top_gc_genes) >= 3) {
  p1 <- VlnPlot(gc_obj, features = top_gc_genes[1],
                group.by = "SR_superfine_pruned", pt.size = 0)
  p2 <- VlnPlot(gc_obj, features = top_gc_genes[2],
                group.by = "SR_superfine_pruned", pt.size = 0)
  p3 <- VlnPlot(gc_obj, features = top_gc_genes[3],
                group.by = "SR_superfine_pruned", pt.size = 0)
  p1 + p2 + p3
}

DefaultAssay(gc_obj) <- "Xenium"

# Feature plot on the granule cell umap
# featuring the top marker gene for each GC superfine cell type
FeaturePlot(
  gc_obj,
  features = top_gc_genes[1:min(5, length(top_gc_genes))],
  ncol = 2,
  min.cutoff = "q10",
  max.cutoff = "q90",
  reduction = "umap"
)

# GC subtype UMAP
DimPlot(gc_obj, reduction = "umap", group.by = "SR_superfine_pruned",
        label = TRUE, repel = TRUE, pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Subtypes") + theme_classic()

# Spatial plot of GC subtypes in one FOV
ImageDimPlot(gc_obj, fov = "P60", group.by = "SR_superfine_pruned",
             size = 0.5, border.size = NA) +
  ggtitle("P60 — GC Subtypes")

# Heatmap of top 5 markers per subtype (save to PDF if plot window is too small)
# Note: not used, not sure how this heatmap is helpful
gc_top5 <- gc_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5)

pdf(file.path(output_dir, "GC_subtype_heatmap.pdf"), width = 14, height = 10)
DoHeatmap(gc_obj, features = gc_top5$gene, group.by = "SR_superfine_pruned") +
  ggtitle("Top 5 Markers per GC Subtype")
dev.off()


###############################################################################
#                                                                             #
#   SECTION 2: PSEUDOBULK LEVEL DEG — COMPARE GC SUBTYPES ACROSS CONDITIONS  #
#                                                                             #
#   For each GC subtype, compare expression between M_group and P_group.      #
#   This asks: does GC.1 behave differently in M samples vs P samples?        #
#                                                                             #
#   NOTE: DESeq2 needs a minimum of 3 samples per group. With 6 M samples    #
#   and 4 P samples, you have enough — but only if a given GC subtype has     #
#   cells in at least 3 samples per group.                                    #
#                                                                             #
###############################################################################

# ── IMPORTANT: Define your experimental conditions ──
gc_obj$condition <- NA
gc_obj$condition[gc_obj$orig.ident %in% c("M70", "P28")] <- "MAID_group"
gc_obj$condition[gc_obj$orig.ident %in% c("P51", "P60", "P71", "M56")] <- "SEPSIS_group"
gc_obj$condition[gc_obj$orig.ident %in% c("M44", "M25", "M41")] <- "HIBI_group"
gc_obj$condition[gc_obj$orig.ident %in% c("M38")] <- "TBI_group"

# Check cell counts per subtype per condition
table(gc_obj$SR_superfine_pruned, gc_obj$condition)

# Stacked bar plot of the results from the table above
ggplot(gc_obj@meta.data, aes(x = condition, fill = SR_superfine_pruned)) +
  geom_bar(position = "fill") +
  theme_classic() +
  labs(x = "Condition", y = "Proportion of Cells", fill = "GC Subtype") +
  scale_fill_manual(values = c(
    "GC.1" = "#E64B35", "GC.2" = "#4DBBD5", "GC.3" = "#00A087",
    "GC.4" = "#3C5488", "GC.5" = "#F39B7F"
  )) +
  ggtitle("GC Subtype Proportions by Condition")

# Stacked bar plot of GC proportions per FOV 
ggplot(gc_obj@meta.data, aes(x = orig.ident, fill = SR_superfine_pruned)) +
  geom_bar(position = "fill") +
  theme_classic() +
  labs(x = "Sample", y = "Proportion of Cells", fill = "GC Subtype") +
  scale_fill_manual(values = c(
    "GC.1" = "#E64B35", "GC.2" = "#4DBBD5", "GC.3" = "#00A087",
    "GC.4" = "#3C5488", "GC.5" = "#F39B7F"
  )) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("GC Subtype Proportions by Sample")

# Check cell counts per sample per subtype (important for pseudobulk power)
message("Cells per sample per GC subtype:")
print(table(gc_obj$orig.ident, gc_obj$SR_superfine_pruned))

# Keep only cells with assigned condition
# Will error if all cells are already assigned to a condition
gc_pseudo_obj <- subset(gc_obj, subset = !is.na(condition))


###############################################################################
# CONSTRUCT THE PSEUDOBULK OBJECT                                             #
###############################################################################

pseudo_gc <- AggregateExpression(
  gc_obj,
  assays = "Xenium",
  return.seurat = TRUE,
  group.by = c("condition", "orig.ident", "SR_superfine_pruned")
)

# Create comparison identity: GC subtype + condition
pseudo_gc$subtype_condition <- paste(pseudo_gc$SR_superfine_pruned, pseudo_gc$condition, sep = "_")
Idents(pseudo_gc) <- "subtype_condition"

# Inspect available groups
message("Available pseudobulk groups:")
print(table(Idents(pseudo_gc)))


###############################################################################
# RUN PSEUDOBULK DEG FOR EACH GC SUBTYPE ACROSS CONDITIONS                    #
###############################################################################

available_groups <- levels(factor(Idents(pseudo_gc)))
gc_subtypes <- c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5")
conditions <- c("SEPSIS-group", "HIBI-group", "MAID-group", "TBI-group")

# Generate all pairwise condition combinations
pairs <- combn(conditions, 2, simplify = FALSE)

all_gc_pb_results <- list()

for (gc_sub in gc_subtypes) {
  for (pair in pairs) {
    id1 <- paste0(gc_sub, "_", pair[1])
    id2 <- paste0(gc_sub, "_", pair[2])
    
    if (id1 %in% available_groups & id2 %in% available_groups) {
      comparison <- paste0(pair[1], "_vs_", pair[2])
      message("Running DESeq2 for: ", gc_sub, " (", comparison, ")")
      tryCatch({
        res <- FindMarkers(
          object = pseudo_gc,
          ident.1 = id1,
          ident.2 = id2,
          assay = "Xenium",
          test.use = "DESeq2",
          min.cells.group = 1
        )
        res$gc_subtype <- gc_sub
        res$comparison <- comparison
        res$gene <- rownames(res)
        all_gc_pb_results[[paste0(gc_sub, "_", comparison)]] <- res
        message("  Found ", nrow(res), " DEGs")
      }, error = function(e) {
        message("  Skipped: ", e$message)
      })
    } else {
      message("  Skipping ", gc_sub, " ", pair[1], " vs ", pair[2], " — missing group")
    }
  }
}

# Combine and save
if (length(all_gc_pb_results) > 0) {
  all_gc_pb_combined <- do.call(rbind, all_gc_pb_results)
  write.csv(all_gc_pb_combined,
            file.path(output_dir, "GC_pseudobulk_DEG_results.csv"),
            row.names = FALSE)
  message("Pseudobulk DEG results saved.")
  
  top_pb <- all_gc_pb_combined %>%
    group_by(gc_subtype, comparison) %>%
    arrange(p_val_adj) %>%
    slice_head(n = 5)
  print(top_pb, n = 60)
}

results <- read.csv(file.path(output_dir, "GC_pseudobulk_DEG_results.csv"))

# 1. How many significant DEGs per subtype per comparison?
results %>%
  filter(p_val_adj < 0.05) %>%
  group_by(gc_subtype, comparison) %>%
  summarise(
    n_up = sum(avg_log2FC > 0),
    n_down = sum(avg_log2FC < 0),
    n_total = n()
  ) %>%
  print(n = 60)

# 2. Top significant genes for the reliable comparison (SEPSIS vs HIBI)
results %>%
  filter(comparison == "SEPSIS-group_vs_HIBI-group" & p_val_adj < 0.05) %>%
  group_by(gc_subtype) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5) %>%
  select(gc_subtype, gene, avg_log2FC, p_val_adj) %>%
  print(n = 30)

# 3. Are the same genes changing across multiple subtypes?
results %>%
  filter(comparison == "SEPSIS-group_vs_HIBI-group" & p_val_adj < 0.05) %>%
  group_by(gene) %>%
  summarise(
    n_subtypes = n(),
    subtypes = paste(gc_subtype, collapse = ", "),
    mean_log2FC = round(mean(avg_log2FC), 2)
  ) %>%
  filter(n_subtypes >= 3) %>%
  arrange(desc(n_subtypes))

# 4. Genes changing in only ONE subtype (subtype-specific response to injury)
results %>%
  filter(comparison == "SEPSIS-group_vs_HIBI-group" & p_val_adj < 0.05) %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  arrange(gc_subtype, p_val_adj) %>%
  select(gc_subtype, gene, avg_log2FC, p_val_adj) %>%
  print(n = 30)

################
### Additional Visualization

# Load results if needed
results <- read.csv(file.path(output_dir, "GC_pseudobulk_DEG_results.csv"))

# Filter to the reliable comparison
sig_results <- results %>%
  filter(comparison == "SEPSIS-group_vs_HIBI-group" & p_val_adj < 0.05)


###############################################################################
# 1. BUBBLE PLOT — log2FC vs gene, sized by significance, colored by subtype  #
###############################################################################

top_bubble <- sig_results %>%
  group_by(gc_subtype) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

ggplot(top_bubble, aes(x = avg_log2FC, y = reorder(gene, avg_log2FC), 
                       size = -log10(p_val_adj), color = gc_subtype)) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  scale_size_continuous(name = "-log10(p_adj)", range = c(2, 8)) +
  scale_color_manual(values = c("GC.1" = "#E64B35", "GC.2" = "#4DBBD5", 
                                "GC.3" = "#00A087", "GC.4" = "#3C5488", 
                                "GC.5" = "#F39B7F")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "avg log2FC (SEPSIS vs HIBI)", y = NULL, color = "GC Subtype") +
  ggtitle("Top DEGs per GC Subtype (SEPSIS vs HIBI)")


###############################################################################
# VISUALIZE PSEUDOBULK GC DEGs                                                #
###############################################################################

Idents(pseudo_gc) <- "subtype_condition"
DefaultAssay(pseudo_gc) <- "Xenium"

# Boxplot function for GC pseudobulk comparisons
# Check the actual group names
print(available_groups)

# Run pairwise comparisons for the two groups with enough samples
all_gc_pb_results <- list()

for (gc_sub in gc_subtypes) {
  id1 <- paste0(gc_sub, "_SEPSIS-group")
  id2 <- paste0(gc_sub, "_HIBI-group")
  
  if (id1 %in% available_groups & id2 %in% available_groups) {
    message("Running DESeq2 for: ", gc_sub, " (SEPSIS vs HIBI)")
    tryCatch({
      res <- FindMarkers(
        object = pseudo_gc,
        ident.1 = id1,
        ident.2 = id2,
        assay = "Xenium",
        test.use = "DESeq2",
        min.cells.group = 2
      )
      res$gc_subtype <- gc_sub
      res$gene <- rownames(res)
      all_gc_pb_results[[gc_sub]] <- res
      message("  Found ", nrow(res), " DEGs")
    }, error = function(e) {
      message("  Skipped ", gc_sub, ": ", e$message)
    })
  } else {
    message("  Skipping ", gc_sub, " — not in both conditions")
  }
}

plot_gc_box <- function(pseudo, gene, gc_subtype, group1 = "SEPSIS-group", group2 = "HIBI-group") {
  ids <- c(paste(gc_subtype, group1, sep = "_"),
           paste(gc_subtype, group2, sep = "_"))
  
  df <- FetchData(
    object = pseudo,
    vars = c("subtype_condition", "condition", gene)
  )
  
  colnames(df)[ncol(df)] <- "expr"
  df <- df %>% filter(subtype_condition %in% ids)
  
  ggplot(df, aes(x = subtype_condition, y = expr, color = condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.35) +
    geom_point(size = 2.8, position = position_jitter(width = 0.08, height = 0)) +
    theme_classic() +
    labs(
      title = paste0(gene, " in ", gc_subtype),
      x = NULL,
      y = "Aggregated expression"
    ) +
    scale_color_manual(values = c("SEPSIS-group" = "#E64B35", "HIBI-group" = "#4DBBD5"))
}

# Example: plot top DEGs for GC.1 (if results exist)
# Uncomment and adjust after seeing your results:
#
if ("GC.1" %in% names(all_gc_pb_results)) {
  gc1_top_genes <- rownames(head(all_gc_pb_results[["GC.1"]], 12))
  boxplots <- lapply(gc1_top_genes, function(g) {
    plot_gc_box(pseudo_gc, gene = g, gc_subtype = "GC.1")
  })
  boxplots[[1]] + boxplots[[2]] + boxplots[[3]] + boxplots[[4]] + boxplots[[5]] + boxplots[[6]] + boxplots[[7]] + boxplots[[8]] + boxplots[[9]] + boxplots[[10]] + boxplots[[11]] + boxplots[[12]]
}


###############################################################################
# OPTIONAL: PAIRWISE DEG BETWEEN GC SUBTYPES (WITHIN ONE CONDITION)           #
#                                                                             #
# Instead of comparing M vs P within one subtype, compare subtypes within     #
# one condition — e.g. what makes GC.1 different from GC.4 in M samples?     #
###############################################################################

# Example: GC.1 vs GC.4 within M_group only
# gc_m <- subset(gc_obj, subset = condition == "M_group" &
#                  SR_superfine_pruned %in% c("GC.1", "GC.4"))
# Idents(gc_m) <- "SR_superfine_pruned"
# gc1_vs_gc4 <- FindMarkers(gc_m, ident.1 = "GC.1", ident.2 = "GC.4",
#                            assay = "SCT", only.pos = FALSE)
# head(gc1_vs_gc4)


###############################################################################
# SAVE GC OBJECT WITH ALL RESULTS                                             #
###############################################################################

saveRDS(gc_obj, file.path(output_dir, "gc_obj_final.rds"))
message("GC object saved to: ", file.path(output_dir, "gc_obj_final.rds"))


################################################################################
# Below are additional DEG comparisons
# This one compares between samples in the sample group

# Run Kruskal-Wallis per gene, per condition group
# Tests: do individual samples within a group differ from each other?

DefaultAssay(gc_obj) <- "Xenium"

kw_within_group <- function(obj, condition_name) {
  # Subset to cells in this condition
  sub_obj <- subset(obj, subset = condition == condition_name)
  
  # Need at least 3 samples to run KW
  n_samples <- length(unique(sub_obj$orig.ident))
  if (n_samples < 3) {
    message("Skipping ", condition_name, " — only ", n_samples, " sample(s)")
    return(NULL)
  }
  
  # Get normalized expression matrix
  expr_mat <- GetAssayData(sub_obj, layer = "data", assay = "Xenium")
  sample_ids <- sub_obj$orig.ident
  
  # Test each gene
  results <- lapply(rownames(expr_mat), function(gene) {
    values <- expr_mat[gene, ]
    # Skip genes with no variance
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

# Run for each condition group
conditions <- unique(gc_obj$condition)
all_kw_results <- list()

for (cond in conditions) {
  message("Running KW for: ", cond)
  res <- kw_within_group(gc_obj, cond)
  if (!is.null(res)) {
    all_kw_results[[cond]] <- res
    sig <- sum(res$p_adj < 0.05)
    message("  Significant genes (FDR < 0.05): ", sig)
  }
}

# View top hits per group
all_kw_combined %>%
  filter(p_adj < 0.05) %>%
  group_by(condition) %>%
  slice_head(n = 10) %>%
  print(n = 40)

####### Visualization #########
library(ggpubr)
library(FSA)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rstatix)

# 1. Set assay
DefaultAssay(gc_obj) <- "Xenium"

# 2. Pick top 3 significant genes for the SEPSIS group
sepsis_condition_name <- "SEPSIS_group" 

top_genes_per_group <- all_kw_combined %>%
  filter(p_adj < 0.05) %>%
  group_by(condition) %>%
  slice_head(n = 3)

sepsis_genes <- top_genes_per_group %>% 
  filter(condition == sepsis_condition_name) %>% 
  pull(gene)

# 3. Helper function to plot Kruskal-Wallis & Dunn's post-hoc
plot_kw_gene <- function(obj, gene, condition_name) {
  sub_obj <- subset(obj, subset = condition == condition_name)
  
  df <- data.frame(
    expr = GetAssayData(sub_obj, layer = "data", assay = "Xenium")[gene, ],
    sample = sub_obj$orig.ident
  )
  
  # Run Dunn's post-hoc test
  dunn_res <- dunnTest(expr ~ sample, data = df, method = "bh")
  dunn_df <- dunn_res$res
  
  # Filter for significant pairwise comparisons
  sig_pairs <- dunn_df %>% filter(P.adj < 0.05)
  
  # Base violin + boxplot
  p <- ggplot(df, aes(x = sample, y = expr, fill = sample)) +
    geom_violin(scale = "width", alpha = 0.7, width = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
    theme_classic() +
    labs(
      title = paste0(gene, " — ", condition_name),
      x = NULL,
      y = "Normalized Expression"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Add significance brackets based on Dunn's ADJUSTED p-values
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
      sig_pairs,
      label = "p.adj.signif",
      y.position = "y.position",
      step.increase = 0.08,
      tip.length = 0.01,
      size = 3.5,
      inherit.aes = FALSE # Prevents inheriting fill = sample from ggplot()
    ) +
      ylim(min(df$expr, na.rm = TRUE), max(sig_pairs$y.position) + y_step)
  }
  
  # Add overall KW adjusted p-value as subtitle
  kw_p <- all_kw_combined %>%
    filter(gene == !!gene, condition == condition_name) %>%
    pull(p_adj)
  
  if (length(kw_p) > 0) {
    p <- p + labs(subtitle = paste0("KW adj. p = ", formatC(kw_p[1], format = "e", digits = 2)))
  }
  
  return(p)
}

# 4. Generate plots for SEPSIS group
sepsis_plots <- lapply(sepsis_genes, function(g) {
  plot_kw_gene(gc_obj, g, sepsis_condition_name)
})

# 5. Display SEPSIS plots
if (length(sepsis_plots) >= 3) {
  sepsis_plots[[1]] + sepsis_plots[[2]] + sepsis_plots[[3]]
} else if (length(sepsis_plots) > 0) {
  wrap_plots(sepsis_plots)
}

####################################
# Oct 2 - scratching the singleR labelling with subtyped GC cells
# Rerun clustering based on gc_sub 

# Use the existing PCA embeddings
gc_obj <- FindNeighbors(gc_obj, dims = 1:15)
gc_obj <- FindClusters(gc_obj, resolution = 0.3)

# Check how many clusters
table(Idents(gc_obj))

# Plot the Seurat clusters
DimPlot(gc_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  coord_fixed() + ggtitle("Seurat Clusters (unsupervised)")

# Compare side by side with SingleR labels
p1 <- DimPlot(gc_obj, reduction = "umap", group.by = "seurat_clusters", 
              label = TRUE) + coord_fixed() + ggtitle("Seurat Clusters")
p2 <- DimPlot(gc_obj, reduction = "umap", group.by = "SR_superfine_pruned", 
              label = TRUE, repel = TRUE) + coord_fixed() + ggtitle("SingleR GC Subtypes")
p1 + p2

# Cross-tabulate to see how they map onto each other
table(gc_obj$seurat_clusters, gc_obj$SR_superfine_pruned)

# Grouped by FOV
DimPlot(gc_obj, reduction = "umap", group.by = "orig.ident",
        pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by FOV")

# Grouped by condition
DimPlot(gc_obj, reduction = "umap", group.by = "condition",
        pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by Condition")


#######################
# This code is used to merge the samples from Mehwish with samples used in
# the epilepsy project

# Subset dat to non-overlapping FOVs
dat_unique <- subset(dat, subset = orig.ident %in% setdiff(unique(dat$orig.ident), overlap))

cols_to_keep <- c("orig.ident", "nCount_Xenium", "nFeature_Xenium",
                  "nCount_BlankCodeword", "nFeature_BlankCodeword",
                  "nCount_ControlCodeword", "nFeature_ControlCodeword",
                  "nCount_ControlProbe", "nFeature_ControlProbe",
                  "SingleR_label_midfine", "SingleR_pruned_midfine",
                  "SingleR_label_superfine", "SingleR_pruned_superfine",
                  "celltype", "cellsubtype", "cellsubtype_fine")

dat_unique@meta.data <- dat_unique@meta.data[, cols_to_keep, drop = FALSE]

colnames(dat_unique@meta.data)

# Merge using merged_obj_clean so your SingleR labels are kept
combined_obj <- merge(merged_obj_clean, dat_unique,
                      add.cell.ids = c("DG", "cortex"))

table(combined_obj$orig.ident)
ncol(combined_obj)

# Save
saveRDS(combined_obj, "X:/CembrowskiLab/Hope Luo/DG analysis/DG_combined_raw_obj.rds")


#############################################
# The following analysis is for creating umap clustering with the combined object
# Combined object contains 10 samples from Mehwish, plus all the samples from LK's 
# Epilepsy paper

####################################
# Subsetting for the correct cells

# Verify they don't overlap
has_cortex <- !is.na(combined_obj$cellsubtype)
has_dg <- !is.na(combined_obj$SR_superfine_pruned)

message("cellsubtype (cortex): ", sum(has_cortex))
message("SR_superfine_pruned (DG): ", sum(has_dg))
message("Both: ", sum(has_cortex & has_dg))
message("Neither: ", sum(!has_cortex & !has_dg))

# If no overlap, merge them
combined_obj$combined_celltype <- NA

# Fill in cortex cells
combined_obj$combined_celltype[has_cortex] <- combined_obj$cellsubtype[has_cortex]

# Fill in DG cells
dg_cells <- !is.na(combined_obj$SR_superfine_pruned) & is.na(combined_obj$combined_celltype)
combined_obj$combined_celltype[dg_cells] <- combined_obj$SR_superfine_pruned[dg_cells]

# Check coverage
message("Total labeled: ", sum(!is.na(combined_obj$combined_celltype)), 
        " / ", ncol(combined_obj))
table(combined_obj$combined_celltype)

gc_labels <- c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5", "ExN_GC")

gc_all <- subset(combined_obj, subset = combined_celltype %in% gc_labels)

table(gc_all$combined_celltype)
table(gc_all$orig.ident)

gc_all$condition <- NA

gc_all$condition[gc_all$orig.ident %in% c("M70", "P28")] <- "MAID"
gc_all$condition[gc_all$orig.ident %in% c("P51", "P60", "P71", "M56")] <- "SEPSIS"
gc_all$condition[gc_all$orig.ident %in% c("M44", "M25", "M41")] <- "HIBI"
gc_all$condition[gc_all$orig.ident %in% c("M38")] <- "TBI"
gc_all$condition[gc_all$orig.ident %in% c("L10")] <- "CARDIAC"
gc_all$condition[gc_all$orig.ident %in% c("L5")] <- "METASTASIS"

# Everything else is EPILEPSY
gc_all$condition[is.na(gc_all$condition)] <- "EPILEPSY"

# Verify
table(gc_all$condition)
table(gc_all$orig.ident, gc_all$condition)

gc_all <- subset(gc_all, subset = orig.ident != "P60")

# Verify it's gone
table(gc_all$orig.ident, gc_all$condition)

gc_all$combined_celltype <- "GC"

# Verify
table(gc_all$combined_celltype)

# Sex
gc_all$sex <- NA
gc_all$sex[gc_all$orig.ident %in% c("B1","D11","D8","I1","C1","C2","E008","D12","D13","D14","D9","E015","M56","M70","P28","P51","P71")] <- "M"
gc_all$sex[gc_all$orig.ident %in% c("F1","H1","J1","L20","L14","L15","L18","L1","L21","L3","L2","L10","L5")] <- "F"

# Patient age
gc_all$patient_age <- NA
gc_all$patient_age[gc_all$orig.ident == "B1"] <- 45
gc_all$patient_age[gc_all$orig.ident == "D11"] <- 46
gc_all$patient_age[gc_all$orig.ident == "D8"] <- 41
gc_all$patient_age[gc_all$orig.ident == "F1"] <- 70
gc_all$patient_age[gc_all$orig.ident == "H1"] <- 74
gc_all$patient_age[gc_all$orig.ident == "I1"] <- 53
gc_all$patient_age[gc_all$orig.ident == "J1"] <- 83
gc_all$patient_age[gc_all$orig.ident == "L20"] <- 52
gc_all$patient_age[gc_all$orig.ident == "C1"] <- 51
gc_all$patient_age[gc_all$orig.ident == "C2"] <- 62
gc_all$patient_age[gc_all$orig.ident == "E008"] <- 44
gc_all$patient_age[gc_all$orig.ident == "L14"] <- 46
gc_all$patient_age[gc_all$orig.ident == "L15"] <- 48
gc_all$patient_age[gc_all$orig.ident == "L18"] <- 50
gc_all$patient_age[gc_all$orig.ident == "D12"] <- 24
gc_all$patient_age[gc_all$orig.ident == "D13"] <- 31
gc_all$patient_age[gc_all$orig.ident == "L1"] <- 27
gc_all$patient_age[gc_all$orig.ident == "L21"] <- 32
gc_all$patient_age[gc_all$orig.ident == "L3"] <- 35
gc_all$patient_age[gc_all$orig.ident == "D14"] <- 32
gc_all$patient_age[gc_all$orig.ident == "D9"] <- 40
gc_all$patient_age[gc_all$orig.ident == "E015"] <- 46
gc_all$patient_age[gc_all$orig.ident == "L2"] <- 44
gc_all$patient_age[gc_all$orig.ident == "L10"] <- 35
gc_all$patient_age[gc_all$orig.ident == "L5"] <- 70
gc_all$patient_age[gc_all$orig.ident == "M56"] <- 68
gc_all$patient_age[gc_all$orig.ident == "M70"] <- 37
gc_all$patient_age[gc_all$orig.ident == "P28"] <- 55
gc_all$patient_age[gc_all$orig.ident == "P51"] <- 66
gc_all$patient_age[gc_all$orig.ident == "P71"] <- 67

# Verify
table(gc_all$sex, gc_all$condition, useNA = "ifany")
table(gc_all$orig.ident, gc_all$patient_age, useNA = "ifany")

# Check which FOVs still have NA
unique(gc_all$orig.ident[is.na(gc_all$sex)])

# Save
saveRDS(gc_all, "X:/CembrowskiLab/Hope Luo/DG analysis/gc_all_combined.rds")

### Check that the subsetting worked:
# Check the combined_celltype column
table(gc_all$combined_celltype, useNA = "ifany")

# Check the fine-level labels
table(gc_all$SR_fine_pruned, useNA = "ifany")

# Check the superfine labels
table(gc_all$SR_superfine_pruned, useNA = "ifany")

# Check the cortex labels if any cortex cells are present
table(gc_all$cellsubtype, useNA = "ifany")

# Count non-GC cells by fine label
non_gc <- gc_all$SR_fine_pruned != "GC" & !is.na(gc_all$SR_fine_pruned)
message("Non-GC cells by fine label: ", sum(non_gc))

######################################

combined_obj <- readRDS("X:/CembrowskiLab/Hope Luo/DG analysis/DG_combined_raw_obj.rds")

##############################################
# Loads the GC cells from all of the combined data
# NOTE: START ANALYSIS HERE

gc_all <- readRDS("X:/CembrowskiLab/Hope Luo/DG analysis/gc_all_combined.rds")

gc_all <- SCTransform(gc_all, assay = "Xenium")
gc_all <- RunPCA(gc_all, features = VariableFeatures(gc_all), npcs = 50, verbose = FALSE)

ElbowPlot(gc_all, ndims = 50)

# Use the PCA embeddings
gc_all <- FindNeighbors(gc_all, dims = 1:15)
gc_all <- FindClusters(gc_all, resolution = 0.2)

# Check how many clusters
table(Idents(gc_all))

set.seed(42)
gc_all <- RunUMAP(gc_all, dims = 1:15)

DimPlot(gc_all, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE) +
  coord_fixed() + ggtitle("Granule Cell Clusters resolution 0.2")

# UMAPS with different labelling

# By FOV
DimPlot(gc_all, reduction = "umap", group.by = "orig.ident", 
        shuffle = TRUE, pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by FOV")

# By Condition
DimPlot(gc_all, reduction = "umap", group.by = "condition", 
        shuffle = TRUE, pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by Condition")

# By Age Group
# Create age groups
gc_all$age_group <- cut(gc_all$patient_age, 
                        breaks = c(0, 30, 40, 50, 60, 70, Inf),
                        labels = c("20-30", "31-40", "41-50", "51-60", "61-70", "70+"),
                        right = TRUE)

# Verify
table(gc_all$age_group, useNA = "ifany")

# Plot
DimPlot(gc_all, reduction = "umap", group.by = "age_group",
        shuffle = TRUE, pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by Age Group")

# By Sex
DimPlot(gc_all, reduction = "umap", group.by = "sex", 
        shuffle = TRUE, pt.size = 0.5) +
  coord_fixed() + ggtitle("GC Cells by Sex")

# Plot feature plot for genes
FeaturePlot(gc_all, features = "PROX1", reduction = "umap",
            min.cutoff = "q10", max.cutoff = "q90")


### Find marker genes for each cluster
# Now prep and find markers

for (model_name in names(gc_all[["SCT"]]@SCTModel.list)) {
  slot(gc_all[["SCT"]]@SCTModel.list[[model_name]], "umi.assay") <- "Xenium"
}

gc_all <- PrepSCTFindMarkers(gc_all)

gc_markers <- FindAllMarkers(gc_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top25 <- gc_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 25) %>%
  ungroup()

top25_list <- split(top25$gene, top25$cluster)
top25_list

gc_all[["Xenium"]] <- JoinLayers(gc_all[["Xenium"]])

# Verify it's now a single layer
Layers(gc_all[["Xenium"]])

# Re-normalize
gc_all <- NormalizeData(gc_all, assay = "Xenium")


# Plot feature plot for specified gene
FeaturePlot(gc_all, features = "PROX1", reduction = "umap",
            min.cutoff = "q10", max.cutoff = "q90")

# Visualziation of hyper vs hypo excitable markers
DefaultAssay(gc_all) <- "Xenium"

# Hypoexcitable markers
FeaturePlot(gc_all, features = c("SLC12A5", "CALB1", "GABRA5"),
            ncol = 3, cols = c("#F0F0F0", "#FF4500")) +
  plot_annotation(title = "Hypoexcitable markers")

# Hyperexcitable markers
FeaturePlot(gc_all, features = c("HOMER1", "SNAP25", "KCNMA1"),
            ncol = 3, cols = c("#F0F0F0", "#FF4500")) +
  plot_annotation(title = "Hyperexcitable markers")


# Bubble plot comparing the expression of hyper versus hypo excitable genes
# for each cluster
DotPlot(gc_all, features = c("SLC12A5", "CALB1", "GABRA5", 
                             "HOMER1", "SNAP25", "KCNMA1"),
        group.by = "seurat_clusters", cols = c("blue", "red")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Hypo vs Hyper Excitability Markers by Cluster")

# Module score approach
hypo_genes <- list(c("SLC12A5", "CALB1", "GABRA5", "GRIN2A", "GABRB3", "SEMA5A"))
hyper_genes <- list(c("HOMER1", "SNAP25", "KCNMA1", "DLG1", "DLG2", "GABRA2"))

gc_all <- AddModuleScore(gc_all, features = hypo_genes, name = "hypo_score", 
                         assay = "Xenium", ctrl = 5, nbin = 10)
gc_all <- AddModuleScore(gc_all, features = hyper_genes, name = "hyper_score", 
                         assay = "Xenium", ctrl = 5, nbin = 10)
# Plot the scores
FeaturePlot(gc_all, features = "hypo_score1", 
            cols = c("#F0F0F0", "blue")) +
  coord_fixed() + ggtitle("Hypoexcitability Score")

FeaturePlot(gc_all, features = "hyper_score1",
            cols = c("#F0F0F0", "red")) +
  coord_fixed() + ggtitle("Hyperexcitability Score")

# Plot the difference (hyper - hypo)
gc_all$excitability_balance <- gc_all$hyper_score1 - gc_all$hypo_score1

FeaturePlot(gc_all, features = "excitability_balance",
            cols = c("blue", "#F0F0F0", "red")) +
  coord_fixed() + ggtitle("Excitability Balance (red = hyper, blue = hypo)")

# Compare the scores calculated above by condition
VlnPlot(gc_all, features = "excitability_balance", 
        group.by = "condition", pt.size = 0) +  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggtitle("Excitability Balance by Condition")

p1 <- VlnPlot(gc_all, features = "hypo_score1", group.by = "condition", pt.size = 0) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3))

p2 <- VlnPlot(gc_all, features = "hyper_score1", group.by = "condition", pt.size = 0) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3))

p1 + p2

# Box plot of excitability balance by condition — this is the key plot
ggplot(gc_all@meta.data, aes(x = condition, y = excitability_balance, fill = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic() +
  labs(x = NULL, y = "Excitability Balance\n(+ hyper / - hypo)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Excitability Balance by Condition")

# Also check per sample to see within-group variability
ggplot(gc_all@meta.data, aes(x = orig.ident, y = excitability_balance, fill = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_classic() +
  labs(x = NULL, y = "Excitability Balance\n(+ hyper / - hypo)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Excitability Balance by Sample")

# Checking nCount and nFeature against the excitability gradient
# 1. Is it driven by QC metrics (more transcripts = more genes detected)?
p1 <- FeaturePlot(gc_all, features = "nCount_Xenium", cols = c("#F0F0F0", "darkblue")) +
  coord_fixed() + ggtitle("Total Transcripts")
p2 <- FeaturePlot(gc_all, features = "nFeature_Xenium", cols = c("#F0F0F0", "darkblue")) +
  coord_fixed() + ggtitle("Genes Detected")
p3 <- FeaturePlot(gc_all, features = "excitability_balance", cols = c("blue", "#F0F0F0", "red")) +
  coord_fixed() + ggtitle("Excitability Balance")
p1 + p2 + p3

# More checks to figure out why there is a gradient in the UMAP
# 2. Is it driven by specific samples?
DimPlot(gc_all, reduction = "umap", group.by = "orig.ident", shuffle = TRUE) +
  coord_fixed() + ggtitle("By Sample")

# 3. Check correlation between excitability and QC
cor.test(gc_all$excitability_balance, gc_all$nCount_Xenium, method = "spearman")
cor.test(gc_all$excitability_balance, gc_all$nFeature_Xenium, method = "spearman")

# 4. Is it the clusters?
ggplot(gc_all@meta.data, aes(x = seurat_clusters, y = excitability_balance)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ggtitle("Excitability Balance by Cluster")

# 5. Check individual hypo and hyper scores separately
p4 <- FeaturePlot(gc_all, features = "hypo_score1", cols = c("#F0F0F0", "blue")) +
  coord_fixed() + ggtitle("Hypo Score")
p5 <- FeaturePlot(gc_all, features = "hyper_score1", cols = c("#F0F0F0", "red")) +
  coord_fixed() + ggtitle("Hyper Score")
p4 + p5

# 6. Check key individual genes that define the gradient
FeaturePlot(gc_all, features = c("HOMER1", "SLC12A5", "SNAP25", "CALB1"),
            ncol = 2, cols = c("#F0F0F0", "darkblue"))

# 7. Is it age or sex?
p6 <- DimPlot(gc_all, reduction = "umap", group.by = "sex", shuffle = TRUE) +
  coord_fixed() + ggtitle("By Sex")
p7 <- DimPlot(gc_all, reduction = "umap", group.by = "age_group", shuffle = TRUE) +
  coord_fixed() + ggtitle("By Age Group")
p6 + p7

# Visual: does nCount follow the same gradient?
p1 <- FeaturePlot(gc_all, features = "nCount_Xenium", cols = c("#F0F0F0", "darkblue")) +
  coord_fixed() + ggtitle("Total Transcripts")
p2 <- FeaturePlot(gc_all, features = "nFeature_Xenium", cols = c("#F0F0F0", "darkblue")) +
  coord_fixed() + ggtitle("Genes Detected")
p3 <- FeaturePlot(gc_all, features = "excitability_balance", cols = c("blue", "#F0F0F0", "red")) +
  coord_fixed() + ggtitle("Excitability Balance")
p1 + p2 + p3

# Check ILAE score — epilepsy severity might drive the gradient
# First add ILAE if not already there
if (!"ILAE_score" %in% colnames(gc_all@meta.data)) {
  gc_all$ILAE_score <- NA
  gc_all$ILAE_score[gc_all$orig.ident %in% c("B1","D11","D8","F1","H1","I1","J1","L20")] <- "0"
  gc_all$ILAE_score[gc_all$orig.ident %in% c("C1","C2","E008","L14","L15","L18")] <- "1"
  gc_all$ILAE_score[gc_all$orig.ident %in% c("D12","D13","L1","L21","L3")] <- "2"
  gc_all$ILAE_score[gc_all$orig.ident %in% c("D14","D9","E015","L2")] <- "3"
  gc_all$ILAE_score[gc_all$condition != "EPILEPSY"] <- "PM"
}

DimPlot(gc_all, reduction = "umap", group.by = "ILAE_score", shuffle = TRUE) +
  coord_fixed() + ggtitle("By ILAE Score")

# Quantify: does ILAE score correlate with excitability?
gc_epi <- subset(gc_all, subset = condition == "EPILEPSY")
gc_epi$ILAE_numeric <- as.numeric(gc_epi$ILAE_score)

cor.test(gc_epi$excitability_balance, gc_epi$ILAE_numeric, method = "spearman")

ggplot(gc_epi@meta.data, aes(x = ILAE_score, y = excitability_balance, fill = ILAE_score)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ggtitle("Excitability Balance by ILAE Score (Epilepsy Only)")

# 6. Check key individual genes that define the gradient
FeaturePlot(gc_all, features = c("HOMER1", "SLC12A5", "SNAP25", "CALB1"),
            ncol = 2, cols = c("#F0F0F0", "darkblue"))

# Further exploration with homer1 and calb1 from trends seen in #6 above
# Check HOMER1 and CALB1 by condition — within postmortem only
gc_pm <- subset(gc_all, subset = tissue_type == "postmortem")
DefaultAssay(gc_pm) <- "Xenium"

p1 <- VlnPlot(gc_pm, features = "HOMER1", group.by = "condition", pt.size = 0) +
  ggtitle("HOMER1 — Postmortem Only")
p2 <- VlnPlot(gc_pm, features = "CALB1", group.by = "condition", pt.size = 0) +
  ggtitle("CALB1 — Postmortem Only")
p1 + p2

# Redoing the excitabiluty gradient based on homer1 and calb1
homer1_expr <- GetAssayData(gc_all, layer = "data", assay = "Xenium")["HOMER1", ]
calb1_expr <- GetAssayData(gc_all, layer = "data", assay = "Xenium")["CALB1", ]

gc_all$homer1_calb1_ratio <- homer1_expr - calb1_expr

FeaturePlot(gc_all, features = "homer1_calb1_ratio",
            cols = c("blue", "#F0F0F0", "red")) +
  coord_fixed() + ggtitle("HOMER1 - CALB1 (red = active, blue = quiescent)")

# Figure out which genes are diving the gradient
# Get UMAP coordinates
umap_coords <- Embeddings(gc_all, reduction = "umap")
umap_x <- umap_coords[, 1]

# Get expression matrix
DefaultAssay(gc_all) <- "Xenium"
expr_mat <- GetAssayData(gc_all, layer = "data", assay = "Xenium")

# Correlate every gene with UMAP x-position
gene_cors <- sapply(rownames(expr_mat), function(gene) {
  cor(expr_mat[gene, ], umap_x, method = "spearman")
})

# Sort by correlation
gene_cors_df <- data.frame(
  gene = names(gene_cors),
  rho = gene_cors
) %>% arrange(desc(abs(rho)))

# Top genes driving the gradient
message("Top 20 genes driving left-right gradient:")
head(gene_cors_df, 20)

# Positive rho = higher on the right (hyperexcitable side)
message("\nTop 10 RIGHT-enriched (active/hyper):")
gene_cors_df %>% filter(rho > 0) %>% head(10)

# Negative rho = higher on the left (quiescent side)
message("\nTop 10 LEFT-enriched (quiescent/hypo):")
gene_cors_df %>% filter(rho < 0) %>% head(10)

# Plot top 3 from each direction
right_genes <- gene_cors_df %>% filter(rho > 0) %>% head(3) %>% pull(gene)
left_genes <- gene_cors_df %>% filter(rho < 0) %>% head(3) %>% pull(gene)

FeaturePlot(gc_all, features = right_genes, ncol = 3, 
            cols = c("#F0F0F0", "red")) +
  plot_annotation(title = "Top RIGHT-side genes (active)")

FeaturePlot(gc_all, features = left_genes, ncol = 3,
            cols = c("#F0F0F0", "blue")) +
  plot_annotation(title = "Top LEFT-side genes (quiescent)")

write.csv(gene_cors_df, file.path(output_dir, "UMAP_gradient_gene_correlations.csv"),
          row.names = FALSE)

# More plotting
# Update the ratio to reflect the real axis
gc_all$maturity_score <- calb1_expr - GetAssayData(gc_all, layer = "data", assay = "Xenium")["GRIN3A", ]

FeaturePlot(gc_all, features = "maturity_score",
            cols = c("red", "#F0F0F0", "blue")) +
  coord_fixed() + ggtitle("Maturity Score (blue = mature, red = plastic)")

# Compare across conditions — within postmortem only
gc_pm <- subset(gc_all, subset = tissue_type == "postmortem")

ggplot(gc_pm@meta.data, aes(x = condition, y = maturity_score, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Maturity Score — Postmortem Only")

# Statsical testing for the above
gc_pm@meta.data %>%
  group_by(condition) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(orig.ident),
    mean_maturity = round(mean(maturity_score), 3),
    median_maturity = round(median(maturity_score), 3)
  ) %>%
  arrange(mean_maturity)

# Statistical test
kruskal.test(maturity_score ~ condition, data = gc_pm@meta.data)

# Pairwise: which conditions differ?
pairwise.wilcox.test(gc_pm@meta.data$maturity_score, 
                     gc_pm@meta.data$condition, 
                     p.adjust.method = "BH")

# plot each condition by maturity score
conditions <- unique(gc_all$condition)

plots <- lapply(conditions, function(cond) {
  sub <- subset(gc_all, subset = condition == cond)
  FeaturePlot(sub, features = "maturity_score",
              cols = c("red", "#F0F0F0", "blue")) +
    ggtitle(paste0(cond, " (n=", ncol(sub), ")")) +
    NoLegend()
})

wrap_plots(plots, ncol = 4)

plots_homer <- lapply(conditions, function(cond) {
  sub <- subset(gc_all, subset = condition == cond)
  FeaturePlot(sub, features = "HOMER1",
              cols = c("#F0F0F0", "red")) +
    ggtitle(cond) + NoLegend()
})

plots_calb <- lapply(conditions, function(cond) {
  sub <- subset(gc_all, subset = condition == cond)
  FeaturePlot(sub, features = "CALB1",
              cols = c("#F0F0F0", "blue")) +
    ggtitle(cond) + NoLegend()
})

wrap_plots(plots_homer, ncol = 4) + plot_annotation(title = "HOMER1 by Condition")
wrap_plots(plots_calb, ncol = 4) + plot_annotation(title = "CALB1 by Condition")

# examine spatial localization of specific genes 
# For each HIBI FOV, do the immature/plastic cells 
# sit in the SGZ (inner GCL) vs outer GCL?
# This is something only spatial data can answer
# HIBI samples
p1 <- ImageFeaturePlot(gc_all, fov = "M44", features = "maturity_score", 
                       cols = c("red", "#F0F0F0", "blue")) + ggtitle("M44 (HIBI)")
p2 <- ImageFeaturePlot(gc_all, fov = "M25", features = "maturity_score", 
                       cols = c("red", "#F0F0F0", "blue")) + ggtitle("M25 (HIBI)")
p3 <- ImageFeaturePlot(gc_all, fov = "M41", features = "maturity_score", 
                       cols = c("red", "#F0F0F0", "blue")) + ggtitle("M41 (HIBI)")

p4 <- ImageFeaturePlot(gc_all, fov = "L10", features = "maturity_score", 
                       cols = c("red", "#F0F0F0", "blue")) + ggtitle("L10 (CARDIAC)")

p5 <- ImageFeaturePlot(gc_all, fov = "B1", features = "maturity_score", 
                       cols = c("red", "#F0F0F0", "blue")) + ggtitle("B1 (Epilepsy ILAE 0)")
p6 <- ImageFeaturePlot(gc_all, fov = "D9", features = "maturity_score", 
                       cols = c("red", "#F0F0F0", "blue")) + ggtitle("D9 (Epilepsy ILAE 3)")
print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)
dev.off()

# Quantify the spatial pattern
# Run the spatial correlation for every FOV
spatial_cors <- data.frame()

for (sample_id in unique(gc_all$orig.ident)) {
  sub <- subset(gc_all, subset = orig.ident == sample_id)
  
  fov <- Images(sub)
  if (length(fov) == 0) next
  
  tryCatch({
    coords <- GetTissueCoordinates(sub, image = fov[1])
    coords$maturity_score <- sub$maturity_score
    
    centroid_x <- mean(coords[,1])
    centroid_y <- mean(coords[,2])
    coords$dist <- sqrt((coords[,1] - centroid_x)^2 + 
                          (coords[,2] - centroid_y)^2)
    
    ct <- cor.test(coords$dist, coords$maturity_score, method = "spearman")
    
    spatial_cors <- rbind(spatial_cors, data.frame(
      sample = sample_id,
      condition = unique(sub$condition),
      rho = ct$estimate,
      p_value = ct$p.value,
      n_cells = ncol(sub)
    ))
  }, error = function(e) {
    message("Skipped ", sample_id, ": ", e$message)
  })
}

print(spatial_cors)

# Plot the spatial correlation by condition
ggplot(spatial_cors, aes(x = condition, y = rho, color = condition)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Spatial correlation (rho)\n(+ = mature cells outer, plastic cells inner)",
       x = NULL) +
  ggtitle("Spatial Organization of Maturity Gradient by Condition")

ggplot(spatial_cors, aes(x = reorder(sample, rho), y = rho, fill = condition)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Spatial correlation (rho)") +
  ggtitle("Maturity vs Distance from Center")

# Summarize by condition
spatial_cors %>%
  group_by(condition) %>%
  summarise(
    n_samples = n(),
    mean_rho = round(mean(rho), 3),
    sd_rho = round(sd(rho), 3)
  ) %>%
  arrange(mean_rho)

# Oct 15: realized that M44 might also be a bad sample
# Will need to subset out

ImageDimPlot(gc_all, fov = "M44", group.by = "seurat_clusters",
             size = 0.5, border.size = NA)

gc_all <- subset(gc_all, subset = orig.ident != "M44")

# Verify it's gone
table(gc_all$orig.ident)

# Remove only M25 and L5
cells_to_keep <- colnames(gc_all)[!(gc_all$orig.ident %in% c("M25", "L5"))]
gc_all <- subset(gc_all, cells = cells_to_keep)


saveRDS(gc_all, "X:/CembrowskiLab/Hope Luo/DG analysis/gc_all_combined_updated.rds")

### CONTINUE ANALYSIS HERE

gc_all <- RunPCA(gc_all, features = VariableFeatures(gc_all), npcs = 50, verbose = FALSE)

gc_all <- FindNeighbors(gc_all, dims = 1:15)
gc_all <- FindClusters(gc_all, resolution = 0.2)

# Check how many clusters
table(Idents(gc_all))

set.seed(42)
gc_all <- RunUMAP(gc_all, dims = 1:15)

### PROX1 analysis
# Plot feature plot for PROX1 to see if there is visible enrichment
FeaturePlot(gc_all, features = "PROX1", reduction = "umap",
            min.cutoff = "q10", max.cutoff = "q90")

DefaultAssay(gc_all) <- "Xenium"

# HIBI
p1 <- ImageFeaturePlot(gc_all, fov = "M25", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("M25 (HIBI)")
p2 <- ImageFeaturePlot(gc_all, fov = "M41", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("M41 (HIBI)")
p1 + p2

# SEPSIS
p3 <- ImageFeaturePlot(gc_all, fov = "M56", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("M56 (SEPSIS)")
p4 <- ImageFeaturePlot(gc_all, fov = "P51", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("P51 (SEPSIS)")
p5 <- ImageFeaturePlot(gc_all, fov = "P71", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("P71 (SEPSIS)")
p3 + p4 + p5

# MAID
p6 <- ImageFeaturePlot(gc_all, fov = "M70", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("M70 (MAID)")
p7 <- ImageFeaturePlot(gc_all, fov = "P28", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("P28 (MAID)")
p6 + p7

# TBI + CARDIAC + METASTASIS
p8 <- ImageFeaturePlot(gc_all, fov = "M38", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("M38 (TBI)")
p9 <- ImageFeaturePlot(gc_all, fov = "L10", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("L10 (CARDIAC)")
p10 <- ImageFeaturePlot(gc_all, fov = "L5", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("L5 (METASTASIS)")
p8 + p9 + p10

# Epilepsy - ILAE 0
p11 <- ImageFeaturePlot(gc_all, fov = "B1", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("B1 (ILAE 0)")
p12 <- ImageFeaturePlot(gc_all, fov = "F1", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("F1 (ILAE 0)")
p13 <- ImageFeaturePlot(gc_all, fov = "D9", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("D9 (ILAE 3)")
p14 <- ImageFeaturePlot(gc_all, fov = "E015", features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.5) + ggtitle("E015 (ILAE 3)")
p11 + p12 + p13 + p14


### Testing the shift in maturity between different conditions
# Regress out nCount before comparing

DefaultAssay(gc_all) <- "Xenium"
calb1_expr <- GetAssayData(gc_all, layer = "data", assay = "Xenium")["CALB1", ]
grin3a_expr <- GetAssayData(gc_all, layer = "data", assay = "Xenium")["GRIN3A", ]

gc_all$maturity_score <- calb1_expr - grin3a_expr

# Now run the model
model <- lm(maturity_score ~ nCount_Xenium + condition, data = gc_all@meta.data)
summary(model)

# Linear model: maturity ~ nCount + condition
model <- lm(maturity_score ~ nCount_Xenium + condition, data = gc_all@meta.data)
summary(model)

# The condition coefficients tell you the effect AFTER controlling for nCount
# Extract just the condition effects
library(broom)
tidy(model) %>% filter(grepl("condition", term)) %>% print()

# Or store the nCount-corrected maturity score
gc_all$maturity_corrected <- residuals(lm(maturity_score ~ nCount_Xenium, 
                                          data = gc_all@meta.data))

ggplot(gc_all@meta.data, aes(x = condition, y = maturity_corrected, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Maturity Score (corrected for nCount)")

# Per sample view
ggplot(gc_all@meta.data, aes(x = orig.ident, y = maturity_corrected, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Maturity Score (corrected for nCount)")

# Label the above plot by ILAE grades
# First make sure ILAE_score exists
if (!"ILAE_score" %in% colnames(gc_all@meta.data)) {
  gc_all$ILAE_score <- NA
  gc_all$ILAE_score[gc_all$orig.ident %in% c("B1","D11","D8","F1","H1","I1","J1","L20")] <- "0"
  gc_all$ILAE_score[gc_all$orig.ident %in% c("C1","C2","E008","L14","L15","L18")] <- "1"
  gc_all$ILAE_score[gc_all$orig.ident %in% c("D12","D13","L1","L21","L3")] <- "2"
  gc_all$ILAE_score[gc_all$orig.ident %in% c("D14","D9","E015","L2")] <- "3"
  gc_all$ILAE_score[gc_all$condition != "EPILEPSY"] <- "PM"
}

# Create a label combining condition + ILAE for ordering
gc_all$sample_order <- paste0(gc_all$condition, "_", gc_all$ILAE_score, "_", gc_all$orig.ident)

# Define sample order: epilepsy by ILAE, then postmortem by condition
ilae0 <- c("B1","D11","D8","F1","H1","I1","J1","L20")
ilae1 <- c("C1","C2","E008","L14","L15","L18")
ilae2 <- c("D12","D13","L1","L21","L3")
ilae3 <- c("D14","D9","E015","L2")
pm_samples <- c("L10","M41","M70","P28","M56","P51","P71","M38")

sample_order <- c(ilae0, ilae1, ilae2, ilae3, pm_samples)
# Keep only samples that exist in the data
sample_order <- sample_order[sample_order %in% unique(gc_all$orig.ident)]

gc_all$orig.ident <- factor(gc_all$orig.ident, levels = sample_order)

# Plot
p <- ggplot(gc_all@meta.data, aes(x = orig.ident, y = maturity_corrected, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(80, 10, 10, 10)) +
  ggtitle("Maturity Score (corrected for nCount)")

# Add ILAE brackets using annotate
y_bracket <- max(gc_all$maturity_corrected, na.rm = TRUE) * 1.15

# Get x positions (numeric positions of samples)
ilae0_pos <- which(sample_order %in% ilae0)
ilae1_pos <- which(sample_order %in% ilae1)
ilae2_pos <- which(sample_order %in% ilae2)
ilae3_pos <- which(sample_order %in% ilae3)
pm_pos <- which(sample_order %in% pm_samples)

p + 
  annotate("segment", x = min(ilae0_pos), xend = max(ilae0_pos), 
           y = y_bracket, yend = y_bracket) +
  annotate("text", x = mean(ilae0_pos), y = y_bracket * 1.05, 
           label = "ILAE 0", size = 3) +
  annotate("segment", x = min(ilae1_pos), xend = max(ilae1_pos), 
           y = y_bracket, yend = y_bracket) +
  annotate("text", x = mean(ilae1_pos), y = y_bracket * 1.05, 
           label = "ILAE 1", size = 3) +
  annotate("segment", x = min(ilae2_pos), xend = max(ilae2_pos), 
           y = y_bracket, yend = y_bracket) +
  annotate("text", x = mean(ilae2_pos), y = y_bracket * 1.05, 
           label = "ILAE 2", size = 3) +
  annotate("segment", x = min(ilae3_pos), xend = max(ilae3_pos), 
           y = y_bracket, yend = y_bracket) +
  annotate("text", x = mean(ilae3_pos), y = y_bracket * 1.05, 
           label = "ILAE 3", size = 3) +
  annotate("segment", x = min(pm_pos), xend = max(pm_pos), 
           y = y_bracket, yend = y_bracket) +
  annotate("text", x = mean(pm_pos), y = y_bracket * 1.05, 
           label = "Postmortem", size = 3) +
  coord_cartesian(clip = "off")

# stat tests
# Kruskal-Wallis on corrected scores
kruskal.test(maturity_corrected ~ condition, data = gc_all@meta.data)

# Pairwise within postmortem
pm_data <- gc_all@meta.data %>% filter(condition != "EPILEPSY")
pairwise.wilcox.test(pm_data$maturity_corrected, pm_data$condition, p.adjust.method = "BH")

############################################################
# Checking PROX1 expression on the entire dataset (with 10 samples)
# Load merged object
merged_obj <- readRDS("X:/CembrowskiLab/Hope Luo/DG analysis/DG_merged_raw_obj.rds")

# QC filter
merged_obj <- subset(merged_obj, subset = nFeature_Xenium > 5 & nCount_Xenium > 0)

# SCTransform
merged_obj <- SCTransform(merged_obj, assay = "Xenium")

# PCA
merged_obj <- RunPCA(merged_obj, features = VariableFeatures(merged_obj), npcs = 50, verbose = FALSE)

# Clustering
merged_obj <- FindNeighbors(merged_obj, dims = 1:15)
merged_obj <- FindClusters(merged_obj, resolution = 0.2)

# UMAP
set.seed(42)
merged_obj <- RunUMAP(merged_obj, dims = 1:15)

# Join layers and normalize Xenium for plotting
merged_obj[["Xenium"]] <- JoinLayers(merged_obj[["Xenium"]])
merged_obj <- NormalizeData(merged_obj, assay = "Xenium")
DefaultAssay(merged_obj) <- "Xenium"

set.seed(42)
merged_obj <- RunUMAP(merged_obj, dims = 1:15)

FeaturePlot(merged_obj, features = "PROX1", reduction = "umap",
            min.cutoff = "q10", max.cutoff = "q90")
DefaultAssay(merged_obj) <- "Xenium"

# HIBI
p1 <- ImageFeaturePlot(merged_obj, fov = "M25", features = "PROX1", cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("M25 (HIBI)")
p2 <- ImageFeaturePlot(merged_obj, fov = "M41", features = "PROX1", cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("M41 (HIBI)")
p3 <- ImageFeaturePlot(merged_obj, fov = "M44", features = "PROX1", cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("M44 (HIBI)")
p1 + p2 + p3

# SEPSIS
p4 <- ImageFeaturePlot(merged_obj, fov = "M56", features = "PROX1", cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("M56 (SEPSIS)")
p5 <- ImageFeaturePlot(merged_obj, fov = "P51", features = "PROX1", cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("P51 (SEPSIS)")
p6 <- ImageFeaturePlot(merged_obj, fov = "P60", features = "PROX1", cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("P60 (SEPSIS)")
p7 <- ImageFeaturePlot(merged_obj, fov = "P71", features = "PROX1", cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("P71 (SEPSIS)")
p4 + p5 + p6 + p7

# MAID
p8 <- ImageFeaturePlot(merged_obj, fov = "M70", features = "PROX1", cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("M70 (MAID)")
p9 <- ImageFeaturePlot(merged_obj, fov = "P28", features = "PROX1", cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("P28 (MAID)")
p8 + p9

# TBI
p10 <- ImageFeaturePlot(merged_obj, fov = "M38", features = "PROX1", cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("M38 (TBI)")
p10

### Testing PROX1 for the epilepsy dataset as well
dat <- readRDS("X:/CembrowskiLab/Hope Luo/Analysis/MegaData_v8_annotated_smaller_trimmed.rds")

# Join layers and normalize — only steps needed for spatial plotting
dat[["Xenium"]] <- JoinLayers(dat[["Xenium"]])
dat <- NormalizeData(dat, assay = "Xenium")
DefaultAssay(dat) <- "Xenium"

# Check available FOVs
Images(dat)

# Plot 4 at a time — adjust FOV names based on what Images() returns
# Batch 1
p1 <- ImageFeaturePlot(dat, fov = Images(dat)[1], features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.3) + ggtitle(Images(dat)[1])
p2 <- ImageFeaturePlot(dat, fov = Images(dat)[2], features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.3) + ggtitle(Images(dat)[2])
p3 <- ImageFeaturePlot(dat, fov = Images(dat)[3], features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.3) + ggtitle(Images(dat)[3])
p4 <- ImageFeaturePlot(dat, fov = Images(dat)[4], features = "PROX1", cols = c("#F0F0F0", "darkblue"), size = 0.3) + ggtitle(Images(dat)[4])
p1 + p2 + p3 + p4

######################################################

### DEG analysis for gc_all
# Set up
Idents(gc_all) <- "seurat_clusters"
DefaultAssay(gc_all) <- "SCT"

# Fix SCT models
for (i in seq_along(gc_all[["SCT"]]@SCTModel.list)) {
  slot(gc_all[["SCT"]]@SCTModel.list[[i]], "umi.assay") <- "Xenium"
}
gc_all <- PrepSCTFindMarkers(gc_all)

# Find markers
gc_markers <- FindAllMarkers(gc_all, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

# How many per cluster
table(gc_markers$cluster)

# Top 25 per cluster
gc_top25 <- gc_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 25)

gc_top25_list <- split(gc_top25$gene, gc_top25$cluster)
gc_top25_list

# Oct 19:
# Plot low expressing PROX1 GCs versus mean nCount
DefaultAssay(gc_all) <- "Xenium"
prox1_expr <- GetAssayData(gc_all, layer = "data", assay = "Xenium")["PROX1", ]

# Calculate per-FOV stats
fov_stats <- gc_all@meta.data %>%
  mutate(prox1_zero = prox1_expr[colnames(gc_all)] == 0) %>%
  group_by(orig.ident, condition) %>%
  summarise(
    n_cells = n(),
    mean_nCount = mean(nCount_Xenium),
    pct_zero_prox1 = mean(prox1_zero) * 100,
    .groups = "drop"
  )

print(fov_stats)

# Correlation test
cor_result <- cor.test(fov_stats$mean_nCount, fov_stats$pct_zero_prox1, method = "spearman")
message("Spearman rho = ", round(cor_result$estimate, 3), ", p = ", signif(cor_result$p.value, 3))

# Plot
ggplot(fov_stats, aes(x = mean_nCount, y = pct_zero_prox1, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text(aes(label = orig.ident), vjust = -0.8, size = 2.5) +
  theme_classic() +
  labs(x = "Mean nCount per FOV",
       y = "% GC Cells with Zero PROX1",
       color = "Condition") +
  annotate("text", x = max(fov_stats$mean_nCount) * 0.7, 
           y = max(fov_stats$pct_zero_prox1) * 0.9,
           label = paste0("rho = ", round(cor_result$estimate, 3), 
                          "\np = ", signif(cor_result$p.value, 3)),
           size = 4) +
  ggtitle("PROX1 Dropout vs Transcript Capture per FOV")

# Simple test plot first
plot(fov_stats$mean_nCount, fov_stats$pct_zero_prox1, 
     main = "PROX1 Dropout vs nCount",
     xlab = "Mean nCount", ylab = "% Zero PROX1")

p <- ggplot(fov_stats, aes(x = mean_nCount, y = pct_zero_prox1, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text(aes(label = orig.ident), vjust = -0.8, size = 2.5) +
  theme_classic() +
  labs(x = "Mean nCount per FOV", y = "% GC Cells with Zero PROX1") +
  ggtitle("PROX1 Dropout vs Transcript Capture per FOV")

print(p)

# Check if other housekeeping genes show the same dropout pattern
# If ALL genes drop equally with nCount, it's technical
# If PROX1 drops MORE than expected, there's a biological component

genes_to_compare <- c("PROX1", "SNAP25", "SLC17A7", "GRIA1", "CALB1")

dropout_by_gene <- sapply(genes_to_compare, function(g) {
  expr <- GetAssayData(gc_all, layer = "data", assay = "Xenium")[g, ]
  tapply(expr == 0, gc_all$orig.ident, mean) * 100
})

# Correlate each gene's dropout with nCount
ncount_per_fov <- tapply(gc_all$nCount_Xenium, gc_all$orig.ident, mean)

sapply(genes_to_compare, function(g) {
  round(cor(ncount_per_fov, dropout_by_gene[, g], method = "spearman"), 3)
})

### Plotting the principal components for gc_all clutsering (resolution = 0.2, dims = 1:15)
VizDimLoadings(gc_all, dims = 1, reduction = "pca")
VizDimLoadings(gc_all, dims = 2, reduction = "pca")
VizDimLoadings(gc_all, dims = 3, reduction = "pca")

### Plot feature plot for ion channel genes
DefaultAssay(gc_all) <- "Xenium"

# KCNA1 encodes for Kv1.1
FeaturePlot(gc_all, features = "KCNA1", raster = TRUE,
            min.cutoff = "q10", max.cutoff = "q90")

### ANALYSIS OF POSTN

DefaultAssay(gc_all) <- "Xenium"

# 1. Feature plot on UMAP
FeaturePlot(gc_all, features = "POSTN", raster = TRUE,
            min.cutoff = "q10", max.cutoff = "q90") +
  ggtitle("POSTN Expression")

# 2. Which FOVs have the highest POSTN expression?
postn_expr <- GetAssayData(gc_all, layer = "data", assay = "Xenium")["POSTN", ]

postn_by_fov <- gc_all@meta.data %>%
  mutate(postn = postn_expr) %>%
  group_by(orig.ident, condition) %>%
  summarise(
    mean_postn = round(mean(postn), 3),
    pct_expressing = round(mean(postn > 0) * 100, 1),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_postn))
print(postn_by_fov)

# Percent expressing per FOV
ggplot(postn_by_fov, aes(x = reorder(orig.ident, -pct_expressing), y = pct_expressing, fill = condition)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "% Cells Expressing POSTN") +
  ggtitle("POSTN Detection Rate by FOV")

# 3. Expression per cluster — is the difference significant?
VlnPlot(gc_all, features = "POSTN", group.by = "seurat_clusters", pt.size = 0) +
  ggtitle("POSTN by Cluster")

# Statistical test
kruskal.test(postn_expr ~ gc_all$seurat_clusters)

# Pairwise comparisons
pairwise.wilcox.test(postn_expr, gc_all$seurat_clusters, p.adjust.method = "BH")

# Dot plot across clusters
DotPlot(gc_all, features = "POSTN", group.by = "seurat_clusters") +
  ggtitle("POSTN — % Expressing and Average Expression by Cluster")

# 4. Spatial plots — grouped by condition
# NOTE: the spatial plots are really bad right now - need to be fixed before they can be used for anything
# HIBI
ImageFeaturePlot(gc_all, fov = "M41", features = "POSTN", 
                       cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("M41 (HIBI)")

# SEPSIS
p2 <- ImageFeaturePlot(gc_all, fov = "M56", features = "POSTN",
                       cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("M56 (SEPSIS)")
p3 <- ImageFeaturePlot(gc_all, fov = "P51", features = "POSTN",
                       cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("P51 (SEPSIS)")
p2 + p3

# MAID
p4 <- ImageFeaturePlot(gc_all, fov = "M70", features = "POSTN",
                       cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("M70 (MAID)")
p5 <- ImageFeaturePlot(gc_all, fov = "P28", features = "POSTN",
                       cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("P28 (MAID)")
p4 + p5

# TBI
p6 <- ImageFeaturePlot(gc_all, fov = "M38", features = "POSTN",
                       cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("M38 (TBI)")
p6

# Epilepsy — a few representative
p7 <- ImageFeaturePlot(gc_all, fov = "B1", features = "POSTN",
                       cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("B1 (ILAE 0)")
p8 <- ImageFeaturePlot(gc_all, fov = "D9", features = "POSTN",
                       cols = c("white", "red"), size = 0.5, dark.background = FALSE) + ggtitle("D9 (ILAE 3)")
p7 + p8

# 5. Bar plot of mean POSTN per FOV colored by condition
ggplot(postn_by_fov, aes(x = reorder(orig.ident, -mean_postn), y = mean_postn, fill = condition)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Mean POSTN Expression") +
  ggtitle("Mean POSTN Expression by FOV")

# Proportion of cluster 2 per condition
prop_table <- gc_all@meta.data %>%
  group_by(condition) %>%
  summarise(
    total = n(),
    n_cluster2 = sum(seurat_clusters == 2),
    pct_cluster2 = round(n_cluster2 / total * 100, 1)
  ) %>%
  arrange(desc(pct_cluster2))
print(prop_table)


ggplot(prop_table, aes(x = reorder(condition, -pct_cluster2), y = pct_cluster2, fill = condition)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "% Cells in Cluster 2 (POSTN+)") +
  ggtitle("Proportion of POSTN+ Mature GCs by Condition")

## Below is the above plot but corrected for nCount confound
# Logistic regression: is a cell in cluster 2, controlling for nCount?
gc_all$is_cluster2 <- as.numeric(gc_all$seurat_clusters == 2)

model <- glm(is_cluster2 ~ nCount_Xenium + condition, 
             data = gc_all@meta.data, family = binomial)
summary(model)

# Extract condition effects as predicted probabilities at mean nCount
library(broom)

# Predict probability of being cluster 2 at the mean nCount for each condition
mean_ncount <- mean(gc_all$nCount_Xenium)

pred_data <- data.frame(
  condition = unique(gc_all$condition),
  nCount_Xenium = mean_ncount
)

pred_data$predicted_pct <- predict(model, newdata = pred_data, type = "response") * 100

pred_data <- pred_data %>% arrange(desc(predicted_pct))
print(pred_data)

# Plot corrected proportions
ggplot(pred_data, aes(x = reorder(condition, -predicted_pct), y = predicted_pct, fill = condition)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "% Cells in Cluster 2 (POSTN+)\n(corrected for nCount)") +
  ggtitle("Proportion of POSTN+ Mature GCs by Condition (nCount-corrected)")

# Side by side comparison: raw vs corrected
prop_table$type <- "Raw"
pred_table <- pred_data %>% 
  rename(pct_cluster2 = predicted_pct) %>%
  mutate(type = "nCount-corrected")

combined_props <- bind_rows(
  prop_table %>% select(condition, pct_cluster2, type),
  pred_table %>% select(condition, pct_cluster2, type)
)

ggplot(combined_props, aes(x = reorder(condition, -pct_cluster2), y = pct_cluster2, fill = type)) +
  geom_col(position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "% Cells in Cluster 2 (POSTN+)", fill = NULL) +
  ggtitle("POSTN+ Proportion: Raw vs nCount-corrected")

# Find genes correlated with POSTN at the single-cell level
expr_mat <- GetAssayData(gc_all, layer = "data", assay = "Xenium")
postn_expr <- expr_mat["POSTN", ]

postn_cors <- sapply(rownames(expr_mat), function(g) {
  cor(postn_expr, expr_mat[g, ], method = "spearman")
})

postn_cors_df <- data.frame(
  gene = names(postn_cors),
  rho = postn_cors
) %>% arrange(desc(abs(rho)))

# Top co-expressed genes
message("Top 15 genes co-expressed with POSTN:")
head(postn_cors_df, 15)

# Top anti-correlated genes
message("Top 15 genes anti-correlated with POSTN:")
tail(postn_cors_df, 15)

# Plot
ggplot(head(postn_cors_df, 20), aes(x = reorder(gene, rho), y = rho)) +
  geom_col(aes(fill = rho > 0)) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = c("TRUE" = "darkblue", "FALSE" = "red"), guide = "none") +
  labs(x = NULL, y = "Spearman rho with POSTN") +
  ggtitle("Top 20 Genes Co-expressed with POSTN")

# is POSTN expression confounded by nCount? -> no
cor.test(postn_expr, gc_all$nCount_Xenium, method = "spearman")

# 1. Statistical test — is the POSTN+ proportion significantly different across conditions?
chisq_data <- gc_all@meta.data %>%
  mutate(is_cluster2 = seurat_clusters == 2)
chisq.test(table(chisq_data$condition, chisq_data$is_cluster2))

# 2. Within postmortem only — is HIBI significantly lower than CARDIAC?
pm_data <- gc_all@meta.data %>%
  filter(condition != "EPILEPSY") %>%
  mutate(is_cluster2 = seurat_clusters == 2)
chisq.test(table(pm_data$condition, pm_data$is_cluster2))

# 3. Within epilepsy — does ILAE score affect POSTN+ proportion?
epi_data <- gc_all@meta.data %>%
  filter(condition == "EPILEPSY") %>%
  mutate(is_cluster2 = seurat_clusters == 2)

epi_prop <- epi_data %>%
  group_by(ILAE_score) %>%
  summarise(pct_cluster2 = round(mean(is_cluster2) * 100, 1))
print(epi_prop)

ggplot(epi_prop, aes(x = ILAE_score, y = pct_cluster2)) +
  geom_col(fill = "steelblue") +
  theme_classic() +
  labs(x = "ILAE Score", y = "% POSTN+ GCs") +
  ggtitle("POSTN+ Proportion by Epilepsy Severity")

# ILAE 3 vs the rest
epi_data$is_ilae3 <- epi_data$ILAE_score == "3"
chisq.test(table(epi_data$is_ilae3, epi_data$is_cluster2))

# ILAE 2 vs ILAE 3 specifically
epi_2v3 <- epi_data %>% filter(ILAE_score %in% c("2", "3"))
chisq.test(table(epi_2v3$ILAE_score, epi_2v3$is_cluster2))

# Trend test across all ILAE scores
prop.trend.test(
  x = c(sum(epi_data$is_cluster2[epi_data$ILAE_score == "0"]),
        sum(epi_data$is_cluster2[epi_data$ILAE_score == "1"]),
        sum(epi_data$is_cluster2[epi_data$ILAE_score == "2"]),
        sum(epi_data$is_cluster2[epi_data$ILAE_score == "3"])),
  n = c(sum(epi_data$ILAE_score == "0"),
        sum(epi_data$ILAE_score == "1"),
        sum(epi_data$ILAE_score == "2"),
        sum(epi_data$ILAE_score == "3"))
)

epi_by_sample <- epi_data %>%
  group_by(orig.ident, ILAE_score) %>%
  summarise(pct_cluster2 = round(mean(is_cluster2) * 100, 1),
            n_cells = n(), .groups = "drop")

ggplot(epi_by_sample, aes(x = ILAE_score, y = pct_cluster2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(size = n_cells), alpha = 0.6) +
  geom_text(aes(label = orig.ident), vjust = -0.8, size = 2.5) +
  theme_classic() +
  labs(x = "ILAE Score", y = "% POSTN+ GCs") +
  ggtitle("POSTN+ Proportion by Sample (Epilepsy)")

# Add ILAE-aware grouping to postn_by_fov
postn_by_fov$group_label <- postn_by_fov$condition
postn_by_fov$group_label[postn_by_fov$orig.ident %in% c("B1","D11","D8","F1","H1","I1","J1","L20")] <- "Epilepsy ILAE 0"
postn_by_fov$group_label[postn_by_fov$orig.ident %in% c("C1","C2","E008","L14","L15","L18")] <- "Epilepsy ILAE 1"
postn_by_fov$group_label[postn_by_fov$orig.ident %in% c("D12","D13","L1","L21","L3")] <- "Epilepsy ILAE 2"
postn_by_fov$group_label[postn_by_fov$orig.ident %in% c("D14","D9","E015","L2")] <- "Epilepsy ILAE 3"

ggplot(postn_by_fov, aes(x = reorder(orig.ident, -pct_expressing), y = pct_expressing, fill = group_label)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "Epilepsy ILAE 0" = "#4DAF4A",
    "Epilepsy ILAE 1" = "#377EB8",
    "Epilepsy ILAE 2" = "#FF7F00",
    "Epilepsy ILAE 3" = "#E41A1C",
    "CARDIAC" = "#984EA3",
    "HIBI" = "#A65628",
    "MAID" = "#F781BF",
    "METASTASIS" = "#999999",
    "SEPSIS" = "#FFD700",
    "TBI" = "#66C2A5"
  )) +
  labs(x = "Sample", y = "% Cells Expressing POSTN", fill = "Group") +
  ggtitle("POSTN Detection Rate by FOV")

DefaultAssay(gc_all) <- "Xenium"
postn_expr <- GetAssayData(gc_all, layer = "data", assay = "Xenium")["POSTN", ]

# What % of cluster 2 expresses POSTN?
cluster2_cells <- gc_all$seurat_clusters == 2
message("Cluster 2 cells expressing POSTN: ", 
        sum(postn_expr[cluster2_cells] > 0), " / ", sum(cluster2_cells),
        " (", round(mean(postn_expr[cluster2_cells] > 0) * 100, 1), "%)")

# Also check POSTN expression in other clusters for comparison
tapply(postn_expr > 0, gc_all$seurat_clusters, function(x) round(mean(x) * 100, 1))

# Create POSTN status column
gc_all$postn_status <- ifelse(postn_expr > 0, "POSTN+", "POSTN-")

# Stacked bar plots
cluster2_data <- gc_all@meta.data %>% filter(seurat_clusters == 2)

# Calculate percentages for POSTN+
postn_pos_pct <- cluster2_data %>%
  filter(postn_status == "POSTN+") %>%
  count(condition) %>%
  mutate(pct = n / sum(n),
         label = paste0(round(pct * 100, 1), "%"))

# Calculate percentages for POSTN-
postn_neg_pct <- cluster2_data %>%
  filter(postn_status == "POSTN-") %>%
  count(condition) %>%
  mutate(pct = n / sum(n),
         label = paste0(round(pct * 100, 1), "%"))

p1 <- ggplot(postn_pos_pct, aes(x = "POSTN+", y = pct, fill = condition)) +
  geom_col(position = "stack") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  theme_classic() +
  labs(x = NULL, y = "Proportion", fill = "Condition") +
  ggtitle(paste0("POSTN+ cells (n=", sum(postn_pos_pct$n), ")"))

p2 <- ggplot(postn_neg_pct, aes(x = "POSTN-", y = pct, fill = condition)) +
  geom_col(position = "stack") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  theme_classic() +
  labs(x = NULL, y = "Proportion", fill = "Condition") +
  ggtitle(paste0("POSTN- cells (n=", sum(postn_neg_pct$n), ")"))

p1 + p2

# Add ILAE-aware grouping
gc_all$group_label <- gc_all$condition
gc_all$group_label[gc_all$orig.ident %in% c("B1","D11","D8","F1","H1","I1","J1","L20")] <- "Epilepsy ILAE 0"
gc_all$group_label[gc_all$orig.ident %in% c("C1","C2","E008","L14","L15","L18")] <- "Epilepsy ILAE 1"
gc_all$group_label[gc_all$orig.ident %in% c("D12","D13","L1","L21","L3")] <- "Epilepsy ILAE 2"
gc_all$group_label[gc_all$orig.ident %in% c("D14","D9","E015","L2")] <- "Epilepsy ILAE 3"

# Calculate POSTN+ percentage per FOV
postn_by_fov_status <- gc_all@meta.data %>%
  group_by(orig.ident, group_label, postn_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(pct = n / sum(n),
         label = paste0(round(pct * 100, 1), "%"))

# Order FOVs by group label then by POSTN+ percentage
fov_order <- postn_by_fov_status %>%
  filter(postn_status == "POSTN+") %>%
  arrange(group_label, desc(pct)) %>%
  pull(orig.ident)

postn_by_fov_status$orig.ident <- factor(postn_by_fov_status$orig.ident, levels = fov_order)

ggplot(postn_by_fov_status, aes(x = orig.ident, y = pct, fill = postn_status)) +
  geom_col() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2) +
  facet_grid(~ group_label, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text = element_text(size = 7)) +
  scale_fill_manual(values = c("POSTN+" = "#E41A1C", "POSTN-" = "#377EB8")) +
  labs(x = NULL, y = "Proportion", fill = NULL) +
  ggtitle("POSTN Expression by FOV (grouped by condition/ILAE)")

DefaultAssay(gc_all) <- "Xenium"

###############################################################################
# 1. CO-EXPRESSION ANALYSIS                                                   #
###############################################################################

# First check which of these markers are in the Xenium panel
markers_to_check <- c(
  # Neuronal
  "PROX1", "RBFOX3", "SNAP25", "SLC17A7",
  # Astrocyte
  "GFAP", "AQP4", "GJA1", "SOX9", "ALDH1L1",
  # Fibroblast/ECM
  "COL1A1", "PDGFRA", "DCN",
  # Endothelial
  "PECAM1",
  # Oligodendrocyte
  "OPALIN", "MBP",
  # The gene of interest
  "POSTN"
)

available <- markers_to_check[markers_to_check %in% rownames(gc_all[["Xenium"]])]
missing <- markers_to_check[!markers_to_check %in% rownames(gc_all[["Xenium"]])]
message("Available: ", paste(available, collapse = ", "))
message("Not in panel: ", paste(missing, collapse = ", "))

# Get expression for available markers in cluster 2 only
cluster2 <- subset(gc_all, subset = seurat_clusters == 2)
expr_mat <- GetAssayData(cluster2, layer = "data", assay = "Xenium")

# Co-expression matrix — what percentage of POSTN+ cells also express each marker?
postn_positive <- expr_mat["POSTN", ] > 0
postn_negative <- expr_mat["POSTN", ] == 0

coexpr_table <- sapply(available[available != "POSTN"], function(gene) {
  c(
    postn_pos_pct = round(mean(expr_mat[gene, postn_positive] > 0) * 100, 1),
    postn_neg_pct = round(mean(expr_mat[gene, postn_negative] > 0) * 100, 1)
  )
})

coexpr_df <- as.data.frame(t(coexpr_table))
coexpr_df$gene <- rownames(coexpr_df)
coexpr_df$diff <- coexpr_df$postn_pos_pct - coexpr_df$postn_neg_pct
coexpr_df <- coexpr_df %>% arrange(desc(abs(diff)))

message("Co-expression with POSTN in Cluster 2:")
print(coexpr_df)

# Visualize as dot plot — neuronal vs non-neuronal markers
# Split cluster 2 by POSTN status
cluster2$postn_status <- ifelse(expr_mat["POSTN", ] > 0, "POSTN+", "POSTN-")
Idents(cluster2) <- "postn_status"

DotPlot(cluster2, features = available[available != "POSTN"]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marker Expression: POSTN+ vs POSTN- in Cluster 2")

# Correlation heatmap of markers within cluster 2
cor_mat <- cor(t(as.matrix(expr_mat[available, ])), method = "spearman")

heatmap(cor_mat, main = "Marker Correlation in Cluster 2",
        scale = "none", col = colorRampPalette(c("blue", "white", "red"))(50))

###############################################################################
# 2. QC METRIC INSPECTION                                                     #
###############################################################################

# Compare nCount and nFeature across clusters
p1 <- VlnPlot(gc_all, features = "nCount_Xenium", group.by = "seurat_clusters", 
              pt.size = 0) + 
  ggtitle("nCount by Cluster") +
  theme(panel.grid.major.y = element_line(color = "grey80"))

p2 <- VlnPlot(gc_all, features = "nFeature_Xenium", group.by = "seurat_clusters", 
              pt.size = 0) + 
  ggtitle("nFeature by Cluster") +
  theme(panel.grid.major.y = element_line(color = "grey80"))

p1 + p2

# Quantify: are cluster 2 cells doublet-like?
gc_all@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    mean_nCount = round(mean(nCount_Xenium), 1),
    median_nCount = round(median(nCount_Xenium), 1),
    mean_nFeature = round(mean(nFeature_Xenium), 1),
    median_nFeature = round(median(nFeature_Xenium), 1)
  )

# Within cluster 2: compare POSTN+ vs POSTN- QC metrics
cluster2_meta <- cluster2@meta.data
cluster2_meta$postn_status <- ifelse(
  expr_mat["POSTN", match(rownames(cluster2_meta), colnames(expr_mat))] > 0, 
  "POSTN+", "POSTN-"
)

cluster2_meta %>%
  group_by(postn_status) %>%
  summarise(
    n_cells = n(),
    mean_nCount = round(mean(nCount_Xenium), 1),
    mean_nFeature = round(mean(nFeature_Xenium), 1)
  )

# Statistical test: do POSTN+ cells have higher counts (doublet signature)?
wilcox.test(nCount_Xenium ~ postn_status, data = cluster2_meta)
wilcox.test(nFeature_Xenium ~ postn_status, data = cluster2_meta)

# Get expression for cluster 2
cluster2 <- subset(gc_all, subset = seurat_clusters == 2)
expr_mat <- GetAssayData(cluster2, layer = "data", assay = "Xenium")
postn_positive <- expr_mat["POSTN", ] > 0

# Known GC markers — check which are in the panel
gc_markers_check <- c("PROX1", "SNAP25", "SLC17A7", "CALB1", "GRIA1", "GRIN1", 
                      "SYNPR", "SCN1A", "SPHKAP", "COL25A1", "SLC12A5", 
                      "GRIA2", "GRIN2A", "GABRB3", "GABRG2")
gc_markers_avail <- gc_markers_check[gc_markers_check %in% rownames(expr_mat)]

# Calculate co-expression percentage in POSTN+ cells
coexpr_pct <- sapply(gc_markers_avail, function(g) {
  round(mean(expr_mat[g, postn_positive] > 0) * 100, 1)
})

coexpr_plot_df <- data.frame(
  gene = names(coexpr_pct),
  pct = coexpr_pct
) %>% arrange(desc(pct))

coexpr_plot_df$gene <- factor(coexpr_plot_df$gene, levels = coexpr_plot_df$gene)

ggplot(coexpr_plot_df, aes(x = gene, y = pct, fill = pct)) +
  geom_col() +
  geom_text(aes(label = paste0(pct, "%")), vjust = -0.5, size = 3) +
  scale_fill_gradient(low = "#377EB8", high = "#E41A1C", guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "% of POSTN+ Cells Co-expressing") +
  ylim(0, 105) +
  ggtitle("GC Marker Co-expression in POSTN+ Cells (Cluster 2)")

###############
# Known immature GC markers — check which are in the panel
immature_markers_check <- c(
  # Classic immature GC markers
  "DCX",       # Doublecortin — gold standard immature neuron marker
  "CALB2",     # Calretinin — marks immature GCs (vs CALB1 for mature)
  "GRIN2B",    # NR2B NMDA subunit — dominant in immature neurons
  "GRIN3A",    # NR3A — present in immature, downregulated with maturation
  "NEUROD1",   # Neuronal differentiation TF
  "SOX11",     # Immature neuron TF
  "STMN1",     # Stathmin — cytoskeletal dynamics in developing neurons
  "GAP43",     # Growth cone marker — axon outgrowth
  "DNER",      # Notch ligand, neurogenesis
  "FOXG1",     # Forebrain development TF
  "NOTCH1",    # Developmental signaling
  "MKI67",     # Proliferation
  "TOP2A",     # Proliferation
  "VCAN",      # ECM remodeling in developing tissue
  "APP",       # Synapse formation/plasticity
  "NNAT",      # Neuronatin — expressed in immature neurons
  "PCSK1",     # Prohormone convertase
  "HOMER1"     # Activity/immediate early gene
)

immature_avail <- immature_markers_check[immature_markers_check %in% rownames(expr_mat)]
immature_missing <- immature_markers_check[!immature_markers_check %in% rownames(expr_mat)]
message("Available: ", paste(immature_avail, collapse = ", "))
message("Not in panel: ", paste(immature_missing, collapse = ", "))

# Calculate co-expression percentage in POSTN+ cells
immature_coexpr <- sapply(immature_avail, function(g) {
  round(mean(expr_mat[g, postn_positive] > 0) * 100, 1)
})

immature_df <- data.frame(
  gene = names(immature_coexpr),
  pct = immature_coexpr
) %>% arrange(desc(pct))

immature_df$gene <- factor(immature_df$gene, levels = immature_df$gene)

ggplot(immature_df, aes(x = gene, y = pct, fill = pct)) +
  geom_col() +
  geom_text(aes(label = paste0(pct, "%")), vjust = -0.5, size = 3) +
  scale_fill_gradient(low = "#377EB8", high = "#E41A1C", guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "% of POSTN+ Cells Co-expressing") +
  ylim(0, 105) +
  ggtitle("Immature GC Marker Co-expression in POSTN+ Cells (Cluster 2)")

# Side by side: do POSTN+ cells have MORE or LESS immature markers?
immature_comparison <- sapply(immature_avail, function(g) {
  c(
    postn_pos = round(mean(expr_mat[g, postn_positive] > 0) * 100, 1),
    postn_neg = round(mean(expr_mat[g, !postn_positive] > 0) * 100, 1)
  )
})

comp_df <- as.data.frame(t(immature_comparison))
comp_df$gene <- rownames(comp_df)
comp_df$diff <- comp_df$postn_pos - comp_df$postn_neg

comp_long <- comp_df %>%
  pivot_longer(cols = c(postn_pos, postn_neg), names_to = "group", values_to = "pct") %>%
  mutate(group = ifelse(group == "postn_pos", "POSTN+", "POSTN-"))

ggplot(comp_long, aes(x = reorder(gene, -pct), y = pct, fill = group)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = paste0(pct, "%")), position = position_dodge(width = 0.9),
            vjust = -0.5, size = 2.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("POSTN+" = "#E41A1C", "POSTN-" = "#377EB8")) +
  labs(x = NULL, y = "% Cells Expressing", fill = NULL) +
  ylim(0, 110) +
  ggtitle("Immature Markers: POSTN+ vs POSTN- in Cluster 2")

########
DefaultAssay(gc_all) <- "Xenium"
expr_mat <- GetAssayData(gc_all, layer = "data", assay = "Xenium")

# Calculate per-FOV stats
fov_stats <- gc_all@meta.data %>%
  mutate(
    postn_pos = expr_mat["POSTN", colnames(gc_all)] > 0,
    prox1_pos = expr_mat["PROX1", colnames(gc_all)] > 0
  ) %>%
  group_by(orig.ident, condition) %>%
  summarise(
    total_gc = n(),
    n_postn = sum(postn_pos),
    n_prox1 = sum(prox1_pos),
    pct_postn = round(mean(postn_pos) * 100, 1),
    pct_prox1 = round(mean(prox1_pos) * 100, 1),
    mean_nCount = round(mean(nCount_Xenium), 1),
    .groups = "drop"
  )

print(fov_stats)

# 1. Total GCs vs number of POSTN+ cells
cor1 <- cor.test(fov_stats$total_gc, fov_stats$n_postn, method = "spearman")

p1 <- ggplot(fov_stats, aes(x = total_gc, y = n_postn, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text(aes(label = orig.ident), vjust = -0.8, size = 2.5) +
  theme_classic() +
  labs(x = "Total GC Cells per FOV", y = "Number of POSTN+ Cells") +
  annotate("text", x = max(fov_stats$total_gc) * 0.3, y = max(fov_stats$n_postn) * 0.9,
           label = paste0("rho = ", round(cor1$estimate, 3), "\np = ", signif(cor1$p.value, 3))) +
  ggtitle("GC Count vs POSTN+ Count per FOV")

# 2. Total GCs vs number of PROX1+ cells
cor2 <- cor.test(fov_stats$total_gc, fov_stats$n_prox1, method = "spearman")

p2 <- ggplot(fov_stats, aes(x = total_gc, y = n_prox1, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text(aes(label = orig.ident), vjust = -0.8, size = 2.5) +
  theme_classic() +
  labs(x = "Total GC Cells per FOV", y = "Number of PROX1+ Cells") +
  annotate("text", x = max(fov_stats$total_gc) * 0.3, y = max(fov_stats$n_prox1) * 0.9,
           label = paste0("rho = ", round(cor2$estimate, 3), "\np = ", signif(cor2$p.value, 3))) +
  ggtitle("GC Count vs PROX1+ Count per FOV")

p1 + p2

# 3. POSTN+ vs PROX1+ directly
cor3 <- cor.test(fov_stats$n_postn, fov_stats$n_prox1, method = "spearman")

p3 <- ggplot(fov_stats, aes(x = n_prox1, y = n_postn, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text(aes(label = orig.ident), vjust = -0.8, size = 2.5) +
  theme_classic() +
  labs(x = "PROX1+ Cells per FOV", y = "POSTN+ Cells per FOV") +
  annotate("text", x = max(fov_stats$n_prox1) * 0.3, y = max(fov_stats$n_postn) * 0.9,
           label = paste0("rho = ", round(cor3$estimate, 3), "\np = ", signif(cor3$p.value, 3))) +
  ggtitle("PROX1+ vs POSTN+ Cells per FOV")

p3

# 4. Use PERCENTAGES instead of raw counts to control for FOV size
cor4 <- cor.test(fov_stats$pct_postn, fov_stats$pct_prox1, method = "spearman")

p4 <- ggplot(fov_stats, aes(x = pct_prox1, y = pct_postn, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text(aes(label = orig.ident), vjust = -0.8, size = 2.5) +
  theme_classic() +
  labs(x = "% Cells Expressing PROX1", y = "% Cells Expressing POSTN") +
  annotate("text", x = max(fov_stats$pct_prox1) * 0.3, y = max(fov_stats$pct_postn) * 0.9,
           label = paste0("rho = ", round(cor4$estimate, 3), "\np = ", signif(cor4$p.value, 3))) +
  ggtitle("PROX1 vs POSTN Detection Rate per FOV")

p4

# 5. Both detection rates vs nCount to see if it's all technical
cor5 <- cor.test(fov_stats$mean_nCount, fov_stats$pct_postn, method = "spearman")
cor6 <- cor.test(fov_stats$mean_nCount, fov_stats$pct_prox1, method = "spearman")

message("POSTN% vs nCount: rho = ", round(cor5$estimate, 3))
message("PROX1% vs nCount: rho = ", round(cor6$estimate, 3))

# 1. nCount vs POSTN+ count
cor1 <- cor.test(fov_stats$mean_nCount, fov_stats$n_postn, method = "spearman")

p7 <- ggplot(fov_stats, aes(x = mean_nCount, y = n_postn, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text(aes(label = orig.ident), vjust = -0.8, size = 2.5) +
  theme_classic() +
  labs(x = "Mean nCount per FOV", y = "Number of POSTN+ Cells") +
  annotate("text", x = max(fov_stats$mean_nCount) * 0.3, y = max(fov_stats$n_postn) * 0.9,
           label = paste0("rho = ", round(cor1$estimate, 3), "\np = ", signif(cor1$p.value, 3))) +
  ggtitle("nCount vs POSTN+ Count per FOV")

# 2. nCount vs PROX1+ count
cor2 <- cor.test(fov_stats$mean_nCount, fov_stats$n_prox1, method = "spearman")

p8 <- ggplot(fov_stats, aes(x = mean_nCount, y = n_prox1, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text(aes(label = orig.ident), vjust = -0.8, size = 2.5) +
  theme_classic() +
  labs(x = "Mean nCount per FOV", y = "Number of PROX1+ Cells") +
  annotate("text", x = max(fov_stats$mean_nCount) * 0.3, y = max(fov_stats$n_prox1) * 0.9,
           label = paste0("rho = ", round(cor2$estimate, 3), "\np = ", signif(cor2$p.value, 3))) +
  ggtitle("nCount vs PROX1+ Count per FOV")

p7 + p8

### #########################################
# NOTE: Oct 20 start here

library(slingshot)
library(SingleCellExperiment)
library(tidyr)
library(dplyr)

DefaultAssay(gc_all) <- "Xenium"

###############################################################################
# DEFINE TRAJECTORY STAGES FROM LITERATURE                                     #
#                                                                             #
# Established DG granule cell maturation (Hochgerner 2018, Mukherjee 2024):   #
#   Neuroblast → Immature GC → Maturing GC → Mature GC                       #
#                                                                             #
# Stage 1 (Neuroblast/earliest GC): SOX2, NES, PAX6, MKI67, TOP2A, SOX4,    #
#          SOX11, CALB2                                                       #
# Stage 2 (Immature GC): CALB2, GRIN2B, GRIN3A, NNAT, DNER                  #
# Stage 3 (Mature GC): PROX1, CALB1, GRIN2A, SLC12A5                        #
# Stage 4 (ECM-anchored GC): POSTN, COL25A1, VCAN                           #
###############################################################################

# Score cells for each stage
early_genes <- list(c("SOX2", "NES", "PAX6", "MKI67", "TOP2A", "SOX4", "SOX11"))
immature_genes <- list(c("CALB2", "GRIN2B", "GRIN3A", "NNAT", "DNER"))
mature_genes <- list(c("PROX1", "CALB1", "GRIN2A", "SLC12A5"))
ecm_genes <- list(c("POSTN", "COL25A1", "VCAN"))

gc_all <- AddModuleScore(gc_all, features = early_genes, name = "early_score",
                         assay = "Xenium", ctrl = 5, nbin = 10)
gc_all <- AddModuleScore(gc_all, features = immature_genes, name = "immature_score",
                         assay = "Xenium", ctrl = 5, nbin = 10)
gc_all <- AddModuleScore(gc_all, features = mature_genes, name = "mature_score",
                         assay = "Xenium", ctrl = 5, nbin = 10)
gc_all <- AddModuleScore(gc_all, features = ecm_genes, name = "ecm_score",
                         assay = "Xenium", ctrl = 5, nbin = 10)

# Verify scores look sensible on UMAP
p1 <- FeaturePlot(gc_all, features = "early_score1", raster = TRUE, cols = c("#F0F0F0", "red")) + ggtitle("Early/Neuroblast Score")
p2 <- FeaturePlot(gc_all, features = "immature_score1", raster = TRUE, cols = c("#F0F0F0", "orange")) + ggtitle("Immature GC Score")
p3 <- FeaturePlot(gc_all, features = "mature_score1", raster = TRUE, cols = c("#F0F0F0", "blue")) + ggtitle("Mature GC Score")
p4 <- FeaturePlot(gc_all, features = "ecm_score1", raster = TRUE, cols = c("#F0F0F0", "darkgreen")) + ggtitle("ECM Score")
p1 + p2 + p3 + p4


###############################################################################
# RUN SLINGSHOT PSEUDOTIME                                                     #
###############################################################################

# Identify start and end clusters based on literature-defined scores
cluster_early <- tapply(gc_all$early_score1, gc_all$seurat_clusters, mean)
cluster_mature <- tapply(gc_all$mature_score1, gc_all$seurat_clusters, mean)

start_cluster <- names(which.max(cluster_early))
end_cluster <- names(which.max(cluster_mature))

message("Start cluster (highest early/neuroblast score): ", start_cluster)
message("End cluster (highest mature GC score): ", end_cluster)

# Convert to SingleCellExperiment and run slingshot
sce <- as.SingleCellExperiment(gc_all)
# Copy PCA embeddings from Seurat to SCE
reducedDim(sce, "PCA") <- Embeddings(gc_all, "pca")

# Now run slingshot
sce <- slingshot(sce, clusterLabels = "seurat_clusters", reducedDim = "PCA",
                 start.clus = start_cluster, end.clus = end_cluster)

gc_all$pseudotime <- slingPseudotime(sce)[, 1]

# Verify
summary(gc_all$pseudotime)

sce <- slingshot(sce, clusterLabels = "seurat_clusters", reducedDim = "PCA",
                 start.clus = start_cluster, end.clus = end_cluster)

gc_all$pseudotime <- slingPseudotime(sce)[, 1]

# Plot pseudotime on UMAP
# Get the pseudotime range to fix the scale
pt_range <- range(gc_all$pseudotime, na.rm = TRUE)
message("Pseudotime range: ", round(pt_range[1], 2), " to ", round(pt_range[2], 2))

FeaturePlot(gc_all, features = "pseudotime", raster = TRUE,
            cols = c("yellow", "red", "darkblue"),
            min.cutoff = pt_range[1], max.cutoff = pt_range[2]) +
  coord_fixed() + 
  scale_color_gradientn(colors = c("yellow", "red", "darkblue"),
                        limits = pt_range) +
  ggtitle("Pseudotime (yellow = immature, blue = mature)")

###############################################################################
# GENE EXPRESSION ALONG PSEUDOTIME                                            #
###############################################################################

expr_mat <- GetAssayData(gc_all, layer = "data", assay = "Xenium")

all_trajectory_genes <- c("SOX2", "NES", "PAX6", "MKI67", "TOP2A", "SOX4", "SOX11",
                          "CALB2", "GRIN2B", "GRIN3A", "NNAT", "DNER",
                          "PROX1", "CALB1", "GRIN2A", "SLC12A5",
                          "POSTN", "COL25A1", "VCAN")

# Build data frame
pt_data <- data.frame(
  pseudotime = gc_all$pseudotime,
  cell = colnames(gc_all)
)

for (g in all_trajectory_genes) {
  pt_data[[g]] <- expr_mat[g, ]
}

pt_data <- pt_data %>% filter(!is.na(pseudotime))

# Assign literature-based stage categories
gene_stages <- data.frame(gene = all_trajectory_genes) %>%
  mutate(stage = case_when(
    gene %in% c("SOX2", "NES", "PAX6", "MKI67", "TOP2A", "SOX4", "SOX11") ~ "1. Early/Neuroblast",
    gene %in% c("CALB2", "GRIN2B", "GRIN3A", "NNAT", "DNER") ~ "2. Immature GC",
    gene %in% c("PROX1", "CALB1", "GRIN2A", "SLC12A5") ~ "3. Mature GC",
    gene %in% c("POSTN", "COL25A1", "VCAN") ~ "4. ECM/Structural"
  ))

pt_long <- pt_data %>%
  pivot_longer(cols = all_of(all_trajectory_genes), names_to = "gene", values_to = "expression") %>%
  left_join(gene_stages, by = "gene")

# Smoothed curves faceted by stage
ggplot(pt_long, aes(x = pseudotime, y = expression, color = gene)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3, linewidth = 1) +
  facet_wrap(~ stage, ncol = 1, scales = "free_y") +
  theme_classic() +
  labs(x = "Pseudotime (immature → mature)", y = "Expression") +
  ggtitle("Gene Expression Along GC Maturation Trajectory\n(Literature-defined markers)")


###############################################################################
# POSTN TIMING — KEY PLOT                                                      #
###############################################################################

# One representative gene per stage + POSTN
key_comparison <- c("SOX2", "CALB2", "GRIN3A", "CALB1", "PROX1", "POSTN", "COL25A1")

pt_key <- pt_data %>%
  pivot_longer(cols = all_of(key_comparison), names_to = "gene", values_to = "expression") %>%
  left_join(gene_stages, by = "gene")

ggplot(pt_key, aes(x = pseudotime, y = expression, color = gene)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3, linewidth = 1.2) +
  theme_classic() +
  scale_color_manual(values = c(
    "SOX2" = "#984EA3",    # purple — earliest
    "CALB2" = "#FF7F00",   # orange — neuroblast/immature
    "GRIN3A" = "#A65628",  # brown — immature
    "CALB1" = "#377EB8",   # blue — mature
    "PROX1" = "#4DAF4A",   # green — mature
    "POSTN" = "#E41A1C",   # red — target gene
    "COL25A1" = "#F781BF"  # pink — ECM
  )) +
  labs(x = "Pseudotime (immature → mature)", y = "Expression", color = "Gene") +
  ggtitle("POSTN Timing Relative to Established Maturation Markers")


###############################################################################
# BINNED VERSION — cleaner visualization                                       #
###############################################################################

pt_data$pt_bin <- cut(pt_data$pseudotime, breaks = 25, labels = FALSE)

binned <- pt_data %>%
  group_by(pt_bin) %>%
  summarise(across(all_of(key_comparison), mean), .groups = "drop")

binned_long <- binned %>%
  pivot_longer(cols = all_of(key_comparison), names_to = "gene", values_to = "expression")

ggplot(binned_long, aes(x = pt_bin, y = expression, color = gene)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  theme_classic() +
  scale_color_manual(values = c(
    "SOX2" = "#984EA3", "CALB2" = "#FF7F00", "GRIN3A" = "#A65628",
    "CALB1" = "#377EB8", "PROX1" = "#4DAF4A",
    "POSTN" = "#E41A1C", "COL25A1" = "#F781BF"
  )) +
  labs(x = "Pseudotime Bin (1 = immature → 25 = mature)", y = "Mean Expression") +
  ggtitle("Binned Expression Along Maturation Trajectory")


###############################################################################
# COMPARE PSEUDOTIME ACROSS CONDITIONS                                         #
###############################################################################

# Do conditions differ in where their cells sit along the trajectory?
ggplot(gc_all@meta.data %>% filter(!is.na(pseudotime)), 
       aes(x = condition, y = pseudotime, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Pseudotime Distribution by Condition")

# Density plot
ggplot(gc_all@meta.data %>% filter(!is.na(pseudotime)), 
       aes(x = pseudotime, fill = condition)) +
  geom_density(alpha = 0.4) +
  theme_classic() +
  ggtitle("Pseudotime Density by Condition")
