#!/usr/bin/env Rscript

# This script is used to preprocess the visium data and create a Seurat object from it
# The input folder should be the folder path to the visium input data
# The output file is a the seurat object inside the cellchat output folder

# Seurat vignette: https://satijalab.org/seurat/articles/spatial_vignette.html

library(dplyr)
library(Seurat)
library(optparse)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8000 * 1024^2)  # Set to 8GB

# nFeature_Spatial: the number of genes detected in each cell
# nCount_Spatial: the total number of molecules detected within a cell

param_list <- list(
    make_option("--input_dir", type="character", help="Input directory"),
    make_option("--output_path", type="character", help="Output file path"),
    make_option("--default_as`say", type="character", default="Spatial.016um", help="Default assay")
)

params <- parse_args(OptionParser(option_list=param_list))

# if the input_dir is not provided, use a manually set path or raise an error
if (is.null(params$input_dir)) {
    input_dir <- "/Users/nicolasdebie/st-pipeline/Visium_input_data"
}

# if the output_path is not provided, use a manually set path or raise an error
if (is.null(params$output_path)) {
    output_path <- "/Users/nicolasdebie/st-pipeline/seurat_object.rds"
}


####### 1. Load the visium data as a Seurat object ###########################################################
img = Read10X_Image(image.dir = paste0(input_dir, "/binned_outputs/Square_016um/spatial"), image.name = "tissue_lowres_image.png", assay = "Spatial.016um")
data <- Load10X_Spatial(data.dir = input_dir, bin.size=16, image = img)
# Set default assay
DefaultAssay(data) <- params$default_assay


####### 2. Seurat pre-processing and clustering ##############################################################
# Calculate quality metrics
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")

# Filter data (adjust thresholds as needed)
data <- subset(data, nFeature_Spatial.016um > 750 & nFeature_Spatial.016um < 1500 & percent.mt < 20)
print(nrow(data))
# find the max in nFeature_Spatial.016um
max_nFeature <- max(data$nFeature_Spatial.016um)


# Normalization and Feature Selection
# Normalize data using LogNormalize instead of SCTransform
data <- NormalizeData(data, normalization.method = "LogNormalize")
# creates a new layer called "data" with the normalized data (data[["default_assay"]])

# Find variable features
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
# add 2000 variable features to the seurat object, can be accessed with VariableFeatures(data)

# Scale the data
data <- ScaleData(data)
# scales the normalized data layer

# Run PCA and UMAP
data <- RunPCA(data) # data$pca
data <- RunUMAP(data, dims = 1:30) # data$umap

# Find neighbors
data <- FindNeighbors(data, dims = 1:30)

# Find clusters and ensure labels start from 1 instead of 0
data <- FindClusters(data, resolution = 0.2)


#Idents(data) <- factor(as.numeric(Idents(data)) + 1)  # Add 1 to all cluster numbers

# give a summary of each cluster based on the gene markers

# give a summary of each cluster based on the gene markers
# Find marker genes for each cluster
markers <- FindAllMarkers(data, 
                         only.pos = TRUE,    
                         min.pct = 0.25,     
                         logfc.threshold = 0.25)

# Get cluster statistics
cluster_stats <- data[[]] %>%
  group_by(seurat_clusters) %>%
  summarise(
    num_cells = n(),
    avg_features = mean(nFeature_Spatial.016um),
    avg_counts = mean(nCount_Spatial.016um),
    avg_mt_percent = mean(percent.mt)
  )

# Get top markers per cluster
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  ungroup()

# Print detailed summary for each cluster
print("=== Detailed Cluster Analysis ===")
for(clust in unique(top_markers$cluster)) {
  cat(sprintf("\n\n=== CLUSTER %s ===\n", clust))
  
  # Print cluster statistics
  stats <- cluster_stats[cluster_stats$seurat_clusters == clust,]
  cat(sprintf("Number of cells: %d\n", stats$num_cells))
  cat(sprintf("Average features per cell: %.1f\n", stats$avg_features))
  cat(sprintf("Average counts per cell: %.1f\n", stats$avg_counts))
  cat(sprintf("Average mitochondrial %%: %.1f%%\n", stats$avg_mt_percent))
  
  # Print top markers with statistics
  cat("\nTop marker genes:\n")
  cluster_markers <- top_markers %>% 
    filter(cluster == clust) %>%
    select(gene, avg_log2FC, pct.1, pct.2, p_val_adj)
  
  for(i in 1:nrow(cluster_markers)) {
    cat(sprintf("%s:\n", cluster_markers$gene[i]))
    cat(sprintf("  - log2 fold change: %.2f\n", cluster_markers$avg_log2FC[i]))
    cat(sprintf("  - expressed in %.1f%% of cluster (vs %.1f%% in other cells)\n", 
                cluster_markers$pct.1[i]*100, 
                cluster_markers$pct.2[i]*100))
    cat(sprintf("  - adjusted p-value: %.2e\n", cluster_markers$p_val_adj[i]))
  }
}

# Get cluster names from user input
cluster_names <- readline(prompt="Enter the cluster names (comma-separated): ")
cluster_names <- strsplit(cluster_names, ", ")[[1]]  # Extract the vector from the list

# rename the clusters 
new_cluster_ids <- cluster_names
names(new_cluster_ids) <- levels(Idents(data))
data <- RenameIdents(data, new_cluster_ids)

# Add this line to drop unused levels
data@active.ident <- droplevels(data@active.ident)

####### 3. Save the seurat object ###########################################################################

saveRDS(data, file = output_path)
