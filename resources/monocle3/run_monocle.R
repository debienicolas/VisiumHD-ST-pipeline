#!/usr/bin/env Rscript

library(Seurat)
library(monocle3)
library(ggplot2)
# library(tidyverse)
library(optparse)
library(Matrix)
library(igraph)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 64000 * 1024^2)  # Set to 64GB

# Takes as input the annotated anndata object

## set the input parameteres
param_list <- list(
    make_option("--input_path", type="character", help="Path to the annotated anndata object."),
    make_option("--data_dir", type="character", help="Path to the data Visium data folder"),
    make_option("--output_path", type="character", help="Path to the output folder"),
    make_option("--config", type="character", help="Path to the config file")
)

params <- parse_args(OptionParser(option_list=param_list))

config <- jsonlite::fromJSON(params$config)
if (is.null(params$input_path)) {
    params$input_path = "results/human_breast_cancer_final/infercnvpy/infercnvpy_output.h5ad"
    params$data_dir = "input/Visium_HD_Human_Breast_Cancer_Fresh_Frozen"
    params$output_path = "results/human_breast_cancer_final/monocle3_spat_right"
    params$config$pt_size_factor = 0.5
    params$config$run_factor = FALSE
    params$config$run_clone = TRUE
}

output_path = params$output_path
data_dir = params$data_dir

# if it doesn't exist, create it
if (!dir.exists(output_path)) {
    dir.create(output_path)
}

#### Load the seurat object from .h5ad file using schard
visium_data <- schard::h5ad2seurat_spatial(params$input_path)

print(visium_data)

######### Create a cell_data_set object ##########

# Extract necessary matrices
expression_matrix <- GetAssayData(visium_data, assay = "Spatial", layer = "counts")
cell_metadata <- visium_data@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

# Create cell_data_set object
cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

print("cell_data object created.")


plot_spatial_pseudotime <- function(cds, visium_data, output_path) {
  # Get pseudotime values
  pseudotime_values <- pseudotime(cds)
  print(paste("Found", length(pseudotime_values), "pseudotime values"))
  # Only include cells that are in both the CDS and Visium object
  common_cells <- intersect(names(pseudotime_values), colnames(visium_data))
  print(paste("Found", length(common_cells), "cells in both the CDS and Visium object"))
  visium_data$pseudotime <- pseudotime_values
  print(paste("Added pseudotime values to Visium object"))

  png(file.path(output_path, "spatial_pseudotime.png"), width=1200, height=1000, res=150)
  print(SpatialFeaturePlot(visium_data, features = "pseudotime", alpha = 0.5, pt.size.factor = config$pt_size_factor, 
                          spatial.legend = TRUE, spatial.legend.position = "right"))
  dev.off()
}

#### Running on individual cell types ####

# Function to run Monocle3 on individual cell types
run_monocle3_by_celltype <- function(cds, output_path, cell_type_column = "factor", auto_root = TRUE) {
  # Get unique cell types from the specified column in metadata
  if (!cell_type_column %in% names(colData(cds))) {
    stop(paste("Column", cell_type_column, "not found in cell metadata. Available columns are:", 
               paste(names(colData(cds)), collapse=", ")))
  }
  
  cell_types <- unique(colData(cds)[[cell_type_column]])
  print(paste("Found", length(cell_types), "cell types in column", cell_type_column))
  
  # Create a directory for each cell type's results
  for (cell_type in cell_types) {
    # Create a safe directory name
    safe_cell_type <- gsub("[^a-zA-Z0-9_-]", "_", as.character(cell_type))
    cell_type_dir <- file.path(output_path, safe_cell_type)
    if (!dir.exists(cell_type_dir)) {
      dir.create(cell_type_dir, recursive = TRUE)
    }
    
    # Subset the cds object for this cell type
    cds_subset <- cds[, colData(cds)[[cell_type_column]] == cell_type]
    
    # Only proceed if we have enough cells
    if (ncol(cds_subset) < 20) {
      print(paste("Skipping", cell_type, "- too few cells:", ncol(cds_subset)))
      next
    }
    
    print(paste("Processing", cell_type, "with", ncol(cds_subset), "cells"))
    
    # Run the Monocle3 workflow on this subset
    tryCatch({
      # Preprocess with more PCs for better dimensionality reduction
      cds_subset <- preprocess_cds(cds_subset, method = "PCA", num_dim = 50)
      print(paste(cell_type, "- PCA done."))
      
      # Save a diagnostic plot to check PCA
      png(file.path(cell_type_dir, "pca_variance.png"), width=1200, height=1000, res=150)
      print(plot_pc_variance_explained(cds_subset))
      dev.off()
      
      # Reduce dimension with UMAP
      cds_subset <- reduce_dimension(cds_subset, reduction_method = "UMAP")
      print(paste(cell_type, "- reduce_dimension done."))
      
      # Check if UMAP coordinates exist
      if ("UMAP" %in% reducedDimNames(cds_subset)) {
        # Plot UMAP colored by the cell type column
        png(file.path(cell_type_dir, paste0("umap_by_", cell_type_column, ".png")), width=1200, height=1000, res=150)
        print(plot_cells(cds_subset, reduction_method = "UMAP", color_cells_by = cell_type_column))
        dev.off()
        
        # Cluster cells
        cds_subset <- cluster_cells(cds_subset, resolution = 1e-3)
        print(paste(cell_type, "- cluster_cells done."))
        
        # Plot clusters
        png(file.path(cell_type_dir, "clusters.png"), width=1200, height=1000, res=150)
        print(plot_cells(cds_subset, reduction_method = "UMAP", color_cells_by = "cluster"))
        dev.off()
        
        # Learn graph
        cds_subset <- learn_graph(cds_subset, use_partition = FALSE)
        print(paste(cell_type, "- learn_graph done."))
        
        # Plot graph
        png(file.path(cell_type_dir, "trajectory.png"), width=1200, height=1000, res=150)
        print(plot_cells(cds_subset, reduction_method = "UMAP", 
                   color_cells_by = cell_type_column, 
                   label_groups_by_cluster = FALSE,
                   label_leaves = FALSE, 
                   label_branch_points = FALSE,
                   graph_label_size = 1.5))
        dev.off()

        png(file.path(cell_type_dir, "trajectory_with_branch_points.png"), width=1200, height=1000, res=150)
        print(plot_cells(cds_subset, reduction_method = "UMAP", 
                   color_cells_by = cell_type_column, 
                   label_groups_by_cluster = FALSE,
                   label_leaves = TRUE, 
                   label_branch_points = TRUE,
                   graph_label_size = 1.5))
        dev.off()
        # Order cells
        if (auto_root) {
          # Find a suitable root node
          get_central_principal_node <- function(cds_obj) {
            g <- principal_graph(cds_obj)[["UMAP"]]
            if (igraph::vcount(g) > 0) {
              centrality_scores <- igraph::closeness(g)
              root_node <- names(centrality_scores)[which.max(centrality_scores)]
              return(root_node)
            } else {
              stop("Graph has no vertices")
            }
          }
          
          tryCatch({
            root_node <- get_central_principal_node(cds_subset)
            print(paste(cell_type, "- selected root node:", root_node))
            cds_subset <- order_cells(cds_subset, root_pr_nodes = root_node)
          }, error = function(e) {
            print(paste("Error in root selection:", e$message))
            print("Using default root selection")
            cds_subset <- order_cells(cds_subset)
          })
        } else {
          cds_subset <- order_cells(cds_subset)
        }
        
        print(paste(cell_type, "- order_cells done."))
        
        # Plot pseudotime
        png(file.path(cell_type_dir, "pseudotime.png"), width=1200, height=1000, res=150)
        print(plot_cells(cds_subset, reduction_method = "UMAP",
                   color_cells_by = "pseudotime",
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   graph_label_size = 1.5))
        dev.off()

        # try to plot the spatial pseudotime
        if (cell_type_column == "factor") {
          visium_subset <- subset(visium_data, subset = factor == cell_type)
          plot_spatial_pseudotime(cds_subset, visium_subset, cell_type_dir)
        } else if (cell_type_column == "cnv_leiden") {
          visium_subset <- subset(visium_data, subset = cnv_leiden == cell_type)
          plot_spatial_pseudotime(cds_subset, visium_subset, cell_type_dir)
        }




      } else {
        print(paste("WARNING: UMAP reduction not found for", cell_type))
      }
      
      # Save the cell_data_set object
      saveRDS(cds_subset, file.path(cell_type_dir, "monocle3_cds.rds"))
      
    }, error = function(e) {
      print(paste("ERROR processing", cell_type, ":", e$message))
    })
  }
  
  return(cell_types)
}

# Call the function with automatic root selection
if (as.logical(config$run_factor)) {
  cell_types <- run_monocle3_by_celltype(cds, file.path(output_path, "factor"), cell_type_column = "factor", auto_root = TRUE)
}
if (as.logical(config$run_clone)) {
  cell_types <- run_monocle3_by_celltype(cds, file.path(output_path, "cnv_leiden"), cell_type_column = "cnv_leiden", auto_root = TRUE)
}


# reduce batch effects

print("preprocessing the cell_data object...")
cds <- preprocess_cds(cds, method = "PCA")
print("PCA done.")
cds <- reduce_dimension(cds)
print("reduce_dimension done.")

png(file.path(output_path, "trajectory_factors.png"), width=1200, height=1000, res=150)
plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "factor")
dev.off()

cds <- cluster_cells(cds)
print("cluster_cells done.")

png(file.path(output_path, "trajectory_clusters.png"), width=1200, height=1000, res=150)
plot_cells(cds, color_cells_by="partition")
dev.off()


cds <- learn_graph(cds)
print("learn_graph done.")

png(file.path(output_path, "learn_graph.png"), width=1200, height=1000, res=150)
plot_cells(cds, color_cells_by="factor", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)
dev.off()

png(file.path(output_path, "learn_graph_2.png"), width=1200, height=1000, res=150)
plot_cells(cds, color_cells_by="factor", label_groups_by_cluster=FALSE, label_leaves=TRUE, label_branch_points=TRUE)
dev.off()


## Ordering cells / selecting root cells


# ways to automatically select the root nodes
# 1.Gene marker based approach -> identify spots with high expression of known early-stage markers
# 2.Cell clusters based approach -> prior knowledge of which cell clusters/types are early states
# 3. Graph based approach -> identify the root principal points as the ones with highest centrality in the trajectory graph
# 4. spatial extrema -> not really convincing

# here we use the graph based approach



#cds <- order_cells(cds)
#print("order_cells done.")



print(cds)

print("Finished running monocle3")


