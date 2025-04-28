#!/usr/bin/env Rscript

# Vignette: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.html
# Loading the required libraries

library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(extrafont)
library(optparse)
library(Seurat)
library(CellChat)
library(gridExtra)
library(grid)
library(pryr)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 148000 * 1024^2)  # Set to 148GB

gc()  # Force garbage collection before starting

## set the input parameteres
param_list <- list(
    make_option("--input_path", type="character", help="Path to the annotated anndata object."),
    make_option("--output_path", type="character", help="Path to the output folder"),
    make_option("--data_dir", type="character", help="Path to the data directory"),
    make_option("--scale_factor_path", type="character", help="Path to the scale factor file"),
    make_option("--species", type="character", help="Species"),
    make_option("--config", type="character", help="dictionary of parameters")
)

params <- parse_args(OptionParser(option_list=param_list))

## parse the config param
config <- jsonlite::fromJSON(params$config)

if (is.null(params$input_path)) {
    params$input_path = "results/run_4/ficture/annotated_anndata.h5ad"
    params$input_path = "results/8_um_nF_10/ficture/annotated_anndata.h5ad"
    params$input_path = "results/multiple_factors/ficture/annotated_anndata.h5ad"
}

if (is.null(params$output_path)) {
    params$output_path = "results/run_4/cellchat"
    params$output_path = "results/8_um_nF_10/cellchat"
    params$output_path = "results/multiple_factors/cellchat"
}

if (is.null(params$data_dir)) {
    params$data_dir = "input/TG23-0301_VGEX_results_VISIUM_HD"
}


output_path = params$output_path
data_dir = params$data_dir

# if it doesn't exist, create it
if (!dir.exists(output_path)) {
    dir.create(output_path)
}

#seurat_obj_path = "seurat_object.rds"

# load the seurat object
#data <- readRDS(seurat_obj_path)

#### Load the seurat object from .h5ad file using schard

data <- schard::h5ad2seurat_spatial(params$input_path)

# Check matrix type
print("Matrix type:")
print(class(GetAssayData(data, slot = "counts")))
print(class(GetAssayData(data, slot = "data")))

data <- NormalizeData(data, normalization.method = "LogNormalize")

# filter data based on counts
#data <- subset(data, subset = nFeature_Spatial > 10)
# filter out low expression genes


# Downsample the number of cells per identity class
# subset the data 
# print("before downsampling")
# print(data)
# # downsample the data 
# data <- subset(data, downsample = 0000)
# print("after downsampling")
# print(data)


print("Seurat data object: ")
print(data)

# print min and max counts
print(min(data$nFeature_Spatial))
print(max(data$nFeature_Spatial))

# remove the rows with NA in the factor column
#data <- subset(data, subset = !is.na(data$factor))


####### 3. Prepare for CellChat ###########################################################################
# Get normalized data (ensure we maintain sparsity)
data.input <- GetAssayData(data, slot = "data", assay = "Spatial")
print("CellChat input matrix type:")
print(class(data.input))


# Create metadata with adjusted cluster labels
meta <- data.frame(
    labels = data$factor,
    samples = factor("sample1"),  # Make samples a factor explicitly
    row.names = names(Idents(data))
)


# Get spatial coordinates and ensure it's a proper two-column matrix
spatial.locs <- GetTissueCoordinates(data, scale = NULL)
# Check the actual column names
print("Column names of spatial coordinates:")
print(colnames(spatial.locs))
# Use the correct column names (likely 'imagerow' and 'imagecol')
spatial.locs <- as.matrix(spatial.locs[, c("imagerow", "imagecol")])  # Adjust column names
colnames(spatial.locs) <- c("x", "y")  # Rename to expected names

# Set spatial factors for 10X Visium
spot.size = config$spot_size     # theoretical spot size (um)
# Read scalefactors from json file
#scalefactors <- jsonlite::fromJSON(txt = file.path(data_dir, "binned_outputs/square_002um/spatial/scalefactors_json.json"))
scalefactors <- jsonlite::fromJSON(txt = params$scale_factor_path)

# Calculate conversion factor using binned data parameters
# The spot diameter needs to be adjusted based on the bin size
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)




####### 4. Create CellChat object + run CellChat ###########################################################################
cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = spatial.locs,
    spatial.factors = spatial.factors
)
gc()  # Force garbage collection after heavy operations

print("Memory usage after CellChat creation:")
print(pryr::mem_used())

# save a cellchat object to file 
saveRDS(cellchat, file = file.path(output_path, "cellchat_object.rds"))


## Set the ligand-receptor interaction database
if (params$species == "mouse") {
    cellchat@DB <- CellChatDB.mouse
} else if (params$species == "human") {
    cellchat@DB <- CellChatDB.human
} else {
    stop("Species must be either mouse or human")
}
cellchat <- subsetData(cellchat)


##### Set the number of cpus to use equal to the number of workers #####
future::plan("multisession", workers = 3)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)


### Inference of cell-cell communication network

# 250 um is the interaction range for the cellchat analysis
cellchat <- computeCommunProb(cellchat, type = config$type, trim = config$trim, 
            distance.use = TRUE, interaction.range = config$interaction_range, scale.distance = 1.0,
            contact.dependent = F, contact.range = NULL, nboot = config$nboot)

print("Memory usage after computeCommunProb:")
print(pryr::mem_used())

cellchat <- filterCommunication(cellchat, min.cells = config$filter_communication_cells)

print("Memory usage after filterCommunication:")
print(pryr::mem_used())

# extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat) # set slot.name = "netP" to access teh inferred communications at the level of signalling pathways
# infer the cell-cell communication at a signalling pathway leel
cellchat <- computeCommunProbPathway(cellchat)

print("Memory usage after computeCommunProbPathway:")
print(pryr::mem_used())

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)


# visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))

###### Global cell-cell communication network plots ###### 
# global cellchat plots + individual factor cell interaction weight/strength plots
source("resources/cellchat/global_cellchat_plots.R")
plot_cellchat_global(cellchat, output_path)


###### Pathway analysis  & plots ###### 
source("resources/cellchat/pathway_plots.R")
pathway_plots(cellchat, output_path)

###### LR plots for pathways of interest ######
source("resources/cellchat/LR_pair_plots.R")
process_pathway_LR_pairs(cellchat, "COMPLEMENT", output_path)


### Gene expression analysis

# save the cellchat object -> this can be used to load in a different script or in R interactive session to run the cellchat app
saveRDS(cellchat, file.path(output_path, "cellchat_object.rds"))

# run the cellchat app 

#runCellChatApp(cellchat)


