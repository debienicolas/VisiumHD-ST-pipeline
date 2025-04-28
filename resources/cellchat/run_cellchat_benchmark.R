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
library(profmem)
library(peakRAM)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 148000 * 1024^2)  # Set to 128GB

gc()  # Force garbage collection before starting

## set the input parameters
option_list <- list(
    optparse::make_option("--input_path", type="character", default=NULL, help="Path to the annotated anndata object."),
    optparse::make_option("--output_path", type="character", default=NULL, help="Path to the output folder"),
    optparse::make_option("--data_dir", type="character", default=NULL, help="Path to the data directory"),
    optparse::make_option("--memory_file", type="character", default="cellchat_memory.txt", help="Path to write memory usage"),
    optparse::make_option("--run_id", type="character", default=NULL, help="Unique run identifier")
)

# Parse the command line arguments
opt_parser <- optparse::OptionParser(option_list=option_list)
params <- optparse::parse_args(opt_parser)

# Record initial memory
initial_mem <- pryr::mem_used()
print(paste0("Initial memory: ", initial_mem / (1024 * 1024), " MB"))

# Validate required arguments
if (is.null(params$input_path) || is.null(params$output_path) || is.null(params$data_dir)) {
    stop("All arguments (--input_path, --output_path, --data_dir) are required")
}

# Print parsed arguments for debugging
print("Parsed arguments:")
print(params)

output_path <- params$output_path
data_dir <- params$data_dir

# if it doesn't exist, create it
if (!dir.exists(output_path)) {
    dir.create(output_path)
}

# Replace peakRAM tracking with a simpler memory tracking approach
memory_log <- list()

# Function to get current memory stats
get_memory_stats <- function(stage_name) {
    gc()  # Force garbage collection
    gc_stats <- gc()
    mem_used <- pryr::mem_used()
    
    # Get peak memory from gc stats
    peak_mem <- gc_stats[2,6]  # "max used" from Vcells row
    
    print(paste0("=== Memory stats for stage: ", stage_name, " ==="))
    print(paste0("Memory used: ", mem_used / (1024 * 1024), " MB"))
    print(paste0("Peak memory: ", peak_mem, " MB"))
    print("GC Stats:")
    print(gc_stats)
    
    return(list(
        stage = stage_name,
        mem_used = mem_used,
        peak_mem = peak_mem * 1024 * 1024  # Convert to bytes
    ))
}

# Initial memory measurement
memory_log[["initial"]] <- get_memory_stats("initial")

data <- schard::h5ad2seurat_spatial(params$input_path)
data <- NormalizeData(data, normalization.method = "LogNormalize")
# filter data based on counts
data <- subset(data, subset = nFeature_Spatial > 10)
print("Seurat data object: ")
print(data)
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
spot.size = 65  # theoretical spot size (um)
# Read scalefactors from json file
scalefactors <- jsonlite::fromJSON(txt = file.path(data_dir, "binned_outputs/square_008um/spatial/scalefactors_json.json"))

# Calculate conversion factor using binned data parameters
# The spot diameter needs to be adjusted based on the bin size
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

# After loading data
memory_log[["after_data_load"]] <- get_memory_stats("after_data_load")

#### Load the seurat object from .h5ad file using schard
data <- schard::h5ad2seurat_spatial(params$input_path)
memory_log[["after_seurat_load"]] <- get_memory_stats("after_seurat_load")

# After normalization
data <- NormalizeData(data, normalization.method = "LogNormalize")
memory_log[["after_normalization"]] <- get_memory_stats("after_normalization")

# After filtering
data <- subset(data, subset = nFeature_Spatial > 10)
memory_log[["after_filtering"]] <- get_memory_stats("after_filtering")

# After CellChat creation
cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = spatial.locs,
    spatial.factors = spatial.factors
)
memory_log[["after_cellchat_creation"]] <- get_memory_stats("after_cellchat_creation")

# After setting DB
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
memory_log[["after_db_setup"]] <- get_memory_stats("after_db_setup")

# After gene expression analysis
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
memory_log[["after_gene_analysis"]] <- get_memory_stats("after_gene_analysis")

# After network inference
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
            distance.use = TRUE, interaction.range = 300, scale.distance = 1.0,
            contact.dependent = F, contact.range = NULL, nboot = 10)
memory_log[["after_network_inference"]] <- get_memory_stats("after_network_inference")

# After communication filtering
cellchat <- filterCommunication(cellchat, min.cells = 20)
memory_log[["after_filtering_communication"]] <- get_memory_stats("after_filtering_communication")

# Final memory measurement
memory_log[["final"]] <- get_memory_stats("final")

# Create a summary of memory usage
memory_summary <- data.frame(
    Run_ID = params$run_id,
    Stage = sapply(memory_log, function(x) x$stage),
    Memory_Used_MB = sapply(memory_log, function(x) x$mem_used / (1024 * 1024)),
    Peak_RAM_MB = sapply(memory_log, function(x) x$peak_mem / (1024 * 1024))
)

# Ensure output directory exists
if (!dir.exists(params$output_path)) {
    dir.create(params$output_path, recursive = TRUE)
}

# Write detailed memory log
write.csv(memory_summary, file = file.path(params$output_path, "memory_profile.csv"))

# Ensure memory file directory exists
memory_file_dir <- dirname(params$memory_file)
if (!dir.exists(memory_file_dir)) {
    dir.create(memory_file_dir, recursive = TRUE)
}

# Write peak memory to the memory file for the Python script
peak_memory <- max(memory_summary$Peak_RAM_MB)
tryCatch({
    write(as.character(peak_memory), file = params$memory_file)
    print(paste0("Wrote peak memory to file: ", params$memory_file))
}, error = function(e) {
    print(paste0("Error writing to memory file: ", e$message))
})

# Print summary
print("Memory Usage Summary:")
print(memory_summary)
print(paste0("Peak RAM usage: ", peak_memory, " MB"))

### Gene expression analysis

# save the cellchat object -> this can be used to load in a different script or in R interactive session to run the cellchat app
saveRDS(cellchat, file.path(output_path, "cellchat_object.rds"))

# run the cellchat app 

#runCellChatApp(cellchat)

# At the end of the script, just before saving:
final_mem <- pryr::mem_used()
memory_used <- (final_mem - initial_mem) / (1024 * 1024)  # Convert to MB
print(paste0("Final memory: ", final_mem / (1024 * 1024), " MB"))
print(paste0("Memory used: ", memory_used, " MB"))

# At the end, ensure the memory file is written
tryCatch({
    write(as.character(memory_used), file = params$memory_file)
    print(paste0("Wrote final memory usage to file: ", params$memory_file))
}, error = function(e) {
    print(paste0("Error writing final memory to file: ", e$message))
})


