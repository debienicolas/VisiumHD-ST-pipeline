# Load required libraries
library(CellChat)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(ggplotify)

# Function to create all possible CellChat plots for a specific L-R pair
visualize_LR_pair <- function(cellchat, pathway.name, LR.pair, output_dir) {
    
    # Create PDF file
    filename <- file.path(output_dir, paste0(gsub("[/:*?\"<>|]", "_", LR.pair), "_analysis.pdf"))
    pdf(filename, width = 20, height = 16)
    
    # Create and print plots directly
    tryCatch({
        # Circle plot - counts
        netVisual_individual(cellchat,signaling = pathway.name, pairLR.use = LR.pair,
                           layout = "circle",
                           height = 10)
        message("Circle plot - counts: DONE")

        # Convert LR.pair to a data frame with the correct column name
        LR_pair_df <- data.frame(interaction_name = LR.pair)
        n_celltypes <- length(cellchat@idents)
        print(netVisual_bubble(cellchat, 
                             sources.use = c(1:n_celltypes), 
                             targets.use = c(1:n_celltypes), 
                             pairLR.use = LR_pair_df, 
                             remove.isolate = FALSE))
        message("Netvisual bubble: DONE")

        # Gene expression -> doesn't accept pairLR.use
        #plotGeneExpression(cellchat, signaling = pathway.name, pairLR.use = LR.pair)
        #message("Gene expression")
        # spatial feature plot
        spatialPlot <- spatialFeaturePlot(cellchat, pairLR.use = LR.pair, point.size = 0.1, 
                            do.binary = FALSE, cutoff = 0.5, enriched.only = F, 
                            color.heatmap = "Reds", direction = 1)
        print(spatialPlot)
        message("Spatial feature plot: DONE")

        # spatial feature plot binary
        spatialPlotBinary <- spatialFeaturePlot(cellchat, pairLR.use = LR.pair, point.size = 0.3, 
                            do.binary = TRUE, cutoff = 0.5, enriched.only = F, 
                            color.heatmap = "Reds", direction = 1)
        print(spatialPlotBinary)
        message("Spatial feature plot binary: DONE")
        
        
    }, error = function(e) {
        message(sprintf("Warning: Error creating plots for %s: %s", LR.pair, e$message))
    })
    
    # Close the PDF device
    dev.off()
    
    message(sprintf("Created plots for %s", LR.pair))
}

# Main function to process all L-R pairs in a pathway
process_pathway_LR_pairs <- function(cellchat, pathway.name, output_dir) {
    output_dir <- file.path(output_dir, paste0(pathway.name, "_LR_plots"))
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }

    # Get all L-R pairs for the pathway
    LR.pairs <- extractEnrichedLR(cellchat, signaling = pathway.name, geneLR.return = FALSE)
    
    # Process each L-R pair
    for (i in 1:nrow(LR.pairs)) {
        # Create plots for this L-R pair
        LR.show <- LR.pairs[i,]
        visualize_LR_pair(cellchat, pathway.name, LR.show, output_dir)
        
        # Print progress
        message(sprintf("Processed %d/%d L-R pairs", i, nrow(LR.pairs)))
    }
    
    message(sprintf("Complete! PDFs saved to %s", output_dir))
}

