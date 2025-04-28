library(CellChat)
library(dplyr)
library(magrittr)
library(patchwork)
library(extrafont)
library(optparse)
library(Seurat)
library(CellChat)
library(gridExtra)
library(grid)

# function to create cellchat plots over all pathways and celltypes

plot_cellchat_global <- function(cellchat, output_path) {
    # cellchat: cellchat object
    # output_path: path to save the plots (cellchat folder in results folder)

    # create a subfolder where the individual plots will be saved
    if (!dir.exists(file.path(output_path, "global_plots"))) {
        dir.create(file.path(output_path, "global_plots"))
    }

    # Create a circle plot with interaction counts
    png(file.path(output_path, "global_plots", "circle_interaction_counts.png"), width = 2000, height = 2200, res = 200)
    par(mar = c(2, 2, 4, 2), xpd = TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), 
                    weight.scale = T, label.edge= F, title.name = "Number of interactions",
                    vertex.label.cex = 0.8)
    dev.off()

    # Create a circle plot with interaction weights
    png(file.path(output_path, "global_plots", "circle_interaction_weights.png"), width = 2000, height = 2200, res = 200)
    par(mar = c(2, 2, 4, 2), xpd = TRUE)
    netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), 
                    weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",
                    vertex.label.cex = 0.8)
    dev.off()

    # Create the heatmaps
    png(file.path(output_path, "global_plots", "heatmap_interaction_counts.png"), width = 2000, height = 2200, res = 200)
    print(netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues"))
    dev.off()

    png(file.path(output_path, "global_plots", "heatmap_interaction_weights.png"), width = 2000, height = 2200, res = 200)
    print(netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues"))
    dev.off()

    # signalling role scatter plot
    # png(file.path(output_path, "global_plots", "signalling_role_scatter.png"), width = 2000, height = 2200, res = 200)
    # print(netAnalysis_signalingRole_scatter(cellchat))
    # dev.off()


    message("Combining all plots into a single png...")
    #### Combine all plots into a single png ###
    plots <- lapply(c("circle_interaction_counts.png", "circle_interaction_weights.png", 
                    "heatmap_interaction_counts.png", "heatmap_interaction_weights.png"), 
                    function(x) {
                        img <- png::readPNG(file.path(output_path, "global_plots", x))
                        rasterGrob(img)
                    })
    png(file.path(output_path, "cellchat_global_plots.png"), 
        width = 4000, height = 4000, res = 200)
    grid.arrange(grobs = plots, ncol = 2, 
                padding = unit(2, "line"))
    dev.off()


    ### Individual factor/celltype cell interaction weight/strength plots ###
    message("Creating individual factor/celltype cell interaction weight/strength plots...")
    mat <- cellchat@net$weight
    
    # Create directory for individual plots if it doesn't exist
    individual_plots_dir <- file.path(output_path, "individual_network_plots")
    if (!dir.exists(individual_plots_dir)) {
        dir.create(individual_plots_dir)
    }
    
    # Set up the combined plot
    png(file.path(output_path, "cellchat_network_plots_detailed.png"), 
        width = 4000, height = 4000, res = 200)
    oldpar <- par(no.readonly = TRUE)
    layout(matrix(1:12, nrow=3, ncol=4, byrow=TRUE), widths=rep(1,4), heights=rep(1,3))
    par(mar=c(3,3,4,3), xpd=TRUE)
    
    for (i in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        
        # Save individual plot
        png(file.path(individual_plots_dir, paste0("network_plot_", rownames(mat)[i], ".png")),
            width = 2000, height = 2200, res = 200)
        netVisual_circle(mat2, 
                        vertex.weight = groupSize, 
                        weight.scale = T, 
                        edge.weight.max = max(mat), 
                        title.name = rownames(mat)[i],
                        vertex.label.cex = 0.8)
        dev.off()
        
        # Add to combined plot
        netVisual_circle(mat2, 
                        vertex.weight = groupSize, 
                        weight.scale = T, 
                        edge.weight.max = max(mat), 
                        title.name = rownames(mat)[i],
                        vertex.label.cex = 0.8)
    }
    
    par(oldpar)
    dev.off()

}