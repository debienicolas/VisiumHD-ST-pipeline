library(CellChat)


pathway_plots <- function(cellchat, output_path) {
    # initialize variables
    all_pathways <- slot(cellchat, "netP")$pathway
    message("Plotting all pathways...")
    cat("Pathways:", paste(all_pathways, collapse = ", "), "\n")
    n_celltypes <- length(cellchat@idents)

    # compute the network centrality scores
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

    # create a subdirectory for the pathway plots
    pathway_plots_dir <- file.path(output_path, "pathway_plots")
    if (!dir.exists(pathway_plots_dir)) {
        dir.create(pathway_plots_dir)
    }

    # loop over each pathway and create a single pdf file with all the plots for that pathway
    for (pathway in all_pathways) {
        message("Plotting pathway: ", pathway)

        pdf(file.path(pathway_plots_dir, paste0(pathway, ".pdf")), width = 20, height = 16)
        plot1 <- netVisual_aggregate(cellchat, signaling = pathway, layout = "spatial", edge.width.max = 2, vertex.size.max = 6, alpha.image = 0.025, vertex.label.cex = 0.0, vertex.weight = "incoming")
        print(plot1) 
        message("Plot 1 done!")
        plot2 <- netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
        print(plot2)
        message("Plot 2 done!")
        # takes a very long time to calculate for every incoming and outgoing cell type
        #plot3 <- netVisual_bubble(cellchat, sources.use = c(1:n_celltypes), targets.use = c(1:n_celltypes), signaling = pathway, remove.isolate = FALSE)   
        #print(plot3)
        #message("Plot 3 done!")
        plot4 <- plotGeneExpression(cellchat, signaling = pathway)
        print(plot4)   
        message("Plot 4 done!")
        plot5 <- netAnalysis_signalingRole_network(cellchat, signaling = pathway, width = 10, height = 8)
        print(plot5)
        message("Plot 5 done!")
        # LR contribution bar chart
        plot6 <- netAnalysis_contribution(cellchat, signaling = pathway)
        print(plot6)
        message("Plot 6 done!")
        
        dev.off()
    }

    message("Done plotting all pathways!")


}

