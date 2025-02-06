#' Makes proportion plots for a specified gene
#'
#' Makes a set of proportion plots for transcripts within a specified gene.
#' For use within a plotting wrapper function.
#'
#' @importFrom dplyr {{}}
#' @param prop_table The output of `generate_prop_tables`.
#' @param gene_of_interest Transcripts belonging to this gene will be plotted.
#' @param output_folder_gene A path to the folder to save the plots to.
#' @param plot_width,plot_height,plot_dpi Plotting arguments supplied to `ggplot2::ggsave()`. Plots by default
#' are saved with a width and height of 12, at 600 DPI.
#' @return A list containing a ggplot2 plot objects.

make_plots_proportion <- function(prop_table, gene_of_interest, output_folder_gene, plot_width = 12, plot_height = 12, plot_dpi = 600) {
    
    # Check inputs
    # Check that prop_table is a list that includes AT LEAST 2 clusters (2 elements)
    # Check that gene_of_interest exists in gene col
    # Check that output folder exists
    
    # Define discrete colors here for pie chart
    # The following is a mashup of the Set3, Dark2, and Spectral palettes from RColorBrewer (in that order)
    # If there's more transcripts in the gene than colors here, we'll use the viridis discrete colors instead
    colorvec <- c('#8DD3C7','#BEBADA','#FB8072','#BC80BD','#80B1D3','#FDB462','#B3DE69','#FCCDE5','#D9D9D9','#FFFFB3',
                  '#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D','#666666',
                  '#9E0142','#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF','#E6F598','#ABDDA4','#66C2A5','#3288BD','#5E4FA2')
    
    # Make all folders needed here
    # Folders made are of the form: geneID/plotType
    cat('Plotting started\n')
    
    # Pull the gene of interest and reorganise data
    clusters <- names(prop_table)
    for (cluster in clusters) {
        prop_table[[cluster]] <- prop_table[[cluster]][prop_table[[cluster]]$gene_id == gene_of_interest, ]
    }
    
    plotlist <- list()
    
    # Helper plotting functions
    plotting_function <- function(datatab, x_str, y_str, fill_str, title_str) {
        produced_plot <- ggplot2::ggplot(datatab, ggplot2::aes_string(x = x_str, 
                                                                      y = y_str, 
                                                                      fill = fill_str)) +
            ggplot2::geom_col(position = ggplot2::position_dodge2()) +
            ggplot2::xlab('') +
            ggplot2::ylab('') +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = 'none',
                           title = ggplot2::element_text(face = 'bold', size = 12),
                           axis.text = ggplot2::element_text(size = 12),
                           axis.text.x = ggplot2::element_blank(),
                           panel.grid.major.x = ggplot2::element_blank(),
                           axis.ticks = ggplot2::element_blank(),
                           panel.border = ggplot2::element_blank()) +
            ggplot2::ggtitle(title_str)
        return(produced_plot)
    }
    
    # Plot 1: Transcript counts between clusters (absolute counts)
    # One plot per transcript that's arranged into a single overall image
    cat('\tPlotting transcript counts...\n')
    tr_clusters <- prop_table[[clusters[1]]][, 1:2]
    for (cluster in clusters) {
        tr_clusters[[cluster]] <- prop_table[[cluster]][, 6]
    }
    tr_count_plotlist <- list()
    for (x in tr_clusters$transcript_id) {
        temptab <- tr_clusters[tr_clusters$transcript_id == x, ]
        temptab <- suppressMessages(reshape2::melt(temptab, ids = c('transcript_id', 'gene_id')))
        colnames(temptab) <- c('transcript_id', 'gene_id', 'cluster', 'transcript_count')
        temptab$cluster <- factor(temptab$cluster, levels = levels(temptab$cluster)[order(levels(temptab$cluster))])
        tempplot <- plotting_function(temptab, 'cluster', 'transcript_count', 'cluster', x)
        tr_count_plotlist[[x]] <- tempplot
    }
    tr_count_plot <- ggpubr::ggarrange(plotlist = tr_count_plotlist, common.legend = T, legend = 'bottom')
    plotlist$transcript_counts_individual <- tr_count_plotlist
    plotlist$transcript_counts <- tr_count_plot
    
    # Plot 2: Proportion of gene counts within a cluster that the transcript occupies (transcript proportion within cluster)
    # One plot per transcript that's arranged into a single overall image
    # Note that all props FOR A CLUSTER add up to 1, but props per transcript don't.
    cat('\tPlotting transcript proportions...\n')
    trprop_clusters <- prop_table[[clusters[1]]][, 1:2]
    for (cluster in clusters) {
        trprop_clusters[[cluster]] <- prop_table[[cluster]][, 8]
    }
    tr_prop_plotlist <- list()
    for (x in trprop_clusters$transcript_id) {
        temptab <- trprop_clusters[trprop_clusters$transcript_id == x, ]
        temptab <- suppressMessages(reshape2::melt(temptab, ids = c('transcript_id', 'gene_id')))
        colnames(temptab) <- c('transcript_id', 'gene_id', 'cluster', 'transcript_proportion')
        temptab$cluster <- factor(temptab$cluster, levels = levels(temptab$cluster)[order(levels(temptab$cluster))])
        tempplot <- plotting_function(temptab, 'cluster', 'transcript_proportion', 'cluster', x)
        tr_prop_plotlist[[x]] <- tempplot
    }
    tr_prop_plot <- ggpubr::ggarrange(plotlist = tr_prop_plotlist, common.legend = T, legend = 'bottom')
    plotlist$transcript_proportions_individual <- tr_prop_plotlist
    plotlist$transcript_proportions <- tr_prop_plot
    
    # Plot 3: Proportion of gene counts within a cluster that the transcript occupies (transcript proportion within cluster) ALTERNATE VIEW
    # One plot per CLUSTER (instead of transcript) that's arranged into a single overall image
    # An alternate view of plot 2. Might make more sense since props in every plot add up to 1.
    # This means we can make this plot a pie chart (alt) and a stacked bar (alt2) instead
    cat('\tPlotting transcript proportions within cluster (alt)...\n')
    tr_prop_plotlist_alt <- list()
    for (x in clusters[order(clusters)]) {
        temptab <- prop_table[[x]][, c(1, 2, 8)]
        temptab$transcript_id <- factor(temptab$transcript_id, levels = temptab$transcript_id[order(temptab$transcript_id)])
        tempplot <- ggplot2::ggplot(temptab, ggplot2::aes(x = '',
                                                          y = proportion_in_cluster,
                                                          fill = transcript_id)) +
            ggplot2::geom_bar(stat = 'identity', width = 1) +
            ggplot2::coord_polar(theta = 'y', start = 0) +
            ggplot2::xlab('') +
            ggplot2::ylab('') +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = 'none',
                           title = ggplot2::element_text(face = 'bold', size = 16),
                           axis.text = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           axis.ticks = ggplot2::element_blank(),
                           panel.border = ggplot2::element_blank(),
                           plot.margin = ggplot2::unit(c(0,0.2,0,1), 'lines')) +
            ggplot2::ggtitle(x)
        
        if (nrow(temptab) > 29) {
            tempplot <- tempplot + ggplot2::scale_color_viridis_d(option = 'turbo')
        } else {
            tempplot <- tempplot + ggplot2::scale_fill_manual(values = colorvec[1:nrow(temptab)])
        }
        tr_prop_plotlist_alt[[x]] <- tempplot
    }
    tr_prop_plot_alt <- ggpubr::ggarrange(plotlist = tr_prop_plotlist_alt, common.legend = T, legend = 'bottom')
    plotlist$transcript_proportions_alt_individual <- tr_prop_plotlist_alt
    plotlist$transcript_proportions_alt <- tr_prop_plot_alt
    
    cat('\tPlotting transcript proportions within cluster (alt2)...\n')
    temptablist <- list()
    for (x in clusters[order(clusters)]) {
        temptablist[[x]] <- prop_table[[x]][, c(1, 8)]
        temptablist[[x]]$cluster <- x
    }
    temptab <- data.table::rbindlist(temptablist)
    tempplot <- ggplot2::ggplot(temptab, ggplot2::aes(x = cluster,
                                                      y = proportion_in_cluster,
                                                      fill = transcript_id)) +
        ggplot2::geom_bar(stat = 'identity') +
        ggplot2::xlab('') +
        ggplot2::ylab('') +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = 'bottom',
                       axis.text = ggplot2::element_text(size = 16),
                       axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1),
                       panel.grid.major.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       plot.margin = ggplot2::unit(c(0,0.2,0,1), 'lines'))
    
    plotlist$transcript_proportions_alt2 <- tempplot
    
    # TODO: Need to fix this to work with repeats
    # if (nrow(temptab) > 29) {
    #     tempplot <- tempplot + ggplot2::scale_color_viridis_d(option = 'turbo')
    # } else {
    #     tempplot <- tempplot + ggplot2::scale_fill_manual(values = colorvec[1:nrow(temptab)])
    # }
    
    
    # Save all summary plot images
    cat('\tSaving to folder...\n')
    ggplot2::ggsave(file.path(output_folder_gene, 'transcript_counts_between_clusters.png'), plotlist$transcript_counts,
                    device = 'png', width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
    ggplot2::ggsave(file.path(output_folder_gene, 'transcript_proportions.png'), plotlist$transcript_proportions,
                    device = 'png', width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
    ggplot2::ggsave(file.path(output_folder_gene, 'transcript_proportions_per_cluster_piechart.png'), plotlist$transcript_proportions_alt,
                    device = 'png', width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
    ggplot2::ggsave(file.path(output_folder_gene, 'transcript_proportions_per_cluster_columns.png'), plotlist$transcript_proportions_alt2,
                    device = 'png', width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
    
    # End
    cat('\tFinished plotting!\n')
    return(plotlist)
    
}

