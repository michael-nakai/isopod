#' Make and save DIU plots
#'
#' A wrapper function around multiple plotting functions.
#' First runs `generate_prop_tables`, then creates DTU plots, proportion plots, and UMAPs.
#'
#' @importFrom magrittr %>%
#' @param transcript_counts_table A dataframe with a column containing transcript IDs, a column with gene IDs for
#' each transcript, and counts for all transcripts in all cells, with cell IDs as subsequent column names.
#' @param cell_labels_table A dataframe with two columns: one listing all cell IDs, and the other listing the group
#' that the cell ID belongs to.
#' @param transcript_id_colname A string corresponding to the column name in transcript_counts_table where the transcript 
#' IDs are stored. Defaults to 'transcript_id'.
#' @param gene_id_colname A string corresponding to the column name in transcript_counts_table where the gene IDs are 
#' stored. Defaults to 'gene_id'.
#' @param cell_labels_colname A string corresponding to the column name in cell_labels_table where the group information 
#' is stored. Defaults to 'group'.
#' @param pvalue_object The p-value object created from `get_permutation_pvals`.
#' @param transcript_id_colname A string corresponding to the column name in transcript_counts_table where the 
#' transcript IDs are stored.
#' @param gene_of_interest The name of a gene of interest to create additional plots for. If set, proportion UMAPs
#' and transcript proportion plots are created for the gene. Defaults to NA.
#' @param generate_UMAPs A bool determining whether UMAPs should be calculated/plotted. UMAP generation can take 
#' a while for larger datasets, so the option is provided to skip this step if computational resources are limited.
#' Default TRUE.
#' @param save_to_folder A string corresponding to the filepath to the folder that all outputs will be saved to. 
#' Defaults to NA.
#' @return A list containing the plots generated, as well as the data used to create these plots.
#' @export

make_plots <- function(transcript_counts_table, cell_labels_table,
                       transcript_id_colname, gene_id_colname,
                       cell_labels_colname, pvalue_object, 
                       gene_of_interest = NA,
                       generate_UMAPs = TRUE, save_to_folder = NA) {
    
    # Output structure should be within whatever folder is specified.
    # Create the dirs plot_data and plots
    # ggplot2 RDS objects go in plot_data, images are saved into plot_images, 
    # all other intermediate data goes into table_data
    if (!is.na(save_to_folder)) {
        # Shortcuts
        plot_images_folder <- file.path(save_to_folder, 'plot_images')
        plot_data_folder <- file.path(save_to_folder, 'plot_data')
        table_data_folder <- file.path(save_to_folder, 'table_data')
        
        # Create dirs
        dir.create(file.path(save_to_folder), showWarnings = FALSE)
        dir.create(file.path(table_data_folder), showWarnings = FALSE)
        dir.create(file.path(plot_data_folder), showWarnings = FALSE)
        dir.create(file.path(plot_images_folder), showWarnings = FALSE)
        if (!is.na(gene_of_interest)) {
            table_data_folder_gene <- file.path(table_data_folder, gene_of_interest)
            plot_data_folder_gene <- file.path(plot_data_folder, gene_of_interest)
            plot_images_folder_gene <- file.path(plot_images_folder, gene_of_interest)
            dir.create(table_data_folder_gene, showWarnings = FALSE)
            dir.create(plot_data_folder_gene, showWarnings = FALSE)
            dir.create(plot_images_folder_gene, showWarnings = FALSE)
        }
    }
    
    master_plot_list <- list()
    master_table_list <- list()
    
    # Generate proportion tables to pass to plotting functions
    cat('Creating proportion metadata...\n')
    proportions_tables <- generate_prop_tables(transcript_counts_table, cell_labels_table,
                                               transcript_id_colname, gene_id_colname,
                                               cell_labels_colname)
    
    if (!is.na(save_to_folder)) {
        saveRDS(proportions_tables, file.path(table_data_folder, 'proportion_data.rds'))
    }
    
    cat('Starting plotting...\n')
    
    # Helper for replacing div by 0 nans with 0s
    # Credit to Hong Ooi (https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame)
    is.nan.data.frame <- function(x) {
        do.call(cbind, lapply(x, is.nan))
    }
    
    # Generate UMAPs
    if (generate_UMAPs) {
        
        
        if (!is.na(save_to_folder)) {
            dir.create(file.path(plot_data_folder, 'UMAPs'), showWarnings = FALSE)
        }
        
        cat('\tGenerating UMAP...\n')
        cols_to_remove <- match(c(transcript_id_colname, gene_id_colname), colnames(transcript_counts_table))
        cells_by_transcripts <- transcript_counts_table[, -cols_to_remove]
        cells_by_transcripts <- as.data.frame(t(cells_by_transcripts))
        cell_id_vec <- rownames(cells_by_transcripts)
        rownames(cells_by_transcripts) <- NULL
        cells_by_transcripts <- as.matrix(cells_by_transcripts)
        umap_tab <- umap::umap(cells_by_transcripts)
        umap_coords <- as.data.frame(umap_tab$layout)
        colnames(umap_coords) <- c('x', 'y')
        umap_coords$cell_id <- cell_id_vec
        cell_id_colname <- colnames(cell_labels_table)[colnames(cell_labels_table) != cell_labels_colname]
        umap_coords <- merge(umap_coords, cell_labels_table, by.x = 'cell_id', by.y = cell_id_colname)
        group_col <- colnames(umap_coords)[!(colnames(umap_coords) %in% c('cell_id', 'x', 'y'))]
        
        # Plot a normal UMAP
        normal_umap <- ggplot2::ggplot(umap_coords, ggplot2::aes_string(x = 'x', y = 'y', color = group_col)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        
        master_table_list$UMAP <- umap_coords
        master_plot_list$UMAP <- normal_umap
        
        if (!is.na(save_to_folder)) {
            saveRDS(umap_coords, file.path(table_data_folder, 'UMAP.rds'))
            saveRDS(normal_umap, file.path(plot_data_folder, 'UMAP.rds'))
            ggplot2::ggsave(file.path(plot_images_folder, 'UMAP.png'), normal_umap,
                            device = 'png', width = 8, height = 8, dpi = 600)
        }
        
        # Calculate transcript props per cell for all transcripts in the gene of interest (if set)
        # Analogous to Feature Plots from Seurat
        if (!is.na(gene_of_interest)) {
            cat('\tConstructing feature plots for', gene_of_interest, '...\n')
            cut_table <- transcript_counts_table[transcript_counts_table[[gene_id_colname]] == gene_of_interest, ]
            cut_table[[gene_id_colname]] <- NULL
            cut_transcripts <- cut_table[[transcript_id_colname]]
            cut_table[[transcript_id_colname]] <- NULL
            cut_cells <- colnames(cut_table)
            cut_table <- as.data.frame(t(cut_table))
            colnames(cut_table) <- cut_transcripts
            total_counts_per_cell <- rowSums(cut_table)
            cut_table_props <- cut_table / total_counts_per_cell  # Row-wise division of tx counts in cell / total gene counts in cell
            tids <- colnames(cut_table_props)
            tids <- tids[order(tids)]
            cut_table$cell_id <- rownames(cut_table)
            cut_table_props$cell_id <- rownames(cut_table_props)
            rownames(cut_table) <- NULL
            rownames(cut_table_props) <- NULL
            cut_table <- merge(umap_coords, cut_table, by = 'cell_id')
            cut_table_props <- merge(umap_coords, cut_table_props, by = 'cell_id')
            cut_table_props[is.nan(cut_table_props)] <- 0
            featureplot_list <- list()
            propplot_list <- list()
            
            for (tid in tids) {
                temptab <- cut_table[order(cut_table[[tid]]), ]
                countplot <- ggplot2::ggplot(temptab, 
                                             ggplot2::aes_string(x = 'x', y = 'y', color = tid)) +
                    ggplot2::geom_point() +
                    ggplot2::ggtitle(tid) +
                    ggplot2::xlab('UMAP 1') +
                    ggplot2::ylab('UMAP 2') +
                    ggplot2::theme(text = ggplot2::element_text(size = 12),
                                   title = ggplot2::element_text(face = 'bold', size = 12),
                                   axis.text = ggplot2::element_text(size = 12),
                                   axis.ticks = ggplot2::element_blank(),
                                   legend.title = ggplot2::element_blank(),
                                   panel.grid.minor = ggplot2::element_blank())
                featureplot_list[[tid]] <- countplot
                
                temptab <- cut_table_props[order(cut_table_props[[tid]]), ]
                propplot <- ggplot2::ggplot(temptab, 
                                             ggplot2::aes_string(x = 'x', y = 'y', color = tid)) +
                    ggplot2::geom_point() +
                    ggplot2::ggtitle(tid) +
                    ggplot2::xlab('UMAP 1') +
                    ggplot2::ylab('UMAP 2') +
                    ggplot2::theme(text = ggplot2::element_text(size = 12),
                                   title = ggplot2::element_text(face = 'bold', size = 12),
                                   axis.text = ggplot2::element_text(size = 12),
                                   axis.ticks = ggplot2::element_blank(),
                                   legend.title = ggplot2::element_blank(),
                                   panel.grid.minor = ggplot2::element_blank())
                propplot_list[[tid]] <- propplot
            }
            
            # Arrange into combined plots
            featureplot_all <- ggpubr::ggarrange(plotlist = featureplot_list)
            propplot_all <- ggpubr::ggarrange(plotlist = propplot_list)
            
            master_table_list[[paste0('umap_counts_', gene_of_interest)]] <- cut_table
            master_table_list[[paste0('umap_proportions_', gene_of_interest)]] <- cut_table_props
            master_plot_list[[paste0('umap_counts_', gene_of_interest)]] <- featureplot_all
            master_plot_list[[paste0('umap_counts_', gene_of_interest, '_individual')]] <- featureplot_list
            master_plot_list[[paste0('umap_proportions_', gene_of_interest)]] <- propplot_all
            master_plot_list[[paste0('umap_proportions_', gene_of_interest, '_individual')]] <- propplot_list
            
            if (!is.na(save_to_folder)) {
                # Plot data
                saveRDS(cut_table, file.path(table_data_folder_gene, 'UMAP_tx_counts_table.rds'))
                saveRDS(cut_table_props, file.path(table_data_folder_gene, 'UMAP_tx_props_table.rds'))
                saveRDS(featureplot_all, file.path(plot_data_folder_gene, 'UMAP_tx_counts_plot.rds'))
                saveRDS(featureplot_list, file.path(plot_data_folder_gene, 'UMAP_tx_counts_individual_plots.rds'))
                saveRDS(propplot_all, file.path(plot_data_folder_gene, 'UMAP_tx_cell_props_plot.rds'))
                saveRDS(propplot_list, file.path(plot_data_folder_gene, 'UMAP_tx_cell_props_individual_plots.rds'))
                
                # Combined plots
                ggplot2::ggsave(file.path(plot_images_folder_gene, 'UMAP_tx_counts.png'), featureplot_all,
                                device = 'png', width = 14, height = 14, dpi = 600, bg = "white")
                ggplot2::ggsave(file.path(plot_images_folder_gene, 'UMAP_tx_cell_props.png'), propplot_all,
                                device = 'png', width = 14, height = 14, dpi = 600, bg = "white")
            }
        }
    }
    
    # DTU MA-like plot generation
    # 1 per group
    plot_list <- list()
    table_list <- list()
    groups_assessed <- colnames(pvalue_object$permutation_pvalues)
    groups_assessed <- groups_assessed[3:length(groups_assessed)]
    for (n in groups_assessed) {
        cat('\tPlotting transcript usage for', n,'\n')
        contingency_lists <- pvalue_object[['first-loop_contingency_tables']]
        pvalue_table <- pvalue_object[['permutation_pvalues']]
        pval_to_merge <- pvalue_table %>% dplyr::select(transcript_id)
        pval_to_merge <- cbind(pval_to_merge, pvalue_table[[n]])
        colnames(pval_to_merge) <- c(transcript_id_colname, 'permutation_pval')
        contingency_to_merge <- contingency_lists[[n]] %>% dplyr::select({{transcript_id_colname}}, isoform_in, isoform_out, other_in, other_out)

        temp_df <- merge(proportions_tables[[n]], pval_to_merge, by = transcript_id_colname)
        moretemp <-  pvalue_object$`first-loop_pvalues`
        moretemp <- moretemp[ order(match(moretemp[[transcript_id_colname]], temp_df[[transcript_id_colname]])), ]
        temp_df$firstloop_pval <- moretemp[[n]]
        contingency_to_merge <- merge(temp_df, contingency_to_merge, by = transcript_id_colname)
        temp_df <- contingency_to_merge

        # Remove any isoforms with nt == 0 & rt == 0
        temp_df <- temp_df[which((temp_df$isoform_in > 0) | (temp_df$isoform_out > 0)), ]

        # Label pval_sig
        pval_sig <- rep(NA, nrow(temp_df))
        i <- 1
        for (v in temp_df$permutation_pval) {
            if (v < 0.05) {
                pval_sig[i] <- 'p < 0.05'
            } else {
                pval_sig[i] <- 'Not significant'
            }
            i <- i + 1
        }
        temp_df$pval_significance_perm <- pval_sig

        grapht <- paste0(n, ' transcript usage, significant transcripts highlighted')

        plot_list[[n]] <- ggplot2::ggplot(data = temp_df %>% dplyr::arrange(desc(permutation_pval)), mapping = ggplot2::aes_string(x = 'total_gene_count', y = 'difference_in_proportions', color = 'pval_significance_perm')) +
            ggplot2::geom_point() +
            ggplot2::labs(x = 'Total counts within gene',
                          y = paste0('Difference in proportions for cluster ', n, ' (cluster - average)')) + # Because this is inverted, higher = more proportion in cluster compared to in all cells
            ggplot2::scale_x_continuous(trans = 'log2') +
            ggplot2::ggtitle(grapht)

        table_list[[n]] <- temp_df

        if (!is.na(save_to_folder)) {
            saveRDS(table_list, file.path(table_data_folder, paste0(n, '_proportion_differences.rds')))
            saveRDS(plot_list, file.path(plot_data_folder, paste0(n, '_proportion_differences_plots.rds')))
            ggplot2::ggsave(file.path(plot_images_folder, paste0(n, '_proportion_differences.png')),
                            plot_list[[n]],
                            device = 'png',
                            width = 8,
                            height = 8,
                            dpi = 600)
        }
    }
    
    master_plot_list$proportion_differences <- plot_list
    master_table_list$proportion_differences <- table_list

    return(list('plots' = master_plot_list, 'tables' = master_table_list))
}
