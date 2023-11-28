#' Make and save DIU plots
#'
#' Creates a differential isoform usage plot for each group.
#'
#' @param proportions_tables The list of dataframes created using generate_prop_tables().
#' @param pvalue_object The p-value object created from get_permutation_pvals().
#' @param transcript_id_colname A string corresponding to the column name in transcript_counts_table where the transcript IDs are stored.
#' @param save_plots Default NA. If specified, the path to the directory that plots should be saved in.
#' @param save_tables Default NA. If specified, the path to the directory that the tables used to make the plots
#' should be saved in. Useful if you'd like to recreate the plots yourself or in another program.
#' @return A list of ggplot objects, each corresponding to a plot for a cell grouping.
#' @export

make_plots <- function(proportions_tables, pvalue_object, transcript_id_colname, save_plots = NA, save_tables = NA) {

    plot_list <- list()
    table_list <- list()
    for (n in names(proportions_tables)) {

        contingency_lists <- pvalue_object[['first-loop_contingency_tables']]
        pvalue_table <- pvalue_object[['permutation_pvalues']]
        pval_to_merge <- pvalue_table %>% dplyr::select({{transcript_id_colname}})
        pval_to_merge <- cbind(pval_to_merge, pvalue_table[[n]])
        colnames(pval_to_merge) <- c(transcript_id_colname, 'permutation_pval')
        contingency_to_merge <- contingency_lists[[n]] %>% dplyr::select({{transcript_id_colname}}, nt, rt, ng, rg)

        temp_df <- merge(proportions_tables[[n]], pval_to_merge, by = transcript_id_colname)
        temp_df <- merge(temp_df, pvalue_object$`first-loop_pvalues`[[n]], by = transcript_id_colname)
        contingency_to_merge <- merge(temp_df, conteingency_to_merge, by = transcript_id_colname)

        # Remove any isoforms with nt == 0 & rt == 0
        temp_df <- temp_df[which((temp_df$nt > 0) | (temp_df$rt > 0)), ]

        # Calculate pval_sig
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

        grapht <- paste0(n, ' isoform usage, significant isoforms highlighted')

        plot_list[[n]] <- ggplot2::ggplot(data = temp_df %>% dplyr::arrange(desc(permutation_pval)), mapping = ggplot2::aes_string(x = 'total_gene_count', y = 'difference_in_proportions', color = 'pval_significance_perm')) +
            ggplot2::geom_point() +
            ggplot2::labs(x = 'Total counts within gene',
                          y = paste0('Difference in proportions for cluster ', n, ' (cluster - average)')) + # Because this is inverted, higher = more proportion in cluster compared to in all cells
            ggplot2::scale_x_continuous(trans = 'log2') +
            ggplot2::ggtitle(grapht)

        table_list[[n]] <- temp_df

        if (!(is.na(save_tables))) {
            # If the filepath doesn't end with a slash, add it first, then save.
            if (substr(save_tables, nchar(save_tables), nchar(save_tables)) != '/') {
                save_tables <- paste0(save_tables, '/')
            }
            write.table(temp_df,
                        file = paste0(save_plots, n, '.tsv'),
                        sep = '\t',
                        row.names = F,
                        quote = F)
        }

        if (!(is.na(save_plots))) {
            # If the filepath doesn't end with a slash, add it first, then save.
            if (substr(save_plots, nchar(save_plots), nchar(save_plots)) != '/') {
                save_plots <- paste0(save_plots, '/')
            }
            ggsave(paste0(save_plots, n, '.png'),
                   plot_list[[n]],
                   device = 'png',
                   width = 8,
                   height = 8,
                   dpi = 600)
        }
    }

    return(list('plots' = plot_list, 'tables' = table_list))
}
