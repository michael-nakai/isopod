#' Generate a list of transcript proportion tables per group
#'
#' This function generates a list of dataframes that include information regarding each transcript's counts,
#' gene-level counts, proportion within all cells, group-specific proportions, and proportion differences.
#'
#' @importFrom magrittr %>%
#' @param transcript_counts_table A dataframe with a column containing transcript IDs, a column with gene IDs for
#' each transcript, and counts for all transcripts in all cells, with cell IDs as subsequent column names.
#' @param cell_labels_table A dataframe with two columns: one listing all cell IDs, and the other listing the group
#' that the cell ID belongs to.
#' @param transcript_id_colname A string corresponding to the column name in transcript_counts_table where the transcript IDs are stored. Defaults to 'transcript_id'.
#' @param gene_id_colname A string corresponding to the column name in transcript_counts_table where the gene IDs are stored. Defaults to 'gene_id'.
#' @param cell_labels_colname A string corresponding to the column name in cell_labels_table where the group information is stored. Defaults to 'group'.
#' @return A list of dataframes, with each dataframe corresponding to a specific group listed in cell_labels_table.
#' @examples
#' counts_table <- data.frame('transcript_id' = c('t1', 't2', 't3'), 'gene_id' = c('g1', 'g1', 'g2'), 'Cell_1' = c(0, 1, 10), 'Cell_2' = c(10, 2, 5), 'Cell_3' = c(2, 5, 1))
#' labels_table <- data.frame('grouping' = c('Cluster_1', 'Cluster_2', 'Cluster 1'), 'Cells' = c('Cell_1', 'Cell_2', 'Cell_3'))
#' proportion_tables <- generate_prop_tables(counts_table, labels_table, 'transcript_id', 'gene_id', 'grouping')
#' @export


generate_prop_tables <- function(transcript_counts_table, cell_labels_table,
                                 transcript_id_colname='transcript_id', gene_id_colname='gene_id',
                                 cell_labels_colname='group') {
  
    ### Argument sanity checking
    # All strings are ACTUALLY strings
    if (typeof(transcript_id_colname) != 'character' | length(transcript_id_colname) != 1) {
        stop('The transcript_id_colname argument is not a string.')
    }
    if (typeof(gene_id_colname) != 'character' | length(gene_id_colname) != 1) {
        stop('The gene_id_colname argument is not a string.')
    }
    if (typeof(cell_labels_colname) != 'character' | length(cell_labels_colname) != 1) {
        stop('The cell_labels_colname argument is not a string.')
    }
    
    # Tables are full tables
    if (is.null(dim(cell_labels_table))) {
        stop("Your cell_labels_table is an empty dataframe or NULL object.\nPlease fill out your counts table and rerun this function.")
    } else if (dim(cell_labels_table)[1] == 0 | dim(cell_labels_table)[2] == 0) {
        stop("Your cell_labels_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    }
    if (sum(sapply(cell_labels_table, anyNA)) > 0) {
        stop("Your cell_labels_table contains NA values. Please fix these first (either by removing or replacing with zeros) and re-run.")
    }
    if (is.null(dim(cell_labels_table))) {
        stop("Your cell_labels_table is an empty dataframe or NULL object.\nPlease fill out your counts table and rerun this function.")
    } else if (dim(cell_labels_table)[1] == 0 | dim(cell_labels_table)[2] == 0) {
        stop("Your cell_labels_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    }
    if (sum(sapply(cell_labels_table, anyNA)) > 0) {
        stop("Your cell_labels_table contains NA values. Please fix these first (either by removing or replacing with zeros) and re-run.")
    }
    
    # Check that transcript_ and gene_id_colnames exist in colnames(transcript_counts_table), and same for cell_labels_colname
    if (!(transcript_id_colname %in% colnames(transcript_counts_table))) {
        stop(paste0(transcript_id_colname,
                    " could not be found as a column name in your counts table.\nPlease check your column names and rerun this function."))
    } else if (!(gene_id_colname %in% colnames(transcript_counts_table))) {
        stop(paste0(gene_id_colname,
                    " could not be found as a column name in your counts table.\nPlease check your column names and rerun this function."))
    } else if (!(cell_labels_colname %in% colnames(cell_labels_table))) {
        stop(paste0(cell_labels_colname,
                    " could not be found as a column name in your cell labels table.\nPlease check your column names and rerun this function."))
    }
  

    other_colname_for_cell_labels <- grep(cell_labels_colname, colnames(cell_labels_table), invert = TRUE, value = TRUE)
    proplist <- list()
    prop_calc_df <- transcript_counts_table %>% dplyr::select({{transcript_id_colname}}, {{gene_id_colname}})  # For finding proportions later on

    # Find total transcript counts across all clusters
    nolabel_transcript_counts_table <- transcript_counts_table[ , !(colnames(transcript_counts_table) %in% c(transcript_id_colname, gene_id_colname))]
    prop_calc_df[['total_transcript_count']] <- rowSums(nolabel_transcript_counts_table)

    # Find total gene counts across all clusters
    a <- aggregate(prop_calc_df$total_transcript_count, by = list(prop_calc_df$gene_id), FUN=sum)
    colnames(a) <- c(gene_id_colname, 'total_gene_count')
    prop_calc_df <- merge(prop_calc_df, a, by = gene_id_colname, all.x = TRUE)

    # Calculate proportion across all clusters
    prop_calc_df[['all_cells_prop']] <- prop_calc_df$total_transcript_count / prop_calc_df$total_gene_count

    # Add to newdf to be propogated between cluster-specific DFs
    # Sorting first ensures the row ordering is the same between all tables
    prop_calc_df <- prop_calc_df[order(prop_calc_df[[gene_id_colname]], prop_calc_df[[transcript_id_colname]]), ]
    new_transcript_counts_table <- transcript_counts_table[order(transcript_counts_table[[gene_id_colname]], transcript_counts_table[[transcript_id_colname]]), ]
    newdf <- prop_calc_df %>% dplyr::select({{transcript_id_colname}}, {{gene_id_colname}}, total_transcript_count, total_gene_count, all_cells_prop)

    # Make df per cluster, then add to prop_calc_df for later
    for (cluster in unique(cell_labels_table[[cell_labels_colname]])) {

        # Subset transcript counts table to only include cells in cluster
        clusterdf <- new_transcript_counts_table[ , (colnames(new_transcript_counts_table) %in% subset(cell_labels_table,
                                                                                                       cell_labels_table[cell_labels_colname] == cluster)[[other_colname_for_cell_labels]])]

        ### Get cluster-specific stats here first
        proplist_df <- newdf

        # Cluster-specific transcript counts
        if (is.null(nrow(clusterdf))) {
          proplist_df[[paste0(cluster, '_transcript_count')]] <- sum(clusterdf)
          prop_calc_df[[paste0(cluster, '_transcript_count')]] <- proplist_df[[paste0(cluster, '_transcript_count')]]
        } else {
          proplist_df[[paste0(cluster, '_transcript_count')]] <- rowSums(clusterdf)
          prop_calc_df[[paste0(cluster, '_transcript_count')]] <- proplist_df[[paste0(cluster, '_transcript_count')]]
        }

        # Cluster-specific gene counts
        a <- aggregate(proplist_df[[paste0(cluster, '_transcript_count')]], by = list(proplist_df$gene_id), FUN=sum)
        colnames(a) <- c(gene_id_colname, paste0(cluster, '_gene_count'))
        prop_calc_df <- merge(prop_calc_df, a, by = gene_id_colname, all.x = TRUE)
        proplist_df[[paste0(cluster, '_gene_count')]] <- prop_calc_df[[paste0(cluster, '_gene_count')]]

        # Calculate cluster-specific transcript proportion
        proplist_df[['proportion_in_cluster']] <- prop_calc_df[[paste0(cluster, '_transcript_count')]] / prop_calc_df[[paste0(cluster, '_gene_count')]]

        # Find difference in cluster-specific proportion and average proportion
        proplist_df[['difference_in_proportions']] <- proplist_df[['proportion_in_cluster']] - proplist_df[['all_cells_prop']]

        # Calculate proportion in all OTHER cells
        proplist_df[['transcript_count_outside_cluster']] <- proplist_df$total_transcript_count - proplist_df[[paste0(cluster, '_transcript_count')]]
        proplist_df[['gene_count_outside_cluster']] <- proplist_df$total_gene_count - proplist_df[[paste0(cluster, '_gene_count')]]
        proplist_df[['proportion_outside_cluster']] <- proplist_df$transcript_count_outside_cluster / proplist_df$gene_count_outside_cluster

        # Find difference in cluster-specific proportion and proportion in other cells
        proplist_df[['difference_in_vs_out_of_cluster']] <- proplist_df$proportion_in_cluster - proplist_df$proportion_outside_cluster

        # Add to proplist to end the loop
        proplist[[cluster]] <- proplist_df
    }

    return(proplist)
}
