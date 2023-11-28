#' Runs all functions in the pipeline
#'
#' Runs all functions needed to analyse differential isoform usage, and returns outputs in a list.
#' In order: (1) `filter_counts_table()`, (2) `generate_prop_tables()`, (3) `get_permutation_pvals()`, and (4) `make_plots`.
#'
#' @param transcript_counts_table A dataframe with a column containing transcript IDs, a column with gene IDs for
#' each transcript, and counts for all transcripts in all cells, with cell IDs as subsequent column names.
#' @param cell_labels_table A dataframe with two columns: one listing all cell IDs, and the other listing the group
#' that the cell ID belongs to.
#' @param transcript_id_colname A string corresponding to the column name in transcript_counts_table where the transcript IDs are stored.
#' @param gene_id_colname A string corresponding to the column name in transcript_counts_table where the gene IDs are stored.
#' @param cell_labels_colname A string corresponding to the column name in cell_labels_table where the group information is stored.
#' @param filter_genes_with_one_transcript A boolean defining whether genes with one associated transcript should be removed. Defaults to TRUE.
#' @param transcript_count_threshold An integer. If the total counts for a transcript is below this number, the transcript will
#' be filtered out of the table. Defaults to 0, which disables this filtering step. Note that this happens before filtering genes with
#' one transcript, so any genes that contain one transcript after count filtering will also be filtered.
#' @param gene_count_threshold An integer. If the total counts for all transcripts in a gene is below this number, the transcripts will be filtered
#' out of the table. Defaults to 0, which disables this filtering step.
#' @param permutations An integer corresponding to the number of permutations that should be run. Default 10000. Note that the first p-value calculation
#' is performed on the original set of data, and label shuffling occurs from permutation 2 and onwards.
#' @param cores An integer representing the number of cores that the function should attempt to use. If left undefined, set to 0, or set
#' to more cores than are available, defaults to all available cores.
#' @param get_detailed_pvalue_output A boolean, defaults to `FALSE`. If `TRUE`, keeps and returns extra p-value information needed to look at
#' permuted p-value distribution for any given transcripts. Note that this may take up very large amounts (>60GB) of RAM, and is therefore
#' only recommended on cluster computing systems with large amounts of available RAM.
#' @param do_gene_level_comparisons A boolean, defaults to `FALSE`. If `TRUE`, also runs a permutation analysis on gene level differences.
#' A significant value indicates that there is a difference in transcript proportions visible at the gene level, but does not specify
#' which transcripts show the change within the gene.
#' @param save_plots If set to a filepath to a directory, resulting DIU plots will be saved to the directory. Default `NA`.
#' @return A list containing the filtered counts table object, proportions tables, a permutation pvalue object, and the resulting DIU plots.
#' @examples
#' counts_table <- data.frame('transcript_id' = c(1, 2, 3), 'gene_id' = c(1, 1, 2), 'Cell_1' = c(0, 1, 10), 'Cell_2' = c(10, 2, 5), 'Cell_3' = c(2, 5, 1))
#' labels_table <- data.frame('grouping' = c('Cluster_1', 'Cluster_2', 'Cluster 1'), 'Cells' = c('Cell_1', 'Cell_2', 'Cell_3'))
#' results <- run_everything(counts_table, labels_table, 'transcript_id', 'gene_id', 'grouping')
#' @export


run_everything <- function(transcript_counts_table, cell_labels_table,
                           transcript_id_colname, gene_id_colname,
                           cell_labels_colname, filter_genes_with_one_transcript = TRUE,
                           transcript_count_threshold = 0, gene_count_threshold = 0,
                           permutations = 10000, cores = 0, get_detailed_pvalue_output = FALSE,
                           do_gene_level_comparisons = FALSE, save_plots = FALSE) {

    filtered_counts_table <- filter_counts_table(transcript_counts_table, transcript_id_colname, gene_id_colname,
                                                 filter_genes_with_one_transcript, gene_count_threshold)

    proplist <- generate_prop_tables(filtered_counts_table$counts_table, cell_labels_table, transcript_id_colname, gene_id_colname, cell_labels_colname)

    permutation_pval_object <- get_permutation_pvals(filtered_counts_table$counts_table, cell_labels_table, transcript_id_colname, gene_id_colname,
                                                     cell_labels_colname, permutations, cores, do_gene_level_comparisons, get_detailed_pvalue_output)

    plots <- make_plots(proplist, permutation_pval_object, transcript_id_colname, save_plots = save_plots)

    to_return <- list('filtered_counts_table' = filtered_counts_table, 'proportions_tables' = proplist,
                      'permutation_pvalue_object' = permutation_pval_object, 'plots' = plots$plots, 'tables_for_plots' = plots$tables)
    return(to_return)
}
