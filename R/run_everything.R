#' Runs all functions in the pipeline
#'
#' Runs all functions needed to analyse differential isoform usage, and returns outputs in a list.
#' In order: (1) `filter_counts_table()` and (2) `get_permutation_pvals()`.
#' For ease of use, available options are reduced when running the functions here. For more detailed
#' usage and for more low-level control, run each function separately. In additon, the permutation step
#' runs all groups against the rest of the dataset. For comparing two specific groups, or to compare a
#' single group against the rest of the dataset, run `get_permutation_pvals()` separately.
#'
#' @param transcript_counts_table A dataframe with a column containing transcript IDs, a column with gene IDs for
#' each transcript, and all other columns containing counts for all transcripts in all cells. Cell IDs should be
#' the column names for all columns containing these counts.
#' @param cell_labels_table A dataframe with two columns: one with all cell IDs, and the other with the group
#' that each cell ID belongs to.
#' @param transcript_id_colname A string corresponding to the column name in transcript_counts_table where the 
#' transcript IDs are stored.
#' @param gene_id_colname A string corresponding to the column name in transcript_counts_table where the gene IDs 
#' are stored.
#' @param cell_labels_colname A string corresponding to the column name in cell_labels_table where the group 
#' information is stored.
#' @param cell_group_to_analyse A string corresponding to the cell group to compare against the rest of the dataset. 
#' This should be a cell group in cell_labels_table.
#' @param remove_outlier_cells A boolean. If `TRUE`, calls scuttle's `is.outlier()` function using nmads = 5, 
#' log = T, and type = higher to filter out cells with extremely high counts. Defaults to `TRUE`.
#' @param output_folder A filepath to an output folder where data will be saved. If the folder doesn't exist, Isopod 
#' will attempt to create it.
#' @param gene_count_threshold An integer. If the total counts for all transcripts in a gene is below this number, 
#' the transcripts will be filtered out of the table. Defaults to 20.
#' @param collapse_isoforms_with_counts_below An integer. If set to an integer higher than 1, isoforms within a gene 
#' that have n total counts or less will be collapsed together at a per-cell basis. This collapsing step is useful 
#' for reducing significant isoforms with very low counts, which may be biologically insignificant. This also may 
#' help reduce the number of p-value calculations in the permutation step, reducing computational load. Defaults to 6. 
#' Setting this to 0 disables the feature.
#' @param permutations An integer corresponding to the number of permutations that should be run. Default 10000. 
#' Note that the first p-value calculation is performed on the original set of data, and label shuffling occurs 
#' from permutation 2 and onwards.
#' @param cores An integer representing the number of cores that the function should attempt to use. If left undefined, 
#' set to 0, or set to more cores than are available, defaults to all available cores.
#' @param run_on_all_groups A boolean. If TRUE, ignores the cell_group_to_analyse argument and runs every group in a
#' 1 vs rest fashion. This can take a long time for large datasets, so should normally be FALSE. Defaults to FALSE.
#' @param do_gene_level_comparisons A boolean. If TRUE, also runs a permutation analysis on gene level differences.
#' A significant value indicates that there is a difference in transcript proportions visible at the gene level, but 
#' does not specify which transcripts show the change within the gene. Enabling this also allows for other automatic 
#' calculations such as odds-ratio calculations. Defaults to TRUE.
#' @param generate_UMAPs A boolean. Determines whether a UMAP should be calculated and plotted from the 
#' provided data. UMAP calculation can be somewhat computationally expensive, and therefore may be disabled if needed.
#' Defaults to `TRUE`.
#' @param verbose A boolean. Outputs verbose logging in the permutation function. Useful for cluster environments, where
#' the normal output might not be correctly captured in stdout.
#' @return A list containing the output of `get_permutation_pvals()` and `make_plots()`.
#' @examples
#' counts_table <- data.frame('transcript_id' = c(1, 2, 3), 'gene_id' = c(1, 1, 2), 
#'                            'Cell_1' = c(0, 1, 10), 'Cell_2' = c(10, 2, 5), 
#'                            'Cell_3' = c(2, 5, 1))
#' labels_table <- data.frame('grouping' = c('Cluster_1', 'Cluster_2', 'Cluster 1'), 
#'                            'Cells' = c('Cell_1', 'Cell_2', 'Cell_3'))
#' results <- run_everything(counts_table, labels_table, 'transcript_id', 'gene_id', 
#'                           'grouping', output_folder = 'output_folder')
#' @export


run_everything <- function(transcript_counts_table, cell_labels_table,
                           transcript_id_colname, gene_id_colname,
                           cell_labels_colname, output_folder, cell_group_to_analyse = NA,
                           remove_outlier_cells = T, gene_count_threshold = 20, 
                           collapse_isoforms_with_counts_below = 6, permutations = 10000, 
                           cores = 0, run_on_all_groups = FALSE,
                           do_gene_level_comparisons = TRUE, generate_UMAPs = TRUE,
                           verbose = FALSE) {

    
    # Check that cell_group_to_analyse or run_on_all_groups is set
    if (is.na(cell_group_to_analyse) & !run_on_all_groups) {
        stop('Either cell_group_to_analyse or run_on_all_groups must be set. Set ONE of these variables and try running the function again.')
    }
    
    # Quick check/trim of output folder string
    lastchar <- function(x){
        substr(x, nchar(x), nchar(x))
    }

    if ((lastchar(output_folder) == '/') | (lastchar(output_folder) == '\\')) {
        output_folder <- substr(output_folder, 1, nchar(output_folder) - 1)
    }

    # Check if output folder exists. 
    # If so, throw a warning and wait for user confirmation.
    if (dir.exists(output_folder)) {
        warning(paste0(output_folder, ' already exists, and will be overwritten.\n',
                       'This warning can be turned off by setting "disable_overwrite_warning = TRUE" when running the function.\n'))
    }
    
    dir.create(output_folder, showWarnings = FALSE)

    cat('------------------------------------------\n')
    cat('              Filtering dataset\n')
    cat('------------------------------------------\n\n')
    cat('1. Filtering...\n')
    filtered_counts_table <- filter_counts_table(transcript_counts_table, transcript_id_colname, gene_id_colname,
                                                 gene_count_threshold, autofiltering = F,
                                                 remove_outlier_cells = remove_outlier_cells,
                                                 cell_labels_table = cell_labels_table, cell_labels_colname = cell_labels_colname,
                                                 collapse_isoforms_with_counts_below = collapse_isoforms_with_counts_below)

    cat('\n2. Saving filtered counts table to', output_folder, '\n\n')
    filtered_iso_list <- filtered_counts_table$list_of_collapsed_isoforms
    filtered_counts_table <- filtered_counts_table$counts_table
    saveRDS(filtered_iso_list, file.path(output_folder, 'list_of_collapsed_isoforms.rds'))
    saveRDS(filtered_counts_table, file.path(output_folder, 'filtered_counts_table.rds'))

    cat('------------------------------------------\n')
    cat('     Starting permutation analysis...\n')
    cat('------------------------------------------\n\n')
    cat('3. Running the permutation analysis using', permutations, 'permutations...\n')
    if (!run_on_all_groups) {
        cat('\tRunning', cell_group_to_analyse, 'against all other groups...\n')
        permutation_pval_object <- get_permutation_pvals(filtered_counts_table,
                                                         cell_labels_table,
                                                         transcript_id_colname,
                                                         gene_id_colname,
                                                         cell_labels_colname,
                                                         analysis_group_1 = cell_group_to_analyse,
                                                         permutations = permutations,
                                                         cores = cores,
                                                         do_gene_level_comparisons = do_gene_level_comparisons,
                                                         verbose = verbose)
    } else if (run_on_all_groups) {
        cat('\tRunning all groups against the dataset...\n')
        permutation_pval_object <- get_permutation_pvals(filtered_counts_table,
                                                         cell_labels_table,
                                                         transcript_id_colname,
                                                         gene_id_colname,
                                                         cell_labels_colname,
                                                         run_on_all_groups = run_on_all_groups,
                                                         permutations = permutations,
                                                         cores = cores,
                                                         do_gene_level_comparisons = do_gene_level_comparisons,
                                                         verbose = verbose)
    }
    

    cat('\n4. Saving permutation results to', output_folder, '\n\n')
    saveRDS(permutation_pval_object, file.path(output_folder, 'permutation_results.rds'))
    
    cat('------------------------------------------\n')
    cat('       Analysing and creating plots...\n')
    cat('------------------------------------------\n\n')
    cat('5. Finding proportions and creating plots...\n')
    plots <- make_plots(filtered_counts_table, 
                        cell_labels_table,
                        transcript_id_colname,
                        gene_id_colname,
                        cell_labels_colname, 
                        permutation_pval_object, 
                        gene_of_interest = NA,
                        generate_UMAPs = TRUE, 
                        save_to_folder = output_folder)

    cat('\nFinished! Check', output_folder, 'for results.\n')
    return(list(permutation_results = permutation_pval_object,
                plots = plots))
}
