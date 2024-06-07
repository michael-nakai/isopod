#' Filter a transcript counts table
#'
#' Filters a transcript counts table on multiple user-defined variables, including the ability to filter out genes with
#' one associated transcript and genes with counts below a threshold.
#'
#' @param transcript_counts_table A dataframe with a column containing transcript IDs, a column with gene IDs for
#' each transcript, and counts for all transcripts in all cells, with cell IDs as subsequent column names.
#' @param transcript_id_colname A string corresponding to the column name in transcript_counts_table where the transcript IDs are stored.
#' @param gene_id_colname A string corresponding to the column name in transcript_counts_table where the gene IDs are stored.
#' @param gene_count_threshold An integer. If the total counts for all transcripts in a gene is below this number, the transcripts will be filtered
#' out of the table. Defaults to 0, which disables this filtering step.
#' @param autofiltering A boolean. If `TRUE`, removes genes where at least N cells per group do not contain counts for the gene, where N is half the cell
#' count of the smallest group. Defaults to `TRUE`.
#' @param cell_labels_table A dataframe with two columns: one listing all cell IDs, and the other listing the group that the cell ID belongs to. Only is
#' required when autofiltering is `TRUE`. Defaults to `NA`.
#' @param cell_labels_colname A string corresponding to the column name in `cell_labels_table` where the group information is stored. Only is required
#' when `cell_labels_table` is provided. Defaults to NA.
#' @param collapse_isoforms_with_counts_below An integer. If set to an integer higher than 1, isoforms within a gene that have n total counts or less will
#' be collapsed together at a per-cell basis. This collapsing step is useful for reducing significant isoforms with very low counts, which may be biologically
#' insignificant. This also may help reduce the number of p-value calculations in the permutation step, reducing computational load. Defaults to 6. Setting this
#' to 0 disables the feature.
#' @param autofilter_threshold A float between 0 - 1, representing the fraction of the smallest cell group needed to pass the autofiltering step.
#' @return An object containing the dataframe representing the filtered transcript counts table, and a list of isoforms collapsed in each gene if
#' `collapse_isoforms_with_counts_below` is greater than 0. The list of collapsed isoforms is of the form gene_id --> vector of collapsed isoform IDs.
#' If the gene itself was filtered due to all isoforms being collapsed (i.e. if all isoforms in a gene didn't reach the threshold), the gene_id will return
#' a character string explaining that all isoforms associated with the gene was removed.


filter_counts_table <- function(transcript_counts_table,
                                transcript_id_colname, gene_id_colname,
                                gene_count_threshold = 0,
                                autofiltering = TRUE,
                                cell_labels_table = NA,
                                cell_labels_colname = NA,
                                collapse_isoforms_with_counts_below = 6,
                                autofilter_threshold = 0.5) {


    ### Argument sanity checking
    # Check that transcript_counts_table is a dataframe, matrix, or NULL object. If NULL, a more descriptive error will be raised later.
    stopifnot(typeof(transcript_counts_table) == 'NULL' | inherits(transcript_counts_table, "data.frame") | inherits(transcript_counts_table, 'matrix'))

    # Check for empty dataframe, and raise errors if true.
    if (is.null(dim(transcript_counts_table))) {
        stop("Your transcript_counts_table is an empty dataframe or NULL object.\nPlease fill out your counts table and rerun this function.")
    } else if (dim(transcript_counts_table)[1] == 0 | dim(transcript_counts_table)[2] == 0) {
        stop("Your transcript_counts_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    }

    # Check that transcript_id_colname and gene_id_colname are character strings.
    stopifnot(assertthat::is.string(transcript_id_colname) | assertthat::is.string(gene_id_colname))

    # Check that transcript_ and gene_id_colnames exist in colnames(transcript_counts_table)
    if (!(transcript_id_colname %in% colnames(transcript_counts_table))) {
        stop(paste0(transcript_id_colname,
                    " could not be found as a column name in your counts table.\nPlease check your column names and rerun this function."))
    } else if (!(gene_id_colname %in% colnames(transcript_counts_table))) {
        stop(paste0(gene_id_colname,
                    " could not be found as a column name in your counts table.\nPlease check your column names and rerun this function."))
    }

    #TODO: Make this work
    # Stop if autofiltering is TRUE, but a cell_labels_table isn't provided.
    #if ((autofiltering) & ((typeof(cell_labels_table) == 'logical') | !(is.na(cell_labels_colname)))) {
    #    stop('Autofiltering is set to TRUE, but the cell_labels_table or cell_labels_colname is not provided.\nPlease provide these arguments and rerun the function.')
    #}

    ### Function body

    # INTERNAL SETUP
    collapsed_list <- NA

    # Find gene counts
    write('Setting up data', stdout())
    transcript_counts_table <- as.data.frame(transcript_counts_table)

    transcript_counts_table <- transcript_counts_table[order(transcript_counts_table[[gene_id_colname]], transcript_counts_table[[transcript_id_colname]]), ]
    prop_calc_df <- transcript_counts_table %>% dplyr::select({{transcript_id_colname}}, {{gene_id_colname}})
    nolabel_transcript_counts_table <- transcript_counts_table[ , !(colnames(transcript_counts_table) %in% c(transcript_id_colname, gene_id_colname))]
    prop_calc_df[['total_transcript_count']] <- rowSums(nolabel_transcript_counts_table)
    a <- aggregate(prop_calc_df$total_transcript_count, by = list(prop_calc_df$gene_id), FUN=sum)
    colnames(a) <- c(gene_id_colname, 'total_gene_count')
    prop_calc_df <- merge(prop_calc_df, a, by = gene_id_colname, all.x = TRUE)

    # Make sure transcript_counts_table has the tid and gid cols as cols 1 and 2
    transcript_id_idx <- match(transcript_id_colname, colnames(transcript_counts_table))
    gene_id_idx <- match(gene_id_colname, colnames(transcript_counts_table))
    transcript_counts_table <- transcript_counts_table[, c(which(colnames(transcript_counts_table) == transcript_id_colname),
                                                           which(colnames(transcript_counts_table) == gene_id_colname),
                                                           which(!(colnames(transcript_counts_table) %in% c(transcript_id_colname, gene_id_colname))))]


    # Filter genes with one associated transcript if argument is set
    write('Filtering genes with one associated transcript', stdout())

    transcript_counts_table[['temp_gene_id_new']] <- transcript_counts_table[[gene_id_colname]]

    transcript_counts_table <- transcript_counts_table %>%
        dplyr::group_by(temp_gene_id_new) %>%
        dplyr::filter(n() > 1)

    transcript_counts_table[['temp_gene_id_new']] <- NULL
    transcript_counts_table <- transcript_counts_table[order(transcript_counts_table[[gene_id_colname]], transcript_counts_table[[transcript_id_colname]]), ]

    # Remove the genes from prop_calc_df too, just in case we're also filtering by gene_count_threshold
    prop_calc_df <- prop_calc_df[prop_calc_df[[transcript_id_colname]] %in% transcript_counts_table[[transcript_id_colname]], ]


    # Filter genes with counts below a threshold if argument is set
    if (gene_count_threshold > 0) {

        write(paste0('Filtering genes with total counts below ', gene_count_threshold), stdout())

        transcript_counts_table <- transcript_counts_table[prop_calc_df$total_gene_count >= gene_count_threshold, ]

        # Remove the genes from prop_calc_df too, just in case we're also filtering by gene_count_threshold
        prop_calc_df <- prop_calc_df[prop_calc_df[[transcript_id_colname]] %in% transcript_counts_table[[transcript_id_colname]], ]
    }


    # Run the autofiltering step here (heuristics filtering)
    if (autofiltering) {

        write('Running heuristics-based filtering', stdout())

        cell_id_colname <- colnames(cell_labels_table[, !(names(cell_labels_table) %in% cell_labels_colname), drop = FALSE])
        min_cells <- 20  # Set the magic number here for easy debugging later
        group_counts <- sort(table(cell_labels_table[[cell_labels_colname]]))

        # Check that there is a group with at least 20 cells included (so we'd check for 10 cells). If so, autofilter
        if (sum(group_counts >= min_cells) > 0) {
            cell_nums_to_check <- as.integer(group_counts[match(TRUE, (group_counts >= min_cells))] * autofilter_threshold)

            # This whole block is just trying to get the number of cells with a count per group, which is stored in per_group_counts_table
            per_group_counts_table <- transcript_counts_table %>% select({{transcript_id_colname}}, {{gene_id_colname}})
            for (group in unique(cell_labels_table[[cell_labels_colname]])) {
                cell_id_subset <- cell_labels_table[which(cell_labels_table[cell_labels_colname] == group), ][[cell_id_colname]]
                temp <- transcript_counts_table[, (colnames(transcript_counts_table) %in% c(gene_id_colname, cell_id_subset))]
                per_group_counts_table[[group]] <- rowSums(temp[, -1] > 0)
            }
            per_group_counts_table[[transcript_id_colname]] <- NULL
            colnames(per_group_counts_table) <- c('gene_id', colnames(per_group_counts_table)[2:ncol(per_group_counts_table)])
            per_group_counts_table <- aggregate(. ~ gene_id, data = per_group_counts_table, FUN = sum)

            # If a gene has more cells with a count than cell_nums_to_check, we keep the gene. Else, we don't.
            genes_to_keep <- c()
            for (i in 1:nrow(per_group_counts_table)) {
                #if (sum(per_group_counts_table[i, 2:ncol(per_group_counts_table)]) >= cell_nums_to_check) {  #Disregard groups, sum up between groups and compare to threshold
                #if (per_group_counts_table[i, 2] >= cell_nums_to_check) {  #Only looking at smallest group
                #if (sum(per_group_counts_table[i, 2:ncol(per_group_counts_table)] >= cell_nums_to_check) == ncol(per_group_counts_table) - 1) {  #Only keep if ALL groups have >= cell_nums_to_check
                if (sum(per_group_counts_table[i, 2:ncol(per_group_counts_table)] >= cell_nums_to_check) > 0) {  #At least one group above threshold
                    genes_to_keep <- c(genes_to_keep, per_group_counts_table[i, 1])
                }
            }

            # Filter the transcript_counts_table to only keep genes that we want
            transcript_counts_table <- transcript_counts_table[which(transcript_counts_table[[gene_id_colname]] %in% genes_to_keep), ]

            # Remove the genes from prop_calc_df too, just in case we're also filtering by gene_count_threshold
            prop_calc_df <- prop_calc_df[prop_calc_df[[transcript_id_colname]] %in% transcript_counts_table[[transcript_id_colname]], ]

        } else {
            stop(paste0('Autofiltering cannot be ran, as there are no groups containing at least ', min_cells, ' cells.'))
        }

    }


    # Run the collapse step here
    if (collapse_isoforms_with_counts_below > 0) {

        write(paste0('Collapsing isoforms with total counts below ', collapse_isoforms_with_counts_below, ' counts'), stdout())
        write('This may take a while (5-20 minutes) for large datasets', stdout())

        # Pre-make the collapsed_list here, of the form [[gene_id]] --> NULL
        collapsed_list = vector(mode = 'list', length = length(unique(transcript_counts_table[[gene_id_colname]])))
        names(collapsed_list) = unique(transcript_counts_table[[gene_id_colname]])

        # Fill out the list, so it becomes [[gene_id]] --> c('iso1', 'iso2') where iso1 and iso2 are isoforms that have less than the threshold
        for (i in 1:nrow(transcript_counts_table)) {
            tid <- transcript_counts_table[[transcript_id_colname]][i]
            gid <- transcript_counts_table[[gene_id_colname]][i]
            isocounts <- prop_calc_df$total_transcript_count[i]
            if (isocounts < collapse_isoforms_with_counts_below) {
                if (is.null(collapsed_list[[gid]])) {
                    collapsed_list[[gid]] <- c(tid)
                } else {
                    collapsed_list[[gid]] <- c(collapsed_list[[gid]], tid)
                }
            }
        }

        # For each element of the list, if there's more than one tid in the vector, then smush those isoforms together into 'genename_collapsed_isoforms'
        total_counts <- transcript_counts_table[, c(transcript_id_colname, gene_id_colname)]  # Make this to return to the user as some extra data
        total_counts$total_counts <- rowSums(transcript_counts_table[, !(colnames(transcript_counts_table) %in% c(transcript_id_colname, gene_id_colname))])
        table_vec <- table(transcript_counts_table[[gene_id_colname]])

        # Change the table to a dataframe (just in case it wasn't already)
        for_calcs <- transcript_counts_table[, -c(transcript_id_idx, gene_id_idx)]
        to_remove_indexes <- c()
        transcript_vec <- transcript_counts_table[[transcript_id_colname]]
        transcript_counts_table <- data.table::as.data.table(transcript_counts_table)
        replacement_vec <- 3:ncol(transcript_counts_table)

        a <- Sys.time()
        for (i in 1:length(names(collapsed_list))) {
            if (i %% 10 == 0) {
                write(paste0('Collapsing isoforms for gene ', i, ' / ', length(names(collapsed_list)),
                             ' (', round(i/length(names(collapsed_list))*100, 2), '%)'),
                      stdout())
            }

            gid <- names(collapsed_list)[i]

            # Get subset of transcripts that associate with the gene, if there's 2 or more isoforms to collapse and as long as we're not
            # collapsing ALL the isoforms in a gene (if we are, then look at the next block). Important to note here that if there's
            # only one isoform below the threshold for the gene, nothing happens (since collapsing wouldn't reduce rows in this case)
            if ((length(collapsed_list[[gid]]) >= 2) & (table_vec[gid] > length(collapsed_list[[gid]]))) {

                # Only work with vecs of rownums here to minimize CPU time
                rows_with_gene <- which(transcript_counts_table[[gene_id_colname]] == gid)
                rows_to_remove <- rows_with_gene[which(transcript_vec[rows_with_gene] %in% collapsed_list[[gid]])]

                # Get the colSums for the rows that we care about
                colsum_vec <- colSums(for_calcs[rows_to_remove, ])
                names(colsum_vec) <- NULL

                # Determine the rows to remove or replace
                row_to_replace <- rows_to_remove[1]
                rows_to_remove <- rows_to_remove[-1]

                # Update the DT
                set(transcript_counts_table, row_to_replace, transcript_id_colname, paste0(gid, '_collapsed_isoforms'))
                set(transcript_counts_table, row_to_replace, replacement_vec, as.list(colsum_vec))

                # Add the indexes to remove
                to_remove_indexes <- c(to_remove_indexes, rows_to_remove)

            # Else, if ALL isoforms in a gene are under the threshold, then remove the gene
            } else if (length(collapsed_list[[gid]]) == table_vec[gid]) {
                to_remove_indexes <- c(to_remove_indexes, which(transcript_counts_table[[gene_id_colname]] == gid))
            }
        }

        # Remove rows and finalize output
        if (!is.null(to_remove_indexes)) {
            transcript_counts_table <- transcript_counts_table[-to_remove_indexes, ]
            prop_calc_df <- prop_calc_df[-to_remove_indexes, ]
        }
        transcript_counts_table <- as.data.frame(transcript_counts_table)

        # Filter the collapsed_list list to only include entries with collapsed isoforms >= 2.
        for (x in names(collapsed_list)) {
            if ((length(collapsed_list[[x]]) < 2)) {
                collapsed_list[[x]] <- NULL
            } else if (table_vec[x] == length(collapsed_list[[x]])) {
                collapsed_list[[x]] <- paste0('Removed due to all isoforms having total counts < ', collapse_isoforms_with_counts_below)
            }
        }
    } else {
        collapsed_list <- NA
    }

    return(list('counts_table' = transcript_counts_table, 'list_of_collapsed_isoforms' = collapsed_list))
}
