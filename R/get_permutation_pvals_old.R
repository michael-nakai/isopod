#' Perform a permutation analysis
#'
#' Performs a permutation analysis using repeated chi-squared measures of significance, while shuffling cell grouping designations
#' between each permutation. Should be run a minimum of 10,000 times to guarantee consistency in results.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @param transcript_counts_table A dataframe with a column containing transcript IDs, a column with gene IDs for
#' each transcript, and counts for all transcripts in all cells, with cell IDs as subsequent column names.
#' @param cell_labels_table A dataframe with two columns: one listing all cell IDs, and the other listing the group
#' that the cell ID belongs to.
#' @param transcript_id_colname A string corresponding to the column name in `transcript_counts_table` where the transcript IDs are stored.
#' @param gene_id_colname A string corresponding to the column name in `transcript_counts_table` where the gene IDs are stored.
#' @param cell_labels_colname A string corresponding to the column name in `cell_labels_table` where the group information is stored.
#' @param permutations An integer corresponding to the number of permutations that should be run. Default 10000. Note that the first p-value calculation
#' is performed on the original set of data, and label shuffling occurs from permutation 2 and onwards. If only 1 permmutation is chosen, the
#' output p-values will be the chi-squared p-values.
#' @param cores An integer representing the number of cores that the function should attempt to use. If left undefined, set to 0, or set
#' to more cores than are available, defaults to all available cores. A good starting point is to allocate at least N cores, where
#' N = (number of cell groupings * 2). The function will run at varying degrees of reduced speed if this threshold isn't met.
#' @param do_gene_level_comparisons If set to `TRUE`, also runs a permutation analysis on gene level differences. A significant value indicates that
#' there is a difference in transcript proportions visible at the gene level, but does not specify which transcripts show the change within the gene.
#' @param return_detailed_pvalue_tables If set to `TRUE`, two additonal lists of dataframes will be returned containing additional information on
#' p-values generated over the permutations. If you'd like to look into p-value distributions, this must be set to `TRUE`. Note that if this is
#' set to `TRUE`, the resulting object generated will be very large (normally ~100GB), and therefore this option is only recommended for detailed
#' analysis that is usually unnecessary.
#' @param checkpoint_every_n_loops Defaults to 0. If set to an integer, a file will be saved in the current working directory containing the pvalue_storage
#' variable every n permutations.
#' @param checkpoint_file_location Defaults to `NA`. A filepath to a folder that will be saved if `checkpoint_every_n_loops` is set to an integer.
#' @return A `list` of five objects: the first is a list dataframes containing p-values for each permutation, for each group.
#' The second is a list of dataframes that describe the values used in pvalue calculation for every transcript, for each group,
#' only on the initial run (while tables were not permuted). The third and fourth are also `NA` if `return_detailed_pvalue_tables`
#' is set to `FALSE`, otherwise they contain additional information on the p-values over permutations. The fifth is also `NA` if
#' `do_gene_level_comparisons` is set to `FALSE`, otherwise contains a dataframe with permutation p-values for gene level transcript
#' proportion differences (see the description for do_gene_level_commparisons).
#' @examples
#' counts_table <- data.frame('transcript_id' = c(1, 2, 3), 'gene_id' = c(1, 1, 2), 'Cell_1' = c(0, 1, 10), 'Cell_2' = c(10, 2, 5), 'Cell_3' = c(2, 5, 1))
#' labels_table <- data.frame('grouping' = c('Cluster_1', 'Cluster_2', 'Cluster 1'), 'Cells' = c('Cell_1', 'Cell_2', 'Cell_3'))
#' permutation_results <- get_permutation_pvals(counts_table, labels_table, 'transcript_id', 'gene_id', 'grouping', 10000, 0, TRUE, FALSE, 6000, '/user/example/folder/')
#' @export

get_permutation_pvals_old <- function(transcript_counts_table, cell_labels_table,
                                  transcript_id_colname, gene_id_colname,
                                  cell_labels_colname, permutations = 10000,
                                  cores = 0, do_gene_level_comparisons = FALSE,
                                  return_detailed_pvalue_tables = FALSE,
                                  checkpoint_every_n_loops = 0, checkpoint_file_location = NA) {

    # Reduced chisq functions for isoforms and genes
    chisq.slim.test <- function(Nt, Rt, Ng, Rg)
    {
        x <- matrix(c(c(Nt, Ng), c(Rt, Rg)), ncol = 2)
        n <- sum(x)

        nr <- as.integer(nrow(x))
        nc <- as.integer(ncol(x))
        sr <- rowSums(x)
        sc <- colSums(x)
        E <- outer(sr, sc, "*") / n

        YATES <- min(0.5, abs(x-E))

        STATISTIC <- sum((abs(x - E) - YATES)^2 / E)
        PARAMETER <- (nr - 1L) * (nc - 1L)
        PVAL <- stats::pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)

        return(PVAL)
    }

    chisq.slim.gene.test <- function(count_vector_in, count_vector_out)
    {
        x <- matrix(c(count_vector_in, count_vector_out), ncol = 2)
        n <- sum(x)

        nr <- as.integer(nrow(x))
        nc <- as.integer(ncol(x))
        sr <- rowSums(x)
        sc <- colSums(x)
        E <- outer(sr, sc, "*") / n

        YATES <- min(0.5, abs(x-E))

        STATISTIC <- sum((abs(x - E) - YATES)^2 / E)
        PARAMETER <- (nr - 1L) * (nc - 1L)
        PVAL <- stats::pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)

        return(PVAL)
    }

    # Main body
    filtered_tcount <- transcript_counts_table
    designations <- cell_labels_table
    initial_vals <- NA
    newconst <- 0.00000001  # Add a very low number to avoid chi.sq errors due to a zero in the 2x2 table

    meta <- filtered_tcount %>% dplyr::select({{transcript_id_colname}}, {{gene_id_colname}})
    filtered_tcount[[transcript_id_colname]] <- NULL
    filtered_tcount[[gene_id_colname]] <- NULL

    # Get total counts per transcript, add it to meta
    meta$counts_per_transcript_in_all_cells <- rowSums(filtered_tcount)

    # Split cells into cluster-specific counts tables within a list
    clusters <- unique(designations[[cell_labels_colname]])[order(unique(designations[[cell_labels_colname]]))]

    # Figure out which column in designations has the grouping labels, and which has the cell labels
    cell_barcodes_col <- ifelse(colnames(designations)[1] == cell_labels_colname, 2, 1)

    # For resetting purposes
    filtered_tcount_base <- filtered_tcount
    meta_base <- meta

    # Initialize pvalue_storage_new so you don't have to do it in the loop
    a <- as.data.frame(cbind(meta[[transcript_id_colname]], meta[[gene_id_colname]], rep(NA, nrow(meta))))
    colnames(a) <- c('transcript_id', 'gene_id', 'pval1')
    sublist <- rep(list(a), length(clusters))
    names(sublist) <- clusters
    pvalue_storage_new <- rep(list(sublist), permutations)
    names(pvalue_storage_new) <- as.character(1:permutations)

    # If do_gene_level_comparisons is set, initialize pvalue_storage_gene here
    if (do_gene_level_comparisons) {
        a <- data.frame('gene_id' = unique(meta[[gene_id_colname]]), 'pval1' = rep(NA, length(unique(meta[[gene_id_colname]]))))
        colnames(a) <- c('gene_id', 'pval1')
        sublist <- rep(list(a), length(clusters))
        names(sublist) <- clusters
        pvalue_storage_gene <- rep(list(sublist), permutations)
        names(pvalue_storage_gene) <- as.character(1:permutations)
    }

    #### LOOP START

    # For parallelizing foreach later on, make the socket first
    if (cores == 0) {
        cores = parallelly::availableCores()
    }

    # Check that there's enough cores to run everything efficiently.
    # If there's equal or less cores than there are clusters, we need to make sure things can actually run, speed be damned.
    if (cores <= length(clusters)) {
        mcmapply_reserved_cores_total <- 1
        cl <- parallel::makeCluster(cores - mcmapply_reserved_cores_total)
        doParallel::registerDoParallel(cl)
        mcmapply_reserved_cores <- floor(mcmapply_reserved_cores_total / length(clusters))
        if (mcmapply_reserved_cores == 0) {
            mcmapply_reserved_cores = 1
        }
    } else {  # If there's enough cores, allocate them correctly.
        mcmapply_reserved_cores_total <- cores - length(clusters)
        cl <- parallel::makeCluster(cores - mcmapply_reserved_cores_total)
        doParallel::registerDoParallel(cl)
        mcmapply_reserved_cores <- floor(mcmapply_reserved_cores_total / length(clusters))
        if (mcmapply_reserved_cores == 0) {
            mcmapply_reserved_cores = 1
        }
    }

    # Main permutation loop
    for (loopnum in 1:permutations) {

        message(paste('Loop number', loopnum))

        # Estimating time left (time1 is first set during the first loop in the if (loopnum == 1) block further down)
        if ((loopnum %% 100 == 0) & (loopnum != 1)) {
            time2 <- Sys.time()
            estimated_time <- round((as.numeric(difftime(time2, time1, units = 'mins')) / 100 * (permutations - loopnum)), 2)
            message(paste0('Estimated time to completion: ', estimated_time, ' minutes'))
            time1 <- time2
        }

        # Reset filtered_tcount and meta (just in case, shouldn't actually change anything)
        filtered_tcount <- filtered_tcount_base
        meta <- meta_base

        # If we're on cycle 2+, then scramble the grouping column in designations
        if (loopnum != 1) {
            designations[[cell_labels_colname]] <- sample(designations[[cell_labels_colname]])
        }

        cluster_counts_list <- list()
        for (cluster in clusters) {
            cluster_counts_list[[cluster]] <- filtered_tcount[designations[which(designations[[cell_labels_colname]] == cluster), ][[colnames(designations)[cell_barcodes_col]]]]
        }

        # rowSums per cluster in list are added to a column named 'counts_per_transcript_in_cluster' per table in cluster_counts_list
        # This will always be the last column
        cluster_counts_list <- lapply(cluster_counts_list, function(x) {
            x$counts_per_transcript_in_cluster <- rowSums(x)
            return(x)
        })

        # Retrieve the rowSums and add to meta, with cluster-specific colunm names
        for (groupname in names(cluster_counts_list)) {
            column_name <- paste0('counts_per_transcript_in_', groupname)
            meta[[column_name]] <- cluster_counts_list[[groupname]][ncol(cluster_counts_list[[groupname]])]
        }
        #saveRDS(meta, file = "/oshlack_lab/michael.nakai/projects/testing_isopod/datasets/pacbio_pbmc/pbmc_10k_rds_saves/meta_debug_1.rds")

        # Get counts per gene in all cells, then by cluster
        # Turn meta into a data.table here, since it speeds up merging quite a bit
        meta <- data.table::as.data.table(meta, key = gene_id_colname)
        for (column_name in colnames(meta)[3:ncol(meta)]) {
            sums <- aggregate(meta[[column_name]], list(temp_gene_id_here_col = meta[[gene_id_colname]]), sum)
            clus <- paste0('counts_per_gene_in_cluster_', str_split(column_name, '_')[[1]][[5]])
            if (clus == 'counts_per_gene_in_cluster_all') {
                clus <- 'counts_per_gene_in_all_cells'
            }
            colnames(sums) <- c(gene_id_colname, clus)
            sums <- data.table::as.data.table(sums)
            meta <- data.table::merge.data.table(meta, sums, by = gene_id_colname, all.x = TRUE)
        }

        ### Calculate proportions for all and per cluster, add them into meta ###
        # Find column index for 'counts_per_transcript_in_all_cells', until (but not including) 'counts_per_gene_in_all_cells'
        index_with_all_cells <- grep('counts_per_transcript_in_all_cells', colnames(meta))
        index_with_gene_counts <- grep('counts_per_gene_in_all_cells', colnames(meta))

        # Make meta back into a data.frame, since we only need data.table for merging
        meta <- as.data.frame(meta)

        #saveRDS(meta, file = "/oshlack_lab/michael.nakai/projects/testing_isopod/datasets/pacbio_pbmc/pbmc_10k_rds_saves/meta_debug_2.rds")

        # Loop over column names starting from the counts_per_gene_in_cluster directly after all_cells
        colnames_to_look_at <- colnames(meta)[(index_with_all_cells+1):(index_with_gene_counts-1)]
        for (i in 1:length(colnames_to_look_at)) {

            # Calculate counts_per_transcript_in_all_cells / meta[[column_name]], then add that to a new column called paste0(cluster_name, '_transcript_within_gene_proportion')
            column_name <- colnames_to_look_at[i]
            current_cluster <- str_split(column_name, '_')[[1]][5]
            gene_counts_column_name <- paste0('counts_per_gene_in_cluster_', current_cluster)

            newColNameClust <- paste0('proportion_for_', current_cluster)
            newColNameRest <- paste0('rest_proportion_for_', current_cluster) # This is the colname for the rest of the data, since we're comparing between cluster 1, then all other clusters
            newColNameRestCountsTranscript <- paste0('counts_per_transcript_rest_in_cluster_', current_cluster)
            newColNameRestCountsGene <- paste0('counts_for_gene_rest_in_cluster_', current_cluster)

            meta[[newColNameRestCountsTranscript]] <- meta$counts_per_transcript_in_all_cells - meta[[column_name]]
            meta[[newColNameRestCountsGene]] <- meta$counts_per_gene_in_all_cells - meta[[gene_counts_column_name]]

            meta[[newColNameClust]] <- meta[[column_name]] / meta[[gene_counts_column_name]]
            meta[[newColNameRest]] <- meta[[newColNameRestCountsTranscript]] / meta[[newColNameRestCountsGene]]

        }

        # If the current loopnum is 1, make the second list of DFs here (for plotting steps)
        # It's wordy, but only needs to run once, so shouldn't impact speed
        if (loopnum == 1) {
            time1 <- Sys.time()
            initial_vals <- vector('list', length(clusters))
            for (clustername in clusters) {
                ntstr <- paste0('counts_per_transcript_in_', clustername)
                rtstr <- paste0('counts_per_transcript_rest_in_cluster_', clustername)
                ngstr <- paste0('counts_per_gene_in_cluster_', clustername)
                rgstr <- paste0('counts_for_gene_rest_in_cluster_', clustername)

                initial_vals_df <- meta %>% dplyr::select({{transcript_id_colname}}, {{gene_id_colname}})
                initial_vals_df[['nt']] <- meta[[ntstr]]
                initial_vals_df[['rt']] <- meta[[rtstr]]
                initial_vals_df[['ng']] <- meta[[ngstr]] - meta[[ntstr]]
                initial_vals_df[['rg']] <- meta[[rgstr]] - meta[[rtstr]]

                initial_vals[[clustername]] <- initial_vals_df
            }
        }

        # Do the gene-level pval calculations here
        if (do_gene_level_comparisons) {
            tempresultsgene <- foreach::foreach (iter = 1:length(clusters), .combine = 'cbind', .multicombine = T, .packages = 'parallel') %dopar% {
                cluster_name <- clusters[iter]
                ntstr <- paste0('counts_per_transcript_in_', cluster_name)
                rtstr <- paste0('counts_per_transcript_rest_in_cluster_', cluster_name)

                vec_in <- lapply(split(meta[[ntstr]], meta[[gene_id_colname]]), FUN = function(x) c(x) + newconst)
                vec_out <- lapply(split(meta[[rtstr]], meta[[gene_id_colname]]), FUN = function(x) c(x) + newconst)

                pval_vector <-  parallel::mcmapply(chisq.slim.gene.test,
                                                   vec_in,
                                                   vec_out,
                                                   mc.cores = mcmapply_reserved_cores)

                pval_vector
            }
            i <- 1
            for (clust_name in names(pvalue_storage_gene[[as.character(loopnum)]])) {
                pvalue_storage_gene[[as.character(loopnum)]][[clust_name]]$pval1 <- tempresultsgene[, i]
                i <- i + 1
            }
        }

        ### Do the chi-square test ###
        #
        # CONTINGENCY TABLE
        # For each transcript, make a 2x2 contingency table that looks like:
        #
        # Nt      Rt
        # Ng-Nt   Rg-Rt
        #
        # Where:
        # Nt = Transcript count in cluster (meta$counts_per_transcript_in_cluster_x)
        # Rt = Transcript count all other cells (meta$counts_per_transcript_rest_in_cluster_x)
        # Ng = All counts of all transcripts for gene in cluster (meta$counts_for_gene_in_cluster_x)
        # Rg = All counts of all transcripts for gene in rest (meta$counts_for_gene_rest_in_cluster_x)

        # If this is the first loop, we need to find the indices to skip, so we need to store the transcript_id column to relate back later
        if (loopnum == 1) {
            transcript_ids <- meta[[transcript_id_colname]]
        }

        # Temporarily cut off transcript_id, and gene_id for mapply to work
        meta[[transcript_id_colname]] <- NULL
        meta[[gene_id_colname]] <- NULL

        # Loop over all clusters (in the same order as colnames_to_look_at)
        tempresults <- foreach::foreach (iter = 1:length(clusters), .combine = 'cbind', .multicombine = T, .packages = 'parallel') %dopar% {
            cluster_name <- clusters[iter]
            ntstr <- paste0('counts_per_transcript_in_', cluster_name)
            rtstr <- paste0('counts_per_transcript_rest_in_cluster_', cluster_name)
            ngstr <- paste0('counts_per_gene_in_cluster_', cluster_name)
            rgstr <- paste0('counts_for_gene_rest_in_cluster_', cluster_name)

            pval_vector <-  parallel::mcmapply(chisq.slim.test,
                                               meta[[ntstr]] + newconst,
                                               meta[[rtstr]] + newconst,
                                               meta[[ngstr]] - meta[[ntstr]] + newconst,
                                               meta[[rgstr]] - meta[[rtstr]] + newconst,
                                               mc.cores = mcmapply_reserved_cores)

            pval_vector

        }
        i <- 1
        for (clust_name in names(pvalue_storage_new[[as.character(loopnum)]])) {
            pvalue_storage_new[[as.character(loopnum)]][[clust_name]]$pval1 <- tempresults[, i]
            i <- i + 1
        }

        # TODO
        # If this is the first loop, find the indices of the isoforms that have a p < cutoff, and add them to a comparison vector here
        if (loopnum == 1) {
            comparison_vec <- c()
        }

        # Save file if checkpoint is reached
        if ((checkpoint_every_n_loops > 0) & (!(is.na(checkpoint_file_location)))) {
            if (loopnum %% checkpoint_every_n_loops == 0) {
                message('Checkpoint reached, saving transcript pvalues...')
                saveRDS(pvalue_storage_new, paste0(checkpoint_file_location, 'pvalue_storage_new.rds'))
                message(paste0('Saved to: ', checkpoint_file_location, 'pvalue_storage_new.rds'))
                if (do_gene_level_comparisons) {
                    message('Checkpoint reached, saving gene pvalues...')
                    saveRDS(pvalue_storage_gene, paste0(checkpoint_file_location, 'pvalue_storage_gene.rds'))
                    message(paste0('Saved to: ', checkpoint_file_location, 'pvalue_storage_gene.rds'))
                }
            }
        }
    }

    parallel::stopCluster(cl)

    # Start calculating permutation p-values here
    pval_table <- pvalue_storage_new$'1'
    for (loopname in names(pvalue_storage_new)) {
        if (loopname != '1') {
            for (clustername in clusters) {
                pval_table[[clustername]] <- cbind(pval_table[[clustername]], as.data.frame(pvalue_storage_new[[loopname]][[clustername]][['pval1']]))
            }
        }
    }

    for (clus in names(pval_table)) {
        colnames(pval_table[[clus]]) <- c(c(transcript_id_colname, gene_id_colname), c(as.character(1:length(names(pvalue_storage_new)))))
    }

    # Count num of times p-val is equal or less than initial p-val
    initial_pvals <- pvalue_storage_new$'1'  # This is a list of vectors
    for (clustername in names(initial_pvals)) {
        initial_pvals[[clustername]] <- initial_pvals[[clustername]][['pval1']]
    }

    # Initialize df to store
    pval_count <- data.frame(matrix(ncol = length(names(pvalue_storage_new$'1')) + 2, nrow = length(pvalue_storage_new$'1'[[clusters[1]]][[transcript_id_colname]])))
    colnames(pval_count) <- c(transcript_id_colname, gene_id_colname, names(pval_table))
    pval_count[[transcript_id_colname]] <- pvalue_storage_new$'1'[[clusters[1]]][[transcript_id_colname]]
    pval_count[[gene_id_colname]] <- pvalue_storage_new$'1'[[clusters[1]]][[gene_id_colname]]

    # Count here
    if (permutations > 1) {
        for (clustername in names(pval_table)) {
            countvec <- pval_table[[clustername]][['2']] <= initial_pvals[[clustername]]
            for (column in pval_table[[clustername]][, 5:ncol(pval_table[[clustername]])]) {
                countvec <- countvec + (column <= initial_pvals[[clustername]])
            }
            pval_count[[clustername]] <- countvec / permutations
        }
    } else {
        for (clustername in names(pval_table)) {
            pval_count[[clustername]] <- initial_pvals[[clustername]]
        }
    }

    # Do it all again for gene-level pvals if do_gene_level_comparisons is TRUE
    if (do_gene_level_comparisons) {
        pval_table_gene <- pvalue_storage_gene$'1'
        for (loopname in names(pvalue_storage_new)) {
            if (loopname != '1') {
                for (clustername in clusters) {
                    pval_table_gene[[clustername]] <- cbind(pval_table_gene[[clustername]], as.data.frame(pvalue_storage_gene[[loopname]][[clustername]][['pval1']]))
                }
            }
        }

        for (clus in names(pval_table_gene)) {
            colnames(pval_table_gene[[clus]]) <- c(c(gene_id_colname), c(as.character(1:length(names(pvalue_storage_gene)))))
        }

        # Count num of times p-val is equal or less than initial p-val
        initial_pvals_gene <- pvalue_storage_gene$'1'  # This is a list of vectors
        for (clustername in names(initial_pvals)) {
            initial_pvals_gene[[clustername]] <- initial_pvals_gene[[clustername]][['pval1']]
        }

        # Initialize df to store
        pval_count_gene <- data.frame(matrix(ncol = length(names(pvalue_storage_gene$'1')) + 1, nrow = length(pvalue_storage_gene$'1'[[clusters[1]]][[gene_id_colname]])))
        colnames(pval_count_gene) <- c(gene_id_colname, names(pval_table_gene))
        pval_count_gene[[gene_id_colname]] <- pvalue_storage_gene$'1'[[clusters[1]]][[gene_id_colname]]

        # Count here
        if (permutations > 1) {
            for (clustername in names(pval_table_gene)) {
                countvec <- pval_table_gene[[clustername]][['2']] <= initial_pvals_gene[[clustername]]
                for (column in pval_table_gene[[clustername]][, 3:ncol(pval_table_gene[[clustername]])]) {
                    countvec <- countvec + (column <= initial_pvals_gene[[clustername]])
                }
                pval_count_gene[[clustername]] <- countvec / permutations
            }
        } else {
            for (clustername in names(pval_table_gene)) {
                pval_count_gene[[clustername]] <- initial_pvals_gene[[clustername]]
            }
        }
    }

    if (!(return_detailed_pvalue_tables)) {
        pvalue_storage_new <- NA
        pval_table <- NA
    } else if (!(do_gene_level_comparisons)) {
        pval_table_gene <- NA
    }

    to_return <- list('p_values' = pval_count, 'initial_pvalue_calculation_values' = initial_vals,
                      'pvalue_storage_list' = pvalue_storage_new, 'pval_table' = pval_table,
                      'pvalues_gene' = pval_count_gene)
    return(to_return)
}
