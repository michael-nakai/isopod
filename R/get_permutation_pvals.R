#' Perform a permutation analysis
#'
#' Performs a permutation analysis using repeated chi-squared measures of significance, while shuffling cell grouping designations
#' between each permutation. Should be run a minimum of 10,000 times to guarantee consistency in results.
#'
#' @import data.table
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
#' @param return_detailed_pvalue_tables If set to `TRUE`, an additonal list of dataframes will be returned containing additional information on
#' p-values generated over the permutations. If you'd like to look into p-value distributions over the permutations, this must be set to `TRUE`.
#' Note that if this is set to `TRUE`, the resulting object generated will be very large (normally ~100GB), and therefore this option is only
#' recommended for further (usually unnecessary) inspection.
#' @param cutoff A `double` (decimal) indicating the initial p-value cutoff that isoforms (and genes if gene-level comparison is enabled) must pass to be included
#' in subsequent permutations. This cutoff applies for the first permutation only, and potentially increases the speed of each permutation by filtering down
#' data before calculations begin. With higher cutoffs, less data is filtered. To disable filtering, set the cutoff to 1. Defaults to `0.1`.
#' @return A `list` of four objects: the first is a list dataframes containing p-values for each permutation, for each group.
#' The second is a list of dataframes that describe the values used in pvalue calculation for every transcript, for each group,
#' only on the initial run (while tables were not permuted). The third is `NA` if `return_detailed_pvalue_tables`
#' is set to `FALSE`, otherwise it contains additional information on the p-values for each permutation. The fifth is also `NA` if
#' `do_gene_level_comparisons` is set to `FALSE`, but otherwise contains a dataframe with permutation p-values for gene level transcript
#' proportion differences (see the description for do_gene_level_commparisons).
#' @examples
#' counts_table <- data.frame('transcript_id' = c(1, 2, 3), 'gene_id' = c(1, 1, 2), 'Cell_1' = c(0, 1, 10), 'Cell_2' = c(10, 2, 5), 'Cell_3' = c(2, 5, 1))
#' labels_table <- data.frame('grouping' = c('Cluster_1', 'Cluster_2', 'Cluster 1'), 'Cells' = c('Cell_1', 'Cell_2', 'Cell_3'))
#' permutation_results <- get_permutation_pvals(counts_table, labels_table, 'transcript_id', 'gene_id', 'grouping', 10000, 0, TRUE, FALSE, 6000, '/user/example/folder/')
#' @export

get_permutation_pvals <- function(transcript_counts_table, cell_labels_table,
                                  transcript_id_colname, gene_id_colname,
                                  cell_labels_colname, permutations = 10000,
                                  cores = 0, do_gene_level_comparisons = FALSE,
                                  return_detailed_pvalue_tables = FALSE,
                                  cutoff = 0.1) {

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

    ### Main body

    # Import and order counts table here first
    filtered_tcount <- transcript_counts_table[order(transcript_counts_table[[gene_id_colname]], transcript_counts_table[[transcript_id_colname]]), ]
    designations <- cell_labels_table

    # INTERNAL SETTINGS
    newconst <- 0.00000001

    full_ids <- filtered_tcount %>% dplyr::select({{transcript_id_colname}}, {{gene_id_colname}})
    colnames(full_ids) <- c('transcript_id', 'gene_id')

    meta <- filtered_tcount %>% dplyr::select({{transcript_id_colname}}, {{gene_id_colname}})
    filtered_tcount[[transcript_id_colname]] <- NULL
    filtered_tcount[[gene_id_colname]] <- NULL

    # Get total counts per transcript, add it to meta
    meta$counts_per_transcript_in_all_cells <- rowSums(filtered_tcount)

    # Split cells into cluster-specific counts tables within a list
    # Note (2023-08-22): But this is just making a vector of cluster names for later
    clusters <- unique(designations[[cell_labels_colname]])[order(unique(designations[[cell_labels_colname]]))]

    # Figure out which column in designations has the grouping labels, and which has the cell labels
    cell_barcodes_col <- ifelse(colnames(designations)[1] == cell_labels_colname, 2, 1)

    # Initialize pvalue_storage_new so you don't have to do it in the loop
    # Only make the first loop's storage though, since we cut down the isoforms in loops 2+
    a <- as.data.frame(cbind(meta[[transcript_id_colname]], meta[[gene_id_colname]], rep(NA, nrow(meta))))
    colnames(a) <- c('transcript_id', 'gene_id', 'pval1')
    sublist <- rep(list(a), length(clusters))
    names(sublist) <- clusters
    pvalue_storage_new <- list(sublist)
    names(pvalue_storage_new) <- c('1')

    # If do_gene_level_comparisons is set, initialize pvalue_storage_gene here
    if (do_gene_level_comparisons) {
        a <- data.frame('gene_id' = unique(meta[[gene_id_colname]]), 'pval1' = rep(NA, length(unique(meta[[gene_id_colname]]))))
        colnames(a) <- c('gene_id', 'pval1')
        sublist <- rep(list(a), length(clusters))
        names(sublist) <- clusters
        pvalue_storage_gene <- rep(list(sublist), permutations)
        names(pvalue_storage_gene) <- as.character(1:permutations)
    }

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

    filtered_tcount_matrix <- as.matrix(filtered_tcount)
    designations <- data.table::as.data.table(designations)
    meta <- data.table::as.data.table(meta)

    # Main permutation loop
    for (loopnum in 1:permutations) {

        message(paste('Loop number', loopnum))
        time_s <- Sys.time()

        # Estimating time left (time1 is first set during the first loop in the if (loopnum == 1) block further down)
        if ((loopnum %% 100 == 0) & (loopnum != 1)) {
            time2 <- Sys.time()
            estimated_time <- round((as.numeric(difftime(time2, time1, units = 'mins')) / 100 * (permutations - loopnum)), 2)
            message(paste0('Estimated time to completion: ', estimated_time, ' minutes'))
            time1 <- time2
        }

        # If we're on cycle 2+, then scramble the grouping column in designations
        if (loopnum != 1) {
            designations[[cell_labels_colname]] <- sample(designations[[cell_labels_colname]])
        }

        # Before the rowSums are calculated, filter the filtered_tcount_matrix rows by removing any genes that DO NOT:
        #   - Pass the initial cutoff at the gene-level (if gene-level is being calculated)
        #   - Associated with isoforms that pass the initial cutoff
        # For the pacbio dataset, this filters the matrix from 5,637 --> 4,723 genes (~20% decrease of genes, but ~10% decrease of isoforms) at a cutoff of 0.05
        # We only need to do this once, since filtered_tcount_matrix is never reset
        if (loopnum == 2) {
            filtered_tcount_matrix <- filtered_tcount_matrix[meta[[gene_id_colname]] %in% genes_to_keep, ]
            meta <- meta[meta[[gene_id_colname]] %in% genes_to_keep, ]
        }

        # Split the counts table into the counts tables per cluster here
        # TODO: The rowSums here is slow, but honestly idk what to do about it, since this is the fastest rowSums available
        # For the pacbio 6 clusters, takes ~4.5 seconds to finish
        # DF implementation
        cluster_counts_list <- list()
        for (cluster in clusters) {
            cells_in_cluster <- unlist(designations[designations$cell_type == cluster, ][, ..cell_barcodes_col])
            column_indexes <- which(colnames(filtered_tcount_matrix) %in% cells_in_cluster)
            a <- filtered_tcount_matrix[, column_indexes]
            cluster_counts_list[[cluster]] <- rowSums(a)
        }


        # Define calc_list here, which will be used for the isoform-level calcs
        calc_list <- list()

        # Retrieve the rowSums and add to meta, with cluster-specific colunm names
        for (groupname in clusters) {
            calc_list[[groupname]] <- data.table::data.table('transcript_id' = meta[[transcript_id_colname]], 'gene_id' = meta[[gene_id_colname]])
            calc_list[[groupname]][['isoform_in']] <- cluster_counts_list[[groupname]]
            calc_list[[groupname]][['isoform_out']] <- unlist(meta$counts_per_transcript_in_all_cells - calc_list[[groupname]][['isoform_in']])
        }

        # Get gene counts per cluster, then you find ng and rg via: ng = all_gene_counts_in_cluster - nt, rg = all_gene_counts_outside_cluster - rt
        # No-merge implementation, since the data.table merge was slightly more costly (279 milliseconds vs the no-merge 171 milliseconds avg to run the whole calc_list generation)
        # Small savings add up over 10,000 permutations
        for (cluster in clusters) {
            per_gene_total <- calc_list[[cluster]][,list(in_aggregated = sum(isoform_in), out_aggregated = sum(isoform_out)), by = 'gene_id']
            repeatnum_vec <- table(calc_list[[cluster]][['gene_id']])
            newtab <- per_gene_total[rep(seq_along(repeatnum_vec), repeatnum_vec), ]
            calc_list[[cluster]][['other_in']] <- newtab$in_aggregated - calc_list[[cluster]][['isoform_in']]
            calc_list[[cluster]][['other_out']] <- newtab$out_aggregated - calc_list[[cluster]][['isoform_out']]
        }

        # If the current loopnum is 1, record the current time (for use later when calculating time per permutation)
        if (loopnum == 1) {
            time1 <- Sys.time()
        }

        # If we're on loops 2+, subset the gene list to genes that had an initial pval < cutoff
        if (do_gene_level_comparisons) {
            calc_list_gene <- list()
            if (loopnum > 1) {
                for (cluster in clusters) {
                    calc_list_gene[[cluster]] <- calc_list[[cluster]][calc_list[[cluster]][['gene_id']] %in% sig_in_at_least_one_cluster_gene, ]
                }
            } else {
                calc_list_gene <- calc_list
            }
        }

        # Do the gene-level pval calculations here
        if (do_gene_level_comparisons) {
            tempresultsgene <- foreach::foreach (iter = 1:length(clusters), .combine = 'cbind', .multicombine = T, .packages = 'parallel') %dopar% {
                cluster_name <- clusters[iter]

                vec_in <- lapply(split(calc_list_gene[[cluster_name]][["isoform_in"]], calc_list_gene[[cluster_name]][["gene_id"]]), FUN = function(x) c(x) + newconst)
                vec_out <- lapply(split(calc_list_gene[[cluster_name]][["isoform_out"]], calc_list_gene[[cluster_name]][["gene_id"]]), FUN = function(x) c(x) + newconst)

                pval_vector <- parallel::mcmapply(chisq.slim.gene.test,
                                                  vec_in,
                                                  vec_out,
                                                  mc.cores = mcmapply_reserved_cores)

                pval_vector
            }

            # If we're in the first loop OR on loop 3+, we can just assign to the pval1 column for pvalue_storage_gene[[loopnum]][[cluster]]
            if ((loopnum == 1) | (loopnum >= 3)) {
                i <- 1
                for (clust_name in names(pvalue_storage_gene[[as.character(loopnum)]])) {
                    pvalue_storage_gene[[as.character(loopnum)]][[clust_name]]$pval1 <- tempresultsgene[, i]
                    i <- i + 1
                }

                # If we're on loop 2, we make the rest of the tables for other iterations in pvalue_storage_gene since we didn't know which isoforms we were gonna calculate until now
                # We also fill out the second loop's pvals per cluster
            } else if (loopnum == 2) {
                temp <- pvalue_storage_gene$`1`
                a <- data.table::data.table('gene_id' = unique(calc_list_gene[[clusters[1]]][['gene_id']]),
                                            'pval1' = rep(NA, length(unique(calc_list_gene[[clusters[1]]][['gene_id']]))))
                sublist <- rep(list(a), length(clusters))
                names(sublist) <- clusters
                pvalue_storage_gene <- rep(list(sublist), permutations)
                names(pvalue_storage_gene) <- c(as.character(1:permutations))
                pvalue_storage_gene[['1']] <- temp

                i <- 1
                for (clust_name in names(pvalue_storage_gene[[as.character(loopnum)]])) {
                    pvalue_storage_gene[[as.character(loopnum)]][[clust_name]]$pval1 <- tempresultsgene[, i]
                    i <- i + 1
                }
            }
        }

        # If we're on loop 2+, we now cut down the calc_list tables to only include isoforms in the vector sig_in_at_least_one_cluster
        # We need to do this AFTER the gene-level calcs are done, since they rely on having all the rows for all isoforms in a gene
        if (loopnum > 1) {
            for (cluster in clusters) {
                calc_list[[cluster]] <- calc_list[[cluster]][calc_list[[cluster]][['transcript_id']] %in% sig_in_at_least_one_cluster, ]
            }
        }


        ### Do the chi-square test ###
        #
        # CONTINGENCY TABLE
        # For each transcript, make a 2x2 contingency table that looks like:
        #
        # Nt    Rt
        # Ng    Rg
        #
        # Where:
        # Nt = Isoform count in cluster (isoform_in)
        # Rt = Isoform count all other cells (isoform_out)
        # Ng = All counts of all other isoforms for gene in cluster (other_in)
        # Rg = All counts of all other isoforms for gene in rest (other_out)

        # Loop over all clusters
        tempresults <- foreach::foreach (iter = 1:length(clusters), .combine = 'cbind', .multicombine = T, .packages = 'parallel') %dopar% {
            cluster_name <- clusters[iter]

            pval_vector <-  parallel::mcmapply(chisq.slim.test,
                                               calc_list[[cluster_name]][['isoform_in']] + newconst,
                                               calc_list[[cluster_name]][['isoform_out']] + newconst,
                                               calc_list[[cluster_name]][['other_in']] + newconst,
                                               calc_list[[cluster_name]][['other_out']] + newconst,
                                               mc.cores = mcmapply_reserved_cores)

            pval_vector

        }

        # If we're in the first loop OR on loop 3+, we can just assign to the pval1 column for pvalue_storage_new[[loopnum]][[cluster]]
        if ((loopnum == 1) | (loopnum >= 3)) {
            i <- 1
            for (clust_name in names(pvalue_storage_new[[as.character(loopnum)]])) {
                pvalue_storage_new[[as.character(loopnum)]][[clust_name]]$pval1 <- tempresults[, i]
                i <- i + 1
            }

            # If we're on loop 2, we make the rest of the tables for other iterations in pvalue_storage_new since we didn't know which isoforms we were gonna calculate until now
            # We also fill out the second loop's pvals per cluster
        } else if (loopnum == 2) {
            temp <- pvalue_storage_new$`1`
            a <- data.table::data.table('transcript_id' = calc_list[[clusters[1]]][['transcript_id']],
                                        'gene_id' = calc_list[[clusters[1]]][['gene_id']],
                                        'pval1' = rep(NA, nrow(calc_list[[clusters[1]]])))
            sublist <- rep(list(a), length(clusters))
            names(sublist) <- clusters
            pvalue_storage_new <- rep(list(sublist), permutations)
            names(pvalue_storage_new) <- c(as.character(1:permutations))
            pvalue_storage_new[['1']] <- temp

            i <- 1
            for (clust_name in names(pvalue_storage_new[[as.character(loopnum)]])) {
                pvalue_storage_new[[as.character(loopnum)]][[clust_name]]$pval1 <- tempresults[, i]
                i <- i + 1
            }
        }

        # If this is the first loop, find the indices of the isoforms that have a p < cutoff, and add them to a vector here
        # For reference, this eliminates ~60% of comparisons for isoforms and ~30% for genes, speeding everything up
        # Split into two parts, depending on whether gene-level comparisons are being calculated or not
        if (loopnum == 1) {

            # If we're on loop 1, store calc_list for retrieval later
            calc_list_for_later <- calc_list

            if (do_gene_level_comparisons) {
                calc_list_gene_for_later <- calc_list_gene

                sig_in_at_least_one_cluster <- c()
                sig_in_at_least_one_cluster_gene <- c()
                for (cluster in clusters) {
                    sig_in_at_least_one_cluster <- c(sig_in_at_least_one_cluster,
                                                     pvalue_storage_new$`1`[[cluster]][pvalue_storage_new$`1`[[cluster]][['pval1']] < cutoff, 'transcript_id'])
                    sig_in_at_least_one_cluster_gene <- c(sig_in_at_least_one_cluster_gene,
                                                          pvalue_storage_gene$`1`[[cluster]][pvalue_storage_gene$`1`[[cluster]][['pval1']] < cutoff, 'gene_id'])
                }
                sig_in_at_least_one_cluster <- unique(sig_in_at_least_one_cluster)
                sig_in_at_least_one_cluster_gene <- unique(sig_in_at_least_one_cluster_gene)
                temp <- meta[[gene_id_colname]][meta[[transcript_id_colname]] %in% sig_in_at_least_one_cluster]
                genes_to_keep <- unique(c(temp, sig_in_at_least_one_cluster_gene))

            } else {
                sig_in_at_least_one_cluster <- c()
                for (cluster in clusters) {
                    sig_in_at_least_one_cluster <- c(sig_in_at_least_one_cluster,
                                                     pvalue_storage_new$`1`[[cluster]][pvalue_storage_new$`1`[[cluster]][['pval1']] < cutoff, 'transcript_id'])
                }
                sig_in_at_least_one_cluster <- unique(sig_in_at_least_one_cluster)
                genes_to_keep <- unique(meta[[gene_id_colname]][meta[[transcript_id_colname]] %in% sig_in_at_least_one_cluster])
            }
        }

        # time_e <- Sys.time()
        # message(paste0('Loop took ', time_e - time_s, ' seconds'))
    }

    parallel::stopCluster(cl)

    ### PVALUE CALCULATIONS
    # Find the p-values for each isoform that we keep, and store them in a list of vectors (per cluster)
    # Also make a list of vectors where every position is an isoform in the tid column of calc_list (for counting later)
    # NOTE (24/11/23): R inherently is only accurate to ~15 digits. Small differences of 1e-10 or below really don't say anything about the permutation, so just round it.
    initial_pvals <- list()
    transcript_pval_count <- list()
    for (cluster in clusters) {
        initial_pvals[[cluster]] <- round(pvalue_storage_new[['1']][[cluster]][which(pvalue_storage_new[['1']][[cluster]][['transcript_id']] %in% calc_list[[clusters[1]]][['transcript_id']]), 'pval1'], 10)
        transcript_pval_count[[cluster]] <- rep(0, length(initial_pvals[[cluster]]))
    }

    # Go through each loop from 2 onwards, and go through each cluster sequentially.
    # Compare to the initial pval, and if less than or equal to, add 1 to position in vec.
    # NOTE (24/11/23): We round the p-vals for loops 2+ here.
    for (loopname in names(pvalue_storage_new)) {
        if (loopname != '1') {
            for (cluster in clusters) {
                transcript_pval_count[[cluster]] <- transcript_pval_count[[cluster]] +
                    (round(pvalue_storage_new[[loopname]][[cluster]][['pval1']], 10) <= initial_pvals[[cluster]])
            }
        }
    }

    # Design the final table for presentation
    # Should have the columns 'transcript_id', 'gene_id', and one column of p-values per cluster
    permutation_pvals <- full_ids
    colnames(permutation_pvals) <- c('transcript_id', 'gene_id')

    # First find the indexes of the tids we kept throughout, so we know where to put the values
    tid_idx <- which(permutation_pvals$transcript_id %in% pvalue_storage_new$`2`[[clusters[1]]][['transcript_id']])

    # Calculate pvals, then make a vec of NAs and insert pvals into correct positions before appending to final table
    for (cluster in clusters) {
        final_pvals <- transcript_pval_count[[cluster]] / (permutations - 1)
        to_insert <- rep(NA, nrow(permutation_pvals))
        for (i in 1:length(tid_idx)) {
            to_insert[tid_idx[i]] <- final_pvals[i]
        }
        permutation_pvals[[cluster]] <- to_insert
    }

    # If gene-level pvals were also calculated, we make another table of gene-level pvals
    if (do_gene_level_comparisons) {

        initial_pvals <- list()
        gene_pval_count <- list()
        for (cluster in clusters) {
            initial_pvals[[cluster]] <- pvalue_storage_gene[['1']][[cluster]][which(pvalue_storage_gene[['1']][[cluster]][['gene_id']] %in% unique(calc_list_gene[[clusters[1]]][['gene_id']])), 'pval1']
            gene_pval_count[[cluster]] <- rep(0, length(initial_pvals[[cluster]]))
        }

        # Go through each loop from 2 onwards, and go through each cluster sequentially.
        # Compare to the initial pval, and if less than or equal to, add 1 to position in vec.
        for (loopname in names(pvalue_storage_gene)) {
            if (loopname != '1') {
                for (cluster in clusters) {
                    gene_pval_count[[cluster]] <- gene_pval_count[[cluster]] +
                        (pvalue_storage_gene[[loopname]][[cluster]][['pval1']] <= initial_pvals[[cluster]])
                }
            }
        }

        # Design the final table for presentation
        # Should have the column 'gene_id' and one column of p-values per cluster
        permutation_pvals_gene <- pvalue_storage_gene$`1`[[clusters[1]]] %>% dplyr::select(gene_id)

        # First find the indexes of the gids we kept throughout, so we know where to put the values
        gid_idx <- which(permutation_pvals_gene$gene_id %in% pvalue_storage_gene$`2`[[clusters[1]]][['gene_id']])

        # Calculate pvals, then make a vec of NAs and insert pvals into correct positions before appending to final table
        for (cluster in clusters) {
            final_pvals <- gene_pval_count[[cluster]] / permutations
            to_insert <- rep(NA, nrow(permutation_pvals_gene))
            for (i in 1:length(tid_idx)) {
                to_insert[gid_idx[i]] <- final_pvals[i]
            }
            permutation_pvals_gene[[cluster]] <- to_insert
        }
    }

    # Make an initial pvals table to return
    initial_pvals_table <- full_ids
    colnames(initial_pvals_table) <- c('transcript_id', 'gene_id')
    for (cluster in clusters) {
        initial_pvals_table[[cluster]] <- pvalue_storage_new$`1`[[cluster]][['pval1']]
    }

    # Make the initial pvals for genes if needed
    if (do_gene_level_comparisons) {
        initial_pvals_table_gene <- pvalue_storage_gene$`1`[[clusters[1]]] %>% dplyr::select(gene_id)
        for (cluster in clusters) {
            initial_pvals_table_gene[[cluster]] <- pvalue_storage_gene$`1`[[cluster]][['pval1']]
        }
    }

    if (!(return_detailed_pvalue_tables)) {
        pvalue_storage_new <- NA
    }

    if (!(do_gene_level_comparisons)) {
        permutation_pvals_gene <- NA
        initial_pvals_table_gene <- NA
        calc_list_gene_for_later <- NA
    }

    to_return <- list('permutation_pvalues' = permutation_pvals, 'first-loop_pvalues' = initial_pvals_table,
                      'pvalue_storage_list' = pvalue_storage_new, 'permutation_pvalues_gene' = permutation_pvals_gene,
                      'first-loop_pvalues_gene' = initial_pvals_table_gene, 'first-loop_contingency_tables' = calc_list_for_later,
                      'first-loop_contingency_tables_gene' = calc_list_gene_for_later)
    return(to_return)
}
