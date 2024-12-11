#' Perform a permutation analysis
#'
#' Performs a permutation analysis using repeated chi-squared measures of significance, while shuffling cell grouping designations
#' between each permutation.
#'
#' @import data.table
#' @import progress
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom data.table %chin%
#' @param transcript_counts_table A dataframe with a column containing transcript IDs, a column with gene IDs for
#' each transcript, and counts for all transcripts in all cells, with cell IDs as subsequent column names.
#' @param cell_labels_table A dataframe with two columns: one listing all cell IDs, and the other listing the group
#' that the cell ID belongs to.
#' @param transcript_id_colname A string corresponding to the column name in `transcript_counts_table` where the transcript IDs are stored.
#' @param gene_id_colname A string corresponding to the column name in `transcript_counts_table` where the gene IDs are stored.
#' @param cell_labels_colname A string corresponding to the column name in `cell_labels_table` where the group information is stored.
#' @param analysis_group_1 This should be set if you want to compare one cell cluster directly against another, or against the rest of the dataset. Defaults to `NA`.
#' Note that one of `analysis_group_1`, `analysis_group_2`, or `run_on_all_groups` needs to be set.
#' @param analysis_group_2 This should be set if you want to compare one cell cluster directly against another, or against the rest of the dataset. Defaults to `NA`.
#' Note that one of `analysis_group_1`, `analysis_group_2`, or `run_on_all_groups` needs to be set.
#' @param run_on_all_groups This should be set to `TRUE` if you want to compare every cluster against the rest of the dataset. Defaults to `FALSE`.
#' Note that one of `analysis_group_1`, `analysis_group_2`, or `run_on_all_groups` needs to be set.
#' @param permutations An integer corresponding to the number of permutations that should be run. Default 10000. Note that the first p-value calculation
#' is performed on the original set of data, and label shuffling occurs from permutation 2 and onwards. If only 1 permmutation is chosen, the
#' output p-values will be the chi-squared p-values.
#' @param cores An integer representing the number of cores that the function should attempt to use. If left undefined, set to 0, or set
#' to more cores than are available, defaults to all available cores - 1. A good starting point is to allocate at least N cores, where
#' `N = (number of cell groupings * 2)` if `runall = TRUE`, or otherwise `N = all available cores - 1`. The function will run at varying degrees
#' of reduced speed if this threshold isn't met.
#' @param do_gene_level_comparisons If set to `TRUE`, also runs a permutation analysis on gene level differences. A significant value indicates that
#' there is a difference in transcript proportions visible at the gene level, but does not specify which transcripts show the change within the gene.
#' Defaults to `TRUE`. Also adds an additional step to the cutoff filtering (for more detail, see the isopod vignette).
#' @param return_detailed_pvalue_tables If set to `TRUE`, an additonal list of dataframes will be returned containing additional information on
#' p-values generated over the permutations. If you'd like to look into p-value distributions over the permutations, this must be set to `TRUE`.
#' Note that if this is set to `TRUE`, the resulting object generated will be very large (normally ~100GB), and therefore this option is only
#' recommended for further (usually unnecessary) inspection.
#' @param report_adjusted_pvalues If `TRUE`, p-values will be corrected via a Benjamini & Hochberg correction, and unadjusted p-values will
#' be reported in a separate table.
#' @param cutoff A `double` (decimal) indicating the initial p-value cutoff that isoforms (and genes if gene-level comparison is enabled) must pass to be included
#' in subsequent permutations. This cutoff applies for the first permutation only, and potentially increases the speed of each permutation by filtering down
#' data before calculations begin. With higher cutoffs, less data is filtered. To disable filtering, set the cutoff to 1. Defaults to `0.1`.
#' @return A `list` containing:\n
#' - `permutation_pvalues`: A dataframe of permutation p-values.\n
#' - `first-loop_pvalues`: A dataframe of the first-loop chi-squared p-values.\n
#' - `first-loop_contingency_tables`: A dataframe of first-loop contingency tables for each isoform.\n
#' - `first-loop_contingency_tables_gene`: A dataframe of first-loop contingency table values for each gene.\n
#' - `transcripts_filtered_from_cutoff`: A list of transcripts filtered due to the cutoff value.\n
#' - `genes_filtered_from_cutoff`: If `do_gene_level_comparisons = TRUE`, a list of genes filtered due to the cutoff value.\n
#' - `permutation_pvalues`: If `do_gene_level_comparisons = TRUE`, a dataframe of permutation p-values for each gene.\n
#' - `first-loop_pvalues_gene` : If `do_gene_level_comparisons = TRUE`, a dataframe of the first-loop chi-squared p-values for each gene. \n
#' - `odds_ratio_table`: If `do_gene_level_comparisons = TRUE`, the calculated odds-ratio for each isoform. Requires gene-level data to calculate.\n
#' Note that the columns `other_in` and `other_out` are not used for gene-level calculations.\n
#' - `pvalue_storage_list`: If `return_detailed_pvalue_tables = TRUE`, a dataframe containing isoform-level p-values at each step of the permutation.\n
#' - `unadjusted_permutation_pvalues`: If `report_adjusted_pvalues = TRUE`, the unadjusted permutation p-values for each isoform.\n
#' - `unadjusted_first-loop_pvalues`: If `report_adjusted_pvalues = TRUE`, the unadjusted permutation p-values for each isoform.\n
#' - `unadjusted_permutation_pvalues_gene`: If `report_adjusted_pvalues = TRUE` and `do_gene_level_comparisons = TRUE`, the unadjusted permutation p-values for each gene.\n
#' - `unadjusted_first-loop_pvalues_gene`: If `report_adjusted_pvalues = TRUE` and `do_gene_level_comparisons = TRUE`, the first-loop chi-squared p-values for each gene.
#'
#' @examples
#' # Note that transcript and gene IDs can be character strings instead of numbers.
#' counts_table <- data.frame('transcript_id' = c(1, 2, 3), 'gene_id' = c(1, 1, 2), 'Cell_1' = c(0, 1, 10), 'Cell_2' = c(10, 2, 5), 'Cell_3' = c(2, 5, 1))
#' labels_table <- data.frame('grouping' = c('Cluster_1', 'Cluster_2', 'Cluster 1'), 'Cells' = c('Cell_1', 'Cell_2', 'Cell_3'))
#' permutation_results <- get_permutation_pvals(counts_table, labels_table, 'transcript_id', 'gene_id', 'grouping',
#'                                              run_on_all_groups = TRUE, permutations = 10000, cores = 0,
#'                                              do_gene_level_comparisons = TRUE, return_detailed_pvalue_tables = FALSE,
#'                                              report_adjusted_pvalues = TRUE, cutoff = 1)
#' @export

get_permutation_pvals <- function(transcript_counts_table, cell_labels_table,
                                  transcript_id_colname, gene_id_colname,
                                  cell_labels_colname,
                                  analysis_group_1 = NA, analysis_group_2 = NA,
                                  run_on_all_groups = FALSE, permutations = 10000,
                                  cores = 0, do_gene_level_comparisons = TRUE,
                                  return_detailed_pvalue_tables = FALSE, report_adjusted_pvalues = TRUE,
                                  cutoff = 0.1) {

    ### Helper functions
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

    # Check analysis type inputs and determine the runmode of the function
    # The variable runmode can be '1v1', '1vrest', or 'all'
    check_analysis_inputs <- function(analysis_group_1, analysis_group_2, run_on_all_groups)
    {
        g1na <- is.na(analysis_group_1)
        g2na <- is.na(analysis_group_2)
        runall <- run_on_all_groups

        if (g1na & g2na & !runall) {
            stop(paste0('The arguments analysis_group_1, analysis_group_2, and run_on_all_groups have not been set.\n',
                        'Depending on what you want to compare, set the variables as listed below.\n',
                        '-------------------------\n',
                        'Analysing every cell group against all non-group cells (computationally intensive):\n',
                        '\t- Set run_on_all_groups to TRUE\n',
                        'One cell group specifically against another group:\n',
                        '\t- Set BOTH analysis_group_1 and analysis_group_2\n',
                        'One cell group against all non-group cells:\n',
                        '\t- Set EITHER analysis_group_1 or analysis_group_2'))
        } else if (!runall) {  # If run_on_all_groups is FALSE, then check g1na and g2na at once
            if (xor(!g1na, !g2na)) {
                runmode <- '1vrest'
                if (!g2na) {  # If only group 2 is specified, assign it to group 1 instead and make group 2 NA (easier downstream)
                    g1 <- analysis_group_2
                    g2 <- NA
                } else if (!g1na) {
                    g1 <- analysis_group_1
                    g2 <- NA
                }
            } else {
                runmode <- '1v1'
                g1 <- analysis_group_1
                g2 <- analysis_group_2
            }
        } else {  # If runall is TRUE, check whether g1na g2na are true. If so, throw an error.
            if (!g1na & !g2na) {
                stop(paste0('run_on_all_groups AND both analysis groups have been set.\n',
                            'Either disable run_on_all_groups, or remove the analysis group arguments.'))
            } else if (!g1na & g2na) {
                stop(paste0('run_on_all_groups AND analysis_group_1 have been set.\n',
                            'Either disable run_on_all_groups, or remove the analysis group arguments.'))
            } else if (g1na & !g2na) {
                stop(paste0('run_on_all_groups AND analysis_group_2 have been set.\n',
                            'Either disable run_on_all_groups, or remove the analysis group arguments.'))
            } else {
                runmode <- 'all'
                g1 <- NA
                g2 <- NA
            }
        }
        return(list('runmode' = runmode,
                    'g1' = g1,
                    'g2' = g2))
    }


    ### Main body

    results_check <- check_analysis_inputs(analysis_group_1, analysis_group_2, run_on_all_groups)
    runmode <- results_check$runmode
    analysis_group_1 <- results_check$g1
    analysis_group_2 <- results_check$g2

    # Check that we're using a list that we made from the filtering function. If so, retain ONLY the counts table
    if (length(names(transcript_counts_table)) == 2 & identical(names(transcript_counts_table), c("counts_table", "list_of_collapsed_isoforms"))) {
        transcript_counts_table <- transcript_counts_table$counts_table
    }

    # Import and order counts table here first
    filtered_tcount <- transcript_counts_table[order(transcript_counts_table[[gene_id_colname]], transcript_counts_table[[transcript_id_colname]]), ]
    designations <- cell_labels_table

    # INTERNAL SETTINGS
    newconst <- 0.1

    full_ids <- filtered_tcount %>% dplyr::select({{transcript_id_colname}}, {{gene_id_colname}})
    colnames(full_ids) <- c('transcript_id', 'gene_id')

    meta <- filtered_tcount %>% dplyr::select({{transcript_id_colname}}, {{gene_id_colname}})
    filtered_tcount[[transcript_id_colname]] <- NULL
    filtered_tcount[[gene_id_colname]] <- NULL

    # Get total counts per transcript, add it to meta
    meta$counts_per_transcript_in_all_cells <- rowSums(filtered_tcount)

    # Split cells into cluster-specific counts tables within a list
    # Note (2023-08-22): But this is just making a vector of cluster names for later
    # For '1v1' runmode, this is a vector containing the two clusters to analyse.
    if (runmode == 'all') {
        clusters <- unique(designations[[cell_labels_colname]])[order(unique(designations[[cell_labels_colname]]))]
    } else {
        clusters <- c(analysis_group_1)
    }

    # Figure out which column in designations has the grouping labels, and which has the cell labels
    cell_barcodes_col <- ifelse(colnames(designations)[1] == cell_labels_colname, 2, 1)

    # Initialize pvalue_storage_new so you don't have to do it in the loop
    # Only make the first loop's storage though, since we cut down the isoforms in loops 2+
    # IF runmode is 1vrest or 1v1, the pvalue_storage_new has [[name_of_first_clust]]
    # IF runmode is all, the pvalue_storage_new has [[cluster_name]] for each cluster in clusters
    a <- as.data.frame(cbind(meta[[transcript_id_colname]], meta[[gene_id_colname]], rep(NA, nrow(meta))))
    colnames(a) <- c('transcript_id', 'gene_id', 'pval1')
    if (runmode == 'all') {
        sublist <- rep(list(a), length(clusters))
        names(sublist) <- clusters
    } else {
        sublist <- list(a)
        names(sublist) <- clusters
    }
    pvalue_storage_new <- list(sublist)
    names(pvalue_storage_new) <- c('1')

    # If do_gene_level_comparisons is set, initialize pvalue_storage_gene here
    # Naming convention is the same as pvalue_storage_new for isoforms
    if (do_gene_level_comparisons) {
        a <- data.frame('gene_id' = unique(meta[[gene_id_colname]]), 'pval1' = rep(NA, length(unique(meta[[gene_id_colname]]))))
        colnames(a) <- c('gene_id', 'pval1')
        if (runmode == 'all') {
            sublist <- rep(list(a), length(clusters))
            names(sublist) <- clusters
        } else {
            sublist <- list(a)
            names(sublist) <- clusters
        }
        pvalue_storage_gene <- rep(list(sublist), permutations)
        names(pvalue_storage_gene) <- as.character(1:permutations)
    }

    # For parallelizing foreach later on, make the socket first
    if (cores == 0) {
        cores = parallelly::availableCores() - 1
    }

    # Check that there's enough cores to run everything efficiently.
    # If there's equal or less cores than there are clusters, we need to make sure things can actually run, speed be damned.
    # Note (2024-06-05): Added core allocation for '1v1' and '1vrest' runmodes.
    if (runmode == 'all') {
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
    } else {
        mcmapply_reserved_cores <- max(1, cores - 1)
        cl <- parallel::makeCluster(1)
        doParallel::registerDoParallel(cl)
    }

    filtered_tcount_matrix <- as.matrix(filtered_tcount)
    designations <- data.table::as.data.table(designations)
    meta <- data.table::as.data.table(meta)

    pb <- progress::progress_bar$new(
        format = "  Running permutation analysis [:bar] :current/:total (:percent), eta: :eta , elapsed: :elapsed",
        total = permutations, clear = FALSE, width = 100)
    pb$tick(0)

    # Some cluster environments need explicit exporting of .libPaths() to each worker node
    paths <- .libPaths()
    parallel::clusterExport(cl = cl, varlist = c('paths'), envir = environment())
    parallel::clusterEvalQ(cl, .libPaths(paths))

    # Main permutation loop
    for (loopnum in 1:permutations) {

        # If we're on cycle 2+, then scramble the grouping column in designations
        if (loopnum != 1) {
            designations[[cell_labels_colname]] <- sample(designations[[cell_labels_colname]])
        }

        # Split the counts table into the counts tables per cluster here
        # The rowSums here can be slow, but honestly idk what to do about it since this is the fastest rowSums available
        # For the pacbio 6 clusters, takes ~4.5 seconds to finish. (Note (2024-06-05): should be fixed)
        # DF implementation
        # For a 1v1 implementation, only 2 clusters are calculated.
        # For a 1vrest, only 1 cluster is calculated, then counts of all others can be calculated by total - cluster
        cluster_counts_list <- list()
        if (runmode == '1v1') {
            clusters <- c(analysis_group_1, analysis_group_2)
        }
        for (cluster in clusters) {
            cells_in_cluster <- unlist(designations[designations[[cell_labels_colname]] == cluster, ][, ..cell_barcodes_col])
            column_indexes <- which(colnames(filtered_tcount_matrix) %in% cells_in_cluster)
            a <- filtered_tcount_matrix[, column_indexes]
            cluster_counts_list[[cluster]] <- rowSums(a)
        }
        if (runmode == '1v1') {
            clusters <- c(analysis_group_1)
        }

        # Define calc_list here, which will be used for the isoform-level calcs
        calc_list <- list()

        # Retrieve the rowSums and add to meta, with cluster-specific colunm names
        # If we're doing 1v1, the isoform_out is just the counts in the other cluster
        # If we're not doing 1v1, isoform_out is calculated by total transcript count - count in cluster
        for (groupname in clusters) {
            calc_list[[groupname]] <- data.table::data.table('transcript_id' = meta[[transcript_id_colname]], 'gene_id' = meta[[gene_id_colname]])
            calc_list[[groupname]][['isoform_in']] <- cluster_counts_list[[groupname]]
            if (runmode != '1v1') {
                calc_list[[groupname]][['isoform_out']] <- unlist(meta$counts_per_transcript_in_all_cells - calc_list[[groupname]][['isoform_in']], use.names = F)
            } else {
                calc_list[[groupname]][['isoform_out']] <- cluster_counts_list[[analysis_group_2]]
            }
        }

        # Get gene counts per cluster, then you find ng and rg via: ng = all_gene_counts_in_cluster - nt, rg = all_gene_counts_outside_cluster - rt
        # No-merge implementation, since the data.table merge was slightly more costly (279 milliseconds vs the no-merge 171 milliseconds avg to run the whole calc_list generation)
        # Small savings add up over 10,000 permutations
        # This works automatically for the 1v1 runmode as well, since it just sums isoform in and isoform out by gene id.
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

        # If we're on loops 2+, subset the gene list to genes that had an initial pval =< cutoff
        # NOTE (15/12/2023): Changed to filter dynamically for each cluster
        if (do_gene_level_comparisons) {
            calc_list_gene <- list()
            if (loopnum > 1) {
                for (cluster in clusters) {
                    calc_list_gene[[cluster]] <- calc_list[[cluster]][calc_list[[cluster]][['gene_id']] %chin% list_of_genes_to_keep[[cluster]], ]
                }
            } else {
                calc_list_gene <- calc_list
            }
        }

        # See if any genes are kept here
        if (loopnum == 2) {
            templistresult <- lapply(list_of_genes_to_keep, function(x) {identical(x, character(0))})
            if (TRUE %in% templistresult) {
                stop(paste0('No genes passed the first-loop chi-squared cutoff of ', cutoff, '\n',
                            'We recommend setting the cutoff argument to cutoff=1, then rerunning this.'))
            }
        }

        # Do the gene-level pval calculations here
        # NOTE (15/12/2023): We round the pvalues to 13 digits here, both to save space later on and because R doesn't accurately calculate beyond ~15 decimals
        # Honestly if your statstical test depends on decimals beyond 10^-10, is it really a good test for this sort of implementation in the first place?
        if (do_gene_level_comparisons) {
            tempresultsgene <- foreach::foreach (iter = 1:length(clusters), .multicombine = T, .packages = c('parallel')) %dopar% {

                cluster_name <- clusters[iter]

                vec_in <- lapply(split(calc_list_gene[[cluster_name]][["isoform_in"]], calc_list_gene[[cluster_name]][["gene_id"]]), FUN = function(x) c(x) + newconst)
                vec_out <- lapply(split(calc_list_gene[[cluster_name]][["isoform_out"]], calc_list_gene[[cluster_name]][["gene_id"]]), FUN = function(x) c(x) + newconst)

                pval_vector <- parallel::mcmapply(chisq.slim.gene.test,
                                                  vec_in,
                                                  vec_out,
                                                  mc.cores = mcmapply_reserved_cores)

                pval_vector <- round(pval_vector, digits = 13)
                pval_table_gene <- as.data.frame(matrix(NA, nrow = length(pval_vector), ncol = 2))
                colnames(pval_table_gene) <- c('gene_id', 'pval1')
                pval_table_gene$gene_id <- names(pval_vector)
                pval_table_gene$pval1 <- unname(pval_vector)
                pval_table_gene
            }

            # If we're in the first loop OR on loop 3+, we can just assign to the pval1 column for pvalue_storage_gene[[loopnum]][[cluster]]
            if ((loopnum == 1) | (loopnum >= 3)) {
                i <- 1
                for (cluster in clusters) {
                    pvalue_storage_gene[[as.character(loopnum)]][[cluster]]$pval1 <- tempresultsgene[[i]]$pval1
                    i <- i + 1
                }

                # If we're on loop 2, we make the rest of the tables for other iterations in pvalue_storage_gene since we didn't know which isoforms we were gonna calculate until now
                # We also fill out the second loop's pvals per cluster
                # TODO: This either needs to be changed to generate unique tables per cluster from permutation 2+ (since isoforms calculated per cluster changes),
                #       or we need to insert NAs for isoforms filtered out per cluster to the p-value vectors.
                # IMPLEMENTATION NOTE (15/12/2023): Neither of the above. Instead, a unique table for each cluster is made, filled with only the genes kept in that cluster
            } else if (loopnum == 2) {
                temp <- pvalue_storage_gene$`1`
                sublist <- vector(mode = 'list', length(clusters))
                names(sublist) <- clusters
                for (cluster in clusters) {
                    sublist[[cluster]] <- data.table::data.table('gene_id' = unique(calc_list_gene[[cluster]][['gene_id']]),
                                                                 'pval1' = rep(NA, length(unique(calc_list_gene[[cluster]][['gene_id']]))))
                }

                pvalue_storage_gene <- rep(list(sublist), permutations)
                names(pvalue_storage_gene) <- c(as.character(1:permutations))
                pvalue_storage_gene[['1']] <- temp

                i <- 1
                for (cluster in names(pvalue_storage_gene[[as.character(loopnum)]])) {
                    pvalue_storage_gene[[as.character(loopnum)]][[cluster]]$pval1 <- tempresultsgene[[i]]$pval1
                    i <- i + 1
                }
            }
        }

        # If we're on loop 2+, we now cut down the calc_list tables to only include isoforms in the vector sig_in_at_least_one_cluster
        # We need to do this AFTER the gene-level calcs are done, since they rely on having all the rows for all isoforms in a gene
        # This is actually an easy fix! We're already looping over all clusters, we can just force it to filter from a list of vectors instead of one vector
        if (loopnum > 1) {
            for (cluster in clusters) {
                calc_list[[cluster]] <- calc_list[[cluster]][calc_list[[cluster]][['transcript_id']] %in% list_of_isos_to_keep[[cluster]], ]
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

        # Loop over all clusters (RETURNS A LIST OF TABLES)
        # TODO: This only needs to return a list of tables for loop 2. From loop 3+, it can just return a list of vectors that get slotted into the pre-made p-value tables
        tempresults <- foreach::foreach (iter = 1:length(clusters), .multicombine = T, .packages = c('parallel')) %dopar% {

            cluster_name <- clusters[iter]

            pval_vector <-  parallel::mcmapply(chisq.slim.test,
                                               calc_list[[cluster_name]][['isoform_in']] + newconst,
                                               calc_list[[cluster_name]][['isoform_out']] + newconst,
                                               calc_list[[cluster_name]][['other_in']] + newconst,
                                               calc_list[[cluster_name]][['other_out']] + newconst,
                                               mc.cores = mcmapply_reserved_cores)

            pval_vector <- round(pval_vector, digits = 10)
            pval_table <- as.data.frame(matrix(NA, nrow = length(pval_vector), ncol = 3))
            colnames(pval_table) <- c('transcript_id', 'gene_id', 'pval1')
            pval_table$transcript_id <- calc_list[[cluster_name]][['transcript_id']]
            pval_table$gene_id <- calc_list[[cluster_name]][['gene_id']]
            pval_table$pval1 <- pval_vector
            pval_table

        }

        # If we're in the first loop OR on loop 3+, we can just assign to the pval1 column for pvalue_storage_new[[loopnum]][[cluster]]
        if ((loopnum == 1) | (loopnum >= 3)) {
            i <- 1
            for (cluster in clusters) {
                pvalue_storage_new[[as.character(loopnum)]][[cluster]]$pval1 <- tempresults[[i]]$pval1
                if (return_detailed_pvalue_tables) {
                    pvalue_storage_new[[as.character(loopnum)]][[cluster]]$isoform_in <- calc_list[[cluster]][['isoform_in']]
                    pvalue_storage_new[[as.character(loopnum)]][[cluster]]$isoform_out <- calc_list[[cluster]][['isoform_out']]
                    pvalue_storage_new[[as.character(loopnum)]][[cluster]]$other_in <- calc_list[[cluster]][['other_in']]
                    pvalue_storage_new[[as.character(loopnum)]][[cluster]]$other_out <- calc_list[[cluster]][['other_out']]
                }
                i <- i + 1
            }

            # If we're on loop 2, we make the rest of the tables for other iterations in pvalue_storage_new since we didn't know which isoforms we were gonna calculate until now
            # We also fill out the second loop's pvals per cluster
        } else if (loopnum == 2) {
            temp <- pvalue_storage_new$`1`
            sublist <- vector(mode = 'list', length(clusters))
            names(sublist) <- clusters
            for (cluster in clusters) {
                if (return_detailed_pvalue_tables) {
                    sublist[[cluster]] <- data.table::data.table('transcript_id' = calc_list[[cluster]][['transcript_id']],
                                                                 'gene_id' = calc_list[[cluster]][['gene_id']],
                                                                 'pval1' = rep(NA, length(calc_list[[cluster]][['transcript_id']])),
                                                                 'isoform_in' = rep(NA, length(calc_list[[cluster]][['transcript_id']])),
                                                                 'isoform_out' = rep(NA, length(calc_list[[cluster]][['transcript_id']])),
                                                                 'other_in' = rep(NA, length(calc_list[[cluster]][['transcript_id']])),
                                                                 'other_out' = rep(NA, length(calc_list[[cluster]][['transcript_id']])))
                } else {
                    sublist[[cluster]] <- data.table::data.table('transcript_id' = calc_list[[cluster]][['transcript_id']],
                                                                 'gene_id' = calc_list[[cluster]][['gene_id']],
                                                                 'pval1' = rep(NA, length(calc_list[[cluster]][['transcript_id']])))
                }
            }

            pvalue_storage_new <- rep(list(sublist), permutations)
            names(pvalue_storage_new) <- c(as.character(1:permutations))
            pvalue_storage_new[['1']] <- temp

            i <- 1
            for (cluster in names(pvalue_storage_new[[as.character(loopnum)]])) {  # Fill out pvalue_storage_new dataframe for each cluster
                pvalue_storage_new[[as.character(loopnum)]][[cluster]]$pval1 <- tempresults[[i]]$pval1

                if (return_detailed_pvalue_tables) {
                    pvalue_storage_new[[as.character(loopnum)]][[cluster]]$isoform_in <- calc_list[[cluster]][['isoform_in']]
                    pvalue_storage_new[[as.character(loopnum)]][[cluster]]$isoform_out <- calc_list[[cluster]][['isoform_out']]
                    pvalue_storage_new[[as.character(loopnum)]][[cluster]]$other_in <- calc_list[[cluster]][['other_in']]
                    pvalue_storage_new[[as.character(loopnum)]][[cluster]]$other_out <- calc_list[[cluster]][['other_out']]
                }

                i <- i + 1
            }
        }

        pb$tick(1)

        # If this is the first loop, find the indices of the isoforms that have a p < cutoff, and add them to a vector here
        # For reference, this eliminates the majority of isoforms, speeding everything up
        # Split into two parts, depending on whether gene-level comparisons are being calculated or not
        # We need to make a list of vectors, so we can remove isoforms on a per-cluster basis
        if (loopnum == 1) {

            # If we're on loop 1, store calc_list for retrieval later
            calc_list_for_later <- calc_list

            if (do_gene_level_comparisons) {
                calc_list_gene_for_later <- calc_list_gene

                # Isos to keep needs to be as a vector of strings, since it's too slow/complicated to make the T/F vector due to previous cluster-based filtering with genes
                # Genes to keep needs to be as a vector of strings, since it's too slow/complicated to make the T/F vector due to repeats
                # Genes to keep should also include genes that may not have passed the initial cutoff, but are associated with isoforms that need to be kept
                list_of_isos_to_keep <- vector(mode = 'list', length = length(clusters))
                list_of_genes_to_keep <- vector(mode = 'list', length = length(clusters))
                names(list_of_isos_to_keep) <- clusters
                names(list_of_genes_to_keep) <- clusters

                for (cluster in clusters) {
                    list_of_isos_to_keep[[cluster]] <- pvalue_storage_new$`1`[[cluster]][pvalue_storage_new$`1`[[cluster]][['pval1']] <= cutoff, 'transcript_id']
                    list_of_genes_to_keep[[cluster]] <- unique(pvalue_storage_gene$`1`[[cluster]][pvalue_storage_gene$`1`[[cluster]][['pval1']] <= cutoff, 'gene_id'])
                    list_of_genes_to_keep[[cluster]] <- unique(list_of_genes_to_keep[[cluster]],
                                                               pvalue_storage_new$`1`[[cluster]][['transcript_id']][list_of_isos_to_keep[[cluster]]])
                }

                # If we're not calculating gene-level significance, we just keep genes that are associated with the isoforms we keep
            } else {
                list_of_isos_to_keep <- vector(mode = 'list', length = length(clusters))
                list_of_genes_to_keep <- vector(mode = 'list', length = length(clusters))
                names(list_of_isos_to_keep) <- clusters
                names(list_of_genes_to_keep) <- clusters

                for (cluster in clusters) {
                    list_of_isos_to_keep[[cluster]] <- pvalue_storage_new$`1`[[cluster]][pvalue_storage_new$`1`[[cluster]][['pval1']] <= cutoff, 'transcript_id']
                    list_of_genes_to_keep[[cluster]] <- unique(pvalue_storage_new$`1`[[cluster]][['gene_id']][list_of_isos_to_keep[[cluster]]])
                }
            }
        }
    }

    parallel::stopCluster(cl)

    cat('\nFinished permutations\n\nCalculating permutation p-values...\n')

    ### PVALUE CALCULATIONS
    # Find the p-values for each isoform that we keep, and store them in a list of vectors (per cluster)
    # Also make a list of vectors where every position is an isoform in the tid column of calc_list (for counting later)
    initial_pvals <- list()
    transcript_pval_count <- list()
    for (cluster in clusters) {
        initial_pvals[[cluster]] <- pvalue_storage_new[['1']][[cluster]][['pval1']][pvalue_storage_new[['1']][[cluster]][['transcript_id']] %in% pvalue_storage_new[['2']][[cluster]][['transcript_id']]]
        transcript_pval_count[[cluster]] <- rep(0, length(initial_pvals[[cluster]]))
    }

    # Go through each loop from 2 onwards, and go through each cluster sequentially.
    # Compare to the initial pval, and if less than or equal to, add 1 to position in vec.
    for (loopname in names(pvalue_storage_new)) {
        if (loopname != '1') {
            for (cluster in clusters) {
                transcript_pval_count[[cluster]] <- transcript_pval_count[[cluster]] +
                    (pvalue_storage_new[[loopname]][[cluster]][['pval1']] <= initial_pvals[[cluster]])
            }
        }
    }

    # Design the final table for presentation
    # Should have the columns 'transcript_id', 'gene_id', and one column of p-values per cluster
    permutation_pvals <- full_ids
    colnames(permutation_pvals) <- c('transcript_id', 'gene_id')

    # First find the indexes of the tids we kept throughout, so we know where to put the values (per cluster)
    # Calculate pvals, then make a vec of NAs and insert pvals into correct positions before appending to final table
    for (cluster in clusters) {
        tid_idx <- which(permutation_pvals$transcript_id %in% pvalue_storage_new[['2']][[cluster]][['transcript_id']])
        final_pvals <- (transcript_pval_count[[cluster]] + 1) / permutations
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
            initial_pvals[[cluster]] <- pvalue_storage_gene[['1']][[cluster]][which(pvalue_storage_gene[['1']][[cluster]][['gene_id']] %in% unique(calc_list_gene[[cluster]][['gene_id']])), 'pval1']
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
        # Calculate pvals, then make a vec of NAs and insert pvals into correct positions before appending to final table
        for (cluster in clusters) {
            gid_idx <- which(permutation_pvals_gene$gene_id %in% pvalue_storage_gene$`2`[[clusters[1]]][['gene_id']])
            final_pvals <- (gene_pval_count[[cluster]] + 1) / (permutations)
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

    # Make the odds ratio table
    # Note that you can't do this for genes, since you don't make a 2x2 contingency
    ortab <- full_ids
    for (cluster in clusters) {
        temptab <- calc_list_for_later[[cluster]]
        ortab[[cluster]] <- ((temptab$isoform_in + newconst) * (temptab$other_out + newconst)) / ((temptab$isoform_out + newconst) * (temptab$other_in + newconst))
    }


    if (!(return_detailed_pvalue_tables)) {
        pvalue_storage_new <- NA
    }

    if (!(do_gene_level_comparisons)) {
        permutation_pvals_gene <- NA
        initial_pvals_table_gene <- NA
        calc_list_gene_for_later <- NA
        ortab_gene <- NA
    }

    transcripts_filtered_from_cutoff <- list()
    genes_filtered_from_cutoff <- list()
    for (cluster in clusters) {
        transcripts_filtered_from_cutoff[[cluster]] <- permutation_pvals$transcript_id[is.na(permutation_pvals[[cluster]])]

        if (do_gene_level_comparisons) {
            genes_filtered_from_cutoff[[cluster]] <- permutation_pvals_gene$gene_id[is.na(permutation_pvals_gene[[cluster]])]
            permutation_pvals_gene[[cluster]] <- ifelse(is.na(permutation_pvals_gene[[cluster]]), 1, permutation_pvals_gene[[cluster]])
        }

        permutation_pvals[[cluster]] <- ifelse(is.na(permutation_pvals[[cluster]]), 1, permutation_pvals[[cluster]])

        if (report_adjusted_pvalues) {
            permutation_pvals[[cluster]] <- p.adjust(permutation_pvals[[cluster]], method = 'fdr')
            initial_pvals_table[[cluster]] <- p.adjust(initial_pvals_table[[cluster]], method = 'fdr')

            if (do_gene_level_comparisons) {
                permutation_pvals_gene[[cluster]] <- p.adjust(permutation_pvals_gene[[cluster]], method = 'fdr')
                initial_pvals_table_gene[[cluster]] <- p.adjust(initial_pvals_table_gene[[cluster]], method = 'fdr')
            }
        }
    }

    # Replace NAs with p=1 here, and keep a list of transcript_ids and gene_ids kept
    # If report_adjusted_pvalues, then adjust them here and keep a copy of the raw pval tables.
    # FDR-adjust p-values here if adjust is true
    if (report_adjusted_pvalues) {
        permutation_pvals_noadjust <- permutation_pvals
        permutation_pvals_gene_noadjust <- permutation_pvals_gene
        initial_pvals_noadjust <- initial_pvals_table
        initial_pvals_gene_noadjust <- initial_pvals_table_gene
    }

    to_return <- list('permutation_pvalues' = permutation_pvals, 'first-loop_pvalues' = initial_pvals_table,
                      'pvalue_storage_list' = pvalue_storage_new, 'permutation_pvalues_gene' = permutation_pvals_gene,
                      'first-loop_pvalues_gene' = initial_pvals_table_gene, 'first-loop_contingency_tables' = calc_list_for_later,
                      'first-loop_contingency_tables_gene' = calc_list_gene_for_later, 'odds_ratio_table' = ortab,
                      'transcripts_filtered_from_cutoff' = transcripts_filtered_from_cutoff)

    if (do_gene_level_comparisons) {
        to_return[['genes_filtered_from_cutoff']] <- genes_filtered_from_cutoff
    }

    if (report_adjusted_pvalues) {
        to_return[['unadjusted_permutation_pvalues']] <- permutation_pvals_noadjust
        to_return[['unadjusted_first-loop_pvalues']] <- initial_pvals_noadjust
        if (do_gene_level_comparisons) {
            to_return[['unadjusted_permutation_pvalues_gene']] <- permutation_pvals_gene_noadjust
            to_return[['unadjusted_first-loop_pvalues_gene']] <- initial_pvals_gene_noadjust
        }
    }

    cat('Finished DTU permutation analysis')

    return(to_return)
}
