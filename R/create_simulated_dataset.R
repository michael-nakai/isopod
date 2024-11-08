#' @importFrom checkmate checkPathForOutput
create_simulated_dataset <- function(dataset_name,
                                     number_of_genes = 2000,
                                     cells_in_groups = c(200, 200),
                                     isoform_number_distribution = data.frame('isoforms_in_gene' = 1:10, 'proportion' = rep(0.1, 10)),
                                     gamma_shape = 1.5,
                                     gamma_scale = 0.5,
                                     proportion_of_genes_with_DTU = 0.3,
                                     isoforms_to_swap = c(1, 2),
                                     outlier_gene_proportion = 0,
                                     gene_outliers_lognormal_location = 2.5,
                                     gene_outliers_lognormal_scale = 0.3,
                                     proportion_of_genes_with_DGE = 0,
                                     DGE_lognormal_location = 0.8,
                                     DGE_lognormal_scale = 0.8,
                                     counts_distribution_type = 'negative_binomial',
                                     negative_binomial_low_bound = 0.2,
                                     negative_binomial_high_bound = 0.7,
                                     proportion_of_cells_with_inflated_counts = 0,
                                     cell_outliers_lognormal_location = 2.3,
                                     cell_outliers_lognormal_scale = 0.4,
                                     save_filepath = NA,
                                     verbose = F,
                                     skip_distribution_dataframe_check = F) {



    ### ARGUMENT CHECKS
    ### Do ALL the argument checks here first, and error out early if something's wrong

    # 1. save_filepath valid?
    if (!(is.na(save_filepath))) {
        pathcheck <- checkPathForOutput(save_filepath)
        if (is.character(pathcheck)) {
            stop(pathcheck)
        }
    }

    # 2. proportion arguments between 0 and 1?
    if ((proportion_of_genes_with_DTU < 0) | (proportion_of_genes_with_DTU > 1)) {
        stop('The proportion of genes with DTU given was invalid (below 0 or above 1). This should be a number between 0 and 1.')
    }
    if ((outlier_gene_proportion < 0) | (outlier_gene_proportion > 1)) {
        stop('The proportion of outlier genes given was invalid (below 0 or above 1). This should be a number between 0 and 1.')
    }
    if ((proportion_of_genes_with_DGE < 0) | (proportion_of_genes_with_DGE > 1)) {
        stop('The proportion of genes with DGE given was invalid (below 0 or above 1). This should be a number between 0 and 1.')
    }
    if ((proportion_of_cells_with_inflated_counts < 0) | (proportion_of_cells_with_inflated_counts > 1)) {
        stop('The proportion of outlier cells given was invalid (below 0 or above 1). This should be a number between 0 and 1.')
    }

    # 3. counts distribution type correct?
    if (!(counts_distribution_type %in% c('negative_binomial', 'poisson'))) {
        stop("The counts distribution type given was invalid. This should either be 'poisson' or 'negative_binomial'.")
    }

    # 4. isoform_number_distribution dataframe correctly formatted? Adds up to 1? Illegal numbers not included in isospergene?
    correct_df_explanation <- paste0('This should be a dataframe with two columns:\n',
                                     '\t"isoforms_in_gene": The number of different isoforms in the gene.\n',
                                     '\t"proportion": The proportion of genes containing this number of isoforms. Proportions should add up to 1.\n',)
    if (ncol(isoform_number_distribution) != 2) {
        stop(paste0('The isoform_number_distribution dataframe has more than 2 columns. ', correct_df_explanation))
    }
    if (!(sort(colnames(isoform_number_distribution)) == c('isoforms_in_gene', 'proportion'))) {
        stop(paste0('The column names provided in the isoform_number_distribution dataframe are incorrect. ', correct_df_explanation))
    }
    if (!(is.numeric(isoform_number_distribution$isoforms_in_gene))) {
        stop(paste0('The isoforms_in_gene column of the isoform_number_distribution dataframe does not contain numeric characters.\n',
                    'Make sure that everything in the isoforms_in_gene column is an integer greater than or equal to 1.'))
    }
    if (0 %in% isoform_number_distribution$isoforms_in_gene) {
        stop(paste0('The isoforms_in_gene column of the isoform_number_distribution dataframe contains a row with a 0.\n',
                    'There cannot be genes included with 0 isoforms assigned to them.'))
    }
    if (!(is.numeric(isoform_number_distribution$proportion))) {
        stop(paste0('The proportion column of the isoform_number_distribution dataframe does not contain numeric characters.\n',
                    'Make sure that everything in the proportion column is a decimal between 0 (exclusive) and 1 (inclusive).\n'))
    }
    if (!(identical(sum(isoform_number_distribution$proportion), 1)) & !skip_distribution_dataframe_check) {
        stop(paste0('The proportions provided in the isoform_number_distribution dataframe do not add up to 1.\n',
                    'sum(isoform_number_distribution$proportions) should be 1.\n',
                    'If you believe this to be in error, you can skip this check by setting skip_distribution_dataframe_check to TRUE.'))
    }

    ### MAIN
    # Prepare the sciso_object
    sciso_object <- sim_prepare_sciso_data(num_genes = number_of_genes, cells_in_groups = cells_in_groups,
                                           gamma_shape = gamma_shape, gamma_scale = gamma_scale,
                                           isoform_distribution = isoform_number_distribution,
                                           slope_jitter = 0, show_warnings = F)

    # Swap means if DTU proportion is set here
    if (proportion_of_genes_with_DTU != 0) {

        # Repeat swaps for each group
        for (i in 1:length(cells_in_groups)) {

            group_name <- paste0('group_', i + 1)
            sciso_object <- sim_random_swap_means(sciso_object,
                                                  group_name,
                                                  proportion_to_swap = proportion_of_genes_with_DTU,
                                                  swap_from = isoforms_to_swap[1],
                                                  swap_to = isoforms_to_swap[2])
        }
    }

    # Create gene expression outliers if set here
    if (outlier_gene_proportion != 0) {
        sciso_object <- sim_make_gene_expression_outliers(sciso_object,
                                                          lognormal_location = gene_outliers_lognormal_location,
                                                          lognormal_scale = gene_outliers_lognormal_scale,
                                                          proportion_for_outlier_genes = outlier_gene_proportion)
    }

    # Create DGE between groups if set here
    if (proportion_of_genes_with_DGE != 0) {
        sciso_object <- sim_make_differential_gene_expression(sciso_object,
                                                              lognormal_location = DGE_lognormal_location,
                                                              lognormal_scale = DGE_lognormal_scale,
                                                              proportion_for_dge = proportion_of_genes_with_DGE)
    }

    # Fill the counts matrix
    sciso_object <- sim_fill_out_counts(sciso_object,
                                        runmode = counts_distribution_type,
                                        negbinom_low = negative_binomial_low_bound,
                                        negbinom_high = negative_binomial_high_bound)

    # Create cell count outliers if set
    if (proportion_of_cells_with_inflated_counts != 0) {
        sciso_object <- make_cell_outliers(sciso_object,
                                           lognormal_location = cell_outliers_lognormal_location,
                                           lognormal_scale = cell_outliers_lognormal_scale,
                                           proportion_for_outlier_cells = proportion_of_cells_with_inflated_counts)
    }

    # Add the dataset name to the sciso_object
    sciso_object$dataset_name <- dataset_name

    # If the savepath is set, save the sciso_object as a .rds file here
    if (!(is.na(save_filepath))) {
        saveRDS(sciso_object, save_filepath)
    }

    # Finish
    return(sciso_object)
}
