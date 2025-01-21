#' Simulate a transcript-level single-cell counts table
#'
#' This function wraps all other simulation functions to make a user-friendly one-liner that
#' creates a transcript-level single-cell dataset, possibly exhibiting DTU between cell groups.
#'
#' @importFrom checkmate checkPathForOutput
#' @param dataset_name The name of the dataset generated.
#' @param number_of_genes The total number of genes to generate.
#' @param cells_in_groups A vector defining the number of cell groups, as well as the number of cells in each group.
#' For example, c(100, 200) would create two cell groups of 100 and 200 cells respectively.
#' @param isoform_number_distribution A dataframe with two columns. The 'isoforms_in_gene" column contains the number of isoforms per gene,
#' and the 'proportion' column contains a float between 0-1 that determines the proportion of genes to have that isoform number.
#' @param gamma_shape The shape parameter of the gamma distribution to pull gene means from.
#' @param gamma_scale The scale parameter of the gamma distribution to pull gene means from.
#' @param proportion_of_genes_with_DTU A decimal between 0 and 1 that defines the proportion of all genes that should exhibit DTU.
#' @param isoforms_to_swap A two-element vector that defines which isoforms should have their proportions swapped to simulate DTU.
#' Isoforms are numbered by abundance, where isoform 1 contains the highest proportion of gene counts.
#' @param outlier_gene_proportion The proportion of all genes that should exhibit outlier counts (between 0 and 1).
#' @param gene_outliers_lognormal_location The location parameter for the log-normal curve to pull the gene outlier counts multiplier from.
#' @param gene_outliers_lognormal_scale The scale parameter for the log-normal curve to pull the gene outlier counts multiplier from.
#' @param proportion_of_genes_with_DGE The proportion of all genes that should exhibit DGE (between 0 and 1).
#' @param DGE_lognormal_location The location parameter for the log-normal curve to pull the DGE count multiplier from.
#' @param DGE_lognormal_scale The scale parameter for the log-normal curve to pull the DGE count multiplier from.
#' @param counts_distribution_type The curve to generate counts from. Either 'poisson' or 'negative_binomial'.
#' @param negative_binomial_low_bound The low bound for the size parameter for the negative binomial curve to generate counts from. Defaults to 0.2.
#' @param negative_binomial_high_bound The high bound for the size parameter for the negative binomial curve to generate counts from. Defaults to 0.7.
#' @param proportion_of_cells_with_inflated_counts The proportion of all cells that should exhibit outlier counts (between 0 and 1).
#' @param cell_outliers_lognormal_location The location parameter for the log-normal curve to pull the cell outlier counts multiplier from.
#' @param cell_outliers_lognormal_scale The scale parameter for the log-normal curve to pull the cell outlier counts multiplier from.
#' @param save_filepath If set, saves the simulated data as an RDS object at the filepath provided. Defaults to NA.
#' @param skip_distribution_dataframe_check Skips error checking of the dataframe provided in the isoform_number_distribution dataframe. Defaults to FALSE.
#' @return A list object containing all simulation data and parameters.
#' @export

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
                                     '\t"proportion": The proportion of genes containing this number of isoforms. Proportions should add up to 1.\n')
    if (ncol(isoform_number_distribution) != 2) {
        stop(paste0('The isoform_number_distribution dataframe has more than 2 columns. ', correct_df_explanation))
    }
    if (!sum(sort(colnames(isoform_number_distribution)) == c('isoforms_in_gene', 'proportion')) == 2) {
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
    # Repeat swaps for each group
    for (i in 2:length(cells_in_groups)) {

        group_name <- paste0('group_', i)
        sciso_object <- sim_random_swap_means(sciso_object,
                                              group_name,
                                              proportion_to_swap = proportion_of_genes_with_DTU,
                                              swap_from = isoforms_to_swap[1],
                                              swap_to = isoforms_to_swap[2])
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
