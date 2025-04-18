#' Creates the counts table and designations from the sciso_object list.
#' Uses the means table, cell labels vector, groups vector, and cells_per_group vector.
#' Fills out the counts table and designations table.
#' This should be one of the last steps, applied after all outlier and modifications are made (except for cell outliers).
#' Mode can be 'poisson' or 'negative_binomial'.
#' For negative binomial, a default of between 0.2-0.7 is set to reflect values found in real data.
#' 
#' @importFrom magrittr %>%
#' @param sciso_object The list object created from `sim_prepare_sciso_data()`
#' @param runmode The curve to generate counts from. Either 'poisson' or 'negative_binomial'.
#' @param negbinom_low The low bound for the size parameter for the negative binomial curve to generate counts from. Defaults to 0.2
#' @param negbinom_high The high bound for the size parameter for the negative binomial curve to generate counts from. Defaults to 0.7
#' @return An updated sciso_object
#' @noRd

sim_fill_out_counts <- function(sciso_object, runmode = 'poisson', negbinom_low = 0.2, negbinom_high = 0.7) {

    # Check for errors in the object provided
    if (!(length(sciso_object$groups) == length(sciso_object$cells_per_group))) {
        stop(paste0("The length of the object's groups is not equal to the length of object's cells_per_group vector. ",
                    "The length of both must be equal to correctly construct a counts table."))
    }

    # Create the designations table first
    cell_group_vec <- rep(sciso_object$groups, times=sciso_object$cells_per_group)
    designations <- data.frame(list('cell_id' = sciso_object$cell_labels,
                                    'cell_group' = cell_group_vec))
    sciso_object$cell_designations <- designations

    # Initialize the counts table by giving it a transcript_id and gene_id col
    counts_table <- sciso_object$isoform_props_table %>% dplyr::select(transcript_id, gene_id)

    # If the runmode is 'negative_binomial', then randomly select group size (between 1/0.2 and 1/0.7) here:
    if (runmode == 'negative_binomial') {
        negbinom_sizes <- c(rep(NA, length(sciso_object$cells_per_group)))
        for (i in 1:length(negbinom_sizes)) {
            temp <- runif(1, negbinom_low, negbinom_high)
            negbinom_sizes[i] <- 1/temp
        }
        negbinom_recording <- data.frame('groups' = sciso_object$groups,
                                         'negative_binomial_size' = negbinom_sizes)
    }

    # Make a list of lists. The sublist names are each group, and the further sublists of the sublists are rows within that group.
    # Use .subset2 to access columns quickly here (15x faster than normal $ subsetting)
    grouplist <- list()
    meancol_names <- colnames(sciso_object$isoform_means_table)[3:ncol(sciso_object$isoform_means_table)]
    for (i in 1:length(sciso_object$groups)) {
        groupname <- sciso_object$groups[i]
        meansgroupname <- meancol_names[i]
        groupmean_vec <- .subset2(sciso_object$isoform_means_table, meansgroupname)
        cells_in_group <- sciso_object$cells_per_group[i]
        
        if (runmode != 'negative_binomial') {
            grouplist[[groupname]] <- lapply(groupmean_vec, function(x) {
                rpois(cells_in_group, x)
            })
        } else {
            negbinom_size <- negbinom_recording$negative_binomial_size[i]
            grouplist[[groupname]] <- lapply(groupmean_vec, function(x) {
                rnbinom(cells_in_group, mu = x, size = negbinom_size)
            })
        }
    }

    # For each sublist of the sublists, combine them as rows into a table
    bigtable <- as.data.frame(data.table::transpose(data.table::rbindlist(grouplist)))

    # Add the necessary data to the counts table and assign it.
    colnames(bigtable) <- sciso_object$cell_labels
    sciso_object$counts_table <- cbind(counts_table, bigtable)

    if (runmode == 'negative_binomial') {
        sciso_object$other_details[['counts_generation_negative_binomial_details']] <- negbinom_recording
        sciso_object$other_details[['negative_binomial_coeff_low']] <- negbinom_low
        sciso_object$other_details[['negative_binomial_coeff_high']] <- negbinom_high
    }

    # Return all data as a list
    return(sciso_object)
}
