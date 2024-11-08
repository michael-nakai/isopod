#' Creates counts outliers for x% of cells and records them into the sciso_object.
#' Min 1 cell selected, NOT balanced between groups.
#' This CAN cause DTU to appear if counts in a cell are sufficiently inflated.
#' @param sciso_object The sciso_object list object created previously
#' @param lognormal_location The location parameter for the log-normal curve to pull the multiplier from
#' @param lognormal_scale The scale parameter for the log-normal curve to pull the multiplier from
#' @param proportion_for_outlier_genes The proportion of all cells that should exhibit outlier counts (between 0 and 1)
#' @return An updated sciso_object
#' @noRd

make_cell_outliers <- function(sciso_object, lognormal_location = 2.3, lognormal_scale = 0.4, proportion_for_outlier_cells = 0.005) {

    # Select cells
    cell_vec <- sciso_object$cell_designations$cell_id
    to_select <- max((length(cell_vec) * proportion_for_outlier_cells), 1)
    cells <- sample(cell_vec, to_select)

    # Create the record table
    modification_table <- data.frame(list('cell_id' = cells))

    # Pull from the log-normal distribution to get the multipliers.
    # The multipliers are rounded to whole numbers here, as we're modifying the raw counts, not means.
    modification_table$multiplier_applied <- as.integer(rlnorm(nrow(modification_table), lognormal_location, lognormal_scale))

    # Go through and apply the modification.
    temp_counts_tab <- sciso_object$counts_table
    for (i in 1:nrow(modification_table)) {
        outlier_cell <- modification_table$cell_id[i]
        multiplier <- modification_table$multiplier_applied[i]
        temp_counts_tab[[outlier_cell]] <- temp_counts_tab[[outlier_cell]] * multiplier
    }

    # Record everything to the sciso_object
    sciso_object$counts_table <- temp_counts_tab
    sciso_object$cell_outliers <- modification_table
    sciso_object$other_details[['cell_outliers_lognormal_location']] <- lognormal_location
    sciso_object$other_details[['cell_outliers_lognormal_scale']] <- lognormal_scale
    sciso_object$other_details[['cell_outliers_proportion_of_outlier_cells']] <- proportion_for_outlier_cells

    return(sciso_object)
}
