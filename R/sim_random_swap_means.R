#' Randomly swaps 10% of means between isoforms, and creates a new column in the isoform_means_table within the list.
#' This will avoid genes with 1 transcript. It does NOT look at proportions and does NOT filter even if the proportion of
#' two isoforms is 0.5 - 0.5, since this situation is almost impossible with the updated isoform proportion generation.
#' For genes with > 2 isoforms, the two most abundant isoforms are switched. All others are left untouched.
#' @param sciso_object The list object created from `sim_prepare_sciso_data()`
#' @param group_name The name of the group to generate means for
#' @param proportion_to_swap The proportion of genes that should exhibit DTU
#' @param swap_from The isoform number to swap means with, with 1 being the most abundant isoform
#' @param swap_to The isoform number to swap means with, with 1 being the most abundant isoform
#' @return An updated sciso_object list object
#' @noRd

sim_random_swap_means <- function(sciso_object, group_name, proportion_to_swap = 0.1, swap_from = 1, swap_to = 2) {

    # Get a vector of gene ids without the genes with one isoform
    genes_with_at_least_x_isoforms <- rep(NA, nrow(sciso_object$gene_means_table))
    freqs <- table(sciso_object$isoform_props_table$gene_id)
    for (i in 1:nrow(sciso_object$gene_means_table)) {
        gene_id <- sciso_object$gene_means_table$gene_id[i]
        if (freqs[gene_id] >= swap_to) {
            genes_with_at_least_x_isoforms[i] <- gene_id
        }
    }
    genes_with_at_least_x_isoforms <- genes_with_at_least_x_isoforms[!is.na(genes_with_at_least_x_isoforms)]

    # Figure out which to swap
    # Swaps total number of genes * proportion, OR all genes with at least x isoforms if total*prop is more than the eligible genes
    genes_to_swap <- sample(genes_with_at_least_x_isoforms,
                            min(as.integer(length(sciso_object$gene_means_table$gene_id) * proportion_to_swap),
                                length(genes_with_at_least_x_isoforms)))

    # Store the vector of swapped genes into the swapped genes list in the original object
    if (!('swapped_genes' %in% names(sciso_object))) {
        sciso_object$swapped_genes <- list()
    }
    sciso_object$swapped_genes[[group_name]] <- genes_to_swap

    # Do the swaps
    # First split the props_table by isoform, then reverse the genes listed previously
    # Then unsplit by the props_table gene_id col, and reinsert as a new column into the props_table
    props_table <- sciso_object$isoform_props_table
    split_by_iso <- split(props_table$group_1_props, props_table$gene_id)
    for (gene_name in names(split_by_iso)) {
        if (gene_name %in% genes_to_swap) {
            storage_to_swap <- split_by_iso[[gene_name]][swap_from]
            split_by_iso[[gene_name]][swap_from] <- split_by_iso[[gene_name]][swap_to]
            split_by_iso[[gene_name]][swap_to] <- storage_to_swap
        }
    }
    props_col_name <- paste0(group_name, '_props')
    sciso_object$isoform_props_table[[props_col_name]] <- unsplit(split_by_iso, props_table$gene_id)

    # Do the same for the means
    means_table <- sciso_object$isoform_means_table
    split_by_iso <- split(means_table$group_1_means, means_table$gene_id)
    for (gene_name in names(split_by_iso)) {
        if (gene_name %in% genes_to_swap) {
            storage_to_swap <- split_by_iso[[gene_name]][swap_from]
            split_by_iso[[gene_name]][swap_from] <- split_by_iso[[gene_name]][swap_to]
            split_by_iso[[gene_name]][swap_to] <- storage_to_swap
        }
    }
    means_col_name <- paste0(group_name, '_means')
    sciso_object$isoform_means_table[[means_col_name]] <- unsplit(split_by_iso, means_table$gene_id)

    # Add group_name to the group vector
    sciso_object$groups <- c(sciso_object$groups, group_name)
    sciso_object$other_details$swapped_isos <- c(swap_from, swap_to)

    # Return the new full object
    return(sciso_object)
}
