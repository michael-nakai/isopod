#' Simulate differential gene expression between the groups.
#' First pull a DE factor from a log-normal distribution for x% of genes.
#' Which group's counts were multiplied by the factor should also be recorded.
#' Note that this can "stack" with the expression outlier multiplication, resulting in some genes having huge counts and DGE.
#' Note (21/05/2024): lognormal_scale used to be 0.2, changed to 0.8 to widen the DGE expected inflator from 1.75-2.5 to 2-8 to reflect real data
#' @param sciso_object The sciso_object list object created previously
#' @param lognormal_location The location parameter for the log-normal curve to pull the multiplier from
#' @param lognormal_scale The scale parameter for the log-normal curve to pull the multiplier from
#' @param proportion_for_dge The proportion of all genes that should exhibit DGE (between 0 and 1)
#' @return An updated sciso_object list object
#' @noRd

sim_make_differential_gene_expression <- function(sciso_object, lognormal_location = 0.8, lognormal_scale = 0.8, proportion_for_dge = 0.1) {

    # First select x percent of genes. This can select ANY gene, including genes with one isoform. A minimum of 1 gene is selected.
    number_of_genes_to_pull <- max(1, round(nrow(sciso_object$gene_means_table) * proportion_for_dge, 0))
    dge_genes <- sample(sciso_object$gene_means_table$gene_id, number_of_genes_to_pull)
    dge_genes <- dge_genes[order(dge_genes)]

    # Pull from a lognormal distribution to get the DE factor for each selected gene. The table containing dge data is also made here.
    dge_table <- data.frame(list('gene_id' = dge_genes, 'DE_factor' = rlnorm(length(dge_genes), lognormal_location, lognormal_scale)))

    # Choose which group to apply the factor to for each gene
    dge_table$group_with_multiplied_counts <- sample(sciso_object$groups, nrow(dge_table), replace = T)

    # Now multiply the counts of those gene for those groups.
    # The fill_out_counts() function actually DOESN'T look at the sciso_object$gene_means table,
    # so we can just record the updated means into subsequent columns. We also need to update the
    # mean isoform counts to reflect the updated mean gene counts (might need to recalculate here).

    # Update the gene means here
    for (i in 1:length(sciso_object$groups)) {
        g <- sciso_object$groups[i]
        updated_means <- rep(NA, nrow(sciso_object$gene_means_table))
        for (j in 1:nrow(sciso_object$gene_means_table)) {
            gene <- sciso_object$gene_means_table$gene_id[j]
            if (gene %in% dge_genes) {  # This has to be multi-step, since if the gene isn't in DGE genes, it's not in the dge_table
                if (dge_table$group_with_multiplied_counts[dge_table$gene_id == gene] == g) {
                    updated_means[j] <- sciso_object$gene_means_table$gene_means[j] * dge_table$DE_factor[dge_table$gene_id == gene]
                } else {
                    updated_means[j] <- sciso_object$gene_means_table$gene_means[j]
                }
            } else {
                updated_means[j] <- sciso_object$gene_means_table$gene_means[j]
            }
        }
        column_name <- paste0(g, '_means')
        sciso_object$gene_means_table[[column_name]] <- updated_means
    }

    # Re-allocate per-isoform counts here. Proportion split between isoforms is unchanged.
    for (g in sciso_object$groups) {
        genemean_col <- paste0(g, '_means')
        isoprop_col <- paste0(g, '_props')
        isomean_col <- paste0(g, '_means')
        for (i in 1:nrow(sciso_object$isoform_means_table)) {
            gene <- sciso_object$isoform_means_table$gene_id[i]
            sciso_object$isoform_means_table[[isomean_col]][i] <- sciso_object$gene_means_table[[genemean_col]][sciso_object$gene_means_table$gene_id == gene] *
                sciso_object$isoform_props_table[[isoprop_col]][i]
        }
    }

    sciso_object$DGE_details <- dge_table
    sciso_object$other_details[['DGE_lognormal_location']] <- lognormal_location
    sciso_object$other_details[['DGE_lognormal_scale']] <- lognormal_scale
    sciso_object$other_details[['DGE_proportion']] <- proportion_for_dge
    sciso_object$other_details[['DGE_genes']] <- dge_genes

    return(sciso_object)
}
