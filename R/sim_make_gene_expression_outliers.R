#' Creates gene expression outliers by randomly selecting x% of all genes and applying a multiplier to those mean gene counts and mean isoform counts.
#' Selects x% of genes without isoform switching, and x% of genes with isoform switching (balanced).
#' If this function is invoked, a minimum of 1 gene with and without isoform switching will be selected, no matter how low the percentage. (unless there's literally nothing in the gorup)
#' Note that with a default log-normal curve, almost all outliers will have HIGHER counts, since the multiplier will often be > 10.
#' @param sciso_object The sciso_object list object created previously
#' @param lognormal_location The location parameter for the log-normal curve to pull the multiplier from
#' @param lognormal_scale The scale parameter for the log-normal curve to pull the multiplier from
#' @param proportion_for_outlier_genes The proportion of all genes that should exhibit outlier counts (between 0 and 1)
#' @return An updated sciso_object
#' @noRd

sim_make_gene_expression_outliers <- function(sciso_object, lognormal_location = 2.5, lognormal_scale = 0.3, proportion_for_outlier_genes = 0.01) {

    gene_expression_outliers <- list()
    for (x in names(sciso_object$swapped_genes)) {
        number_of_genes_in_swapped_group <- length(sciso_object$swapped_genes[[x]])

        # Guarantees a minimum of 1 gene per group is chosen
        number_to_choose <- ifelse((number_of_genes_in_swapped_group * proportion_for_outlier_genes) >= 1,
                                   as.integer(number_of_genes_in_swapped_group * proportion_for_outlier_genes),
                                   1)
        label_name <- paste0('genes_with_isoform_switching_in_', x)

        if (number_of_genes_in_swapped_group == 0) {
            next
        } else {
            gene_expression_outliers[[label_name]] <- sample(sciso_object$swapped_genes[[x]], number_to_choose)
        }
    }

    # We also need to nominate x% of genes from genes WITHOUT isoform switching
    swapped_genes_vector <- Reduce(c, sciso_object$swapped_genes)
    unswapped_genes <- sciso_object$gene_means_table$gene_id[!(sciso_object$gene_means_table$gene_id %in% swapped_genes_vector)]
    number_to_choose <- ifelse((length(unswapped_genes) * proportion_for_outlier_genes) >= 1,
                               as.integer(length(unswapped_genes) * proportion_for_outlier_genes),
                               1)
    gene_expression_outliers[['genes_without_isoform_switching']] <- sample(unswapped_genes, number_to_choose)

    # Make a table of the genes that will receive outlier count multiplication.
    # Has the columns: gene_id, count_multiplier
    outlier_gene_table <- data.frame(list('gene_id' = Reduce(c, gene_expression_outliers)))
    outlier_gene_table$count_multiplier <- rlnorm(nrow(outlier_gene_table), lognormal_location, lognormal_scale)
    sciso_object$gene_expression_outliers <- outlier_gene_table

    # Apply the multipliers for all mean gene and isoform counts for the selected genes
    for (i in 1:nrow(outlier_gene_table)) {
        gene <- outlier_gene_table$gene_id[i]
        multiplier <- outlier_gene_table$count_multiplier[i]
        sciso_object$gene_means_table$gene_means[sciso_object$gene_means_table$gene_id == gene] <-
            sciso_object$gene_means_table$gene_means[sciso_object$gene_means_table$gene_id == gene] * multiplier

        for (column_name in colnames(sciso_object$isoform_means_table)[3:ncol(sciso_object$isoform_means_table)]) {
            sciso_object$isoform_means_table[[column_name]][which(sciso_object$isoform_means_table$gene_id == gene)] <-
                sciso_object$isoform_means_table[[column_name]][which(sciso_object$isoform_means_table$gene_id == gene)] * multiplier
        }
    }

    sciso_object$other_details[['gene_expression_outlier_lognormal_location']] <- lognormal_location
    sciso_object$other_details[['gene_expression_outlier_lognormal_scale']] <- lognormal_scale
    sciso_object$other_details[['gene_expression_outlier_proportion']] <- proportion_for_outlier_genes

    return(sciso_object)

}
