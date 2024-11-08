#' Creates the start of the means table, as well as other required tables.
#' Notably, this does NOT create the means for any group beyond the first.
#' To create those means, either use the "random_swap_means()" function below, or fill in manually.
#' @param num_genes The number of genes to generate
#' @param cells_in_groups The number of cells in each group to generate, as a vector. Defaults to 200 cells between 2 groups
#' @param gamma_shape The shape parameter of the gamma distribution to pull gene means from
#' @param gamma_scale The scale parameter of the gamma distribution to pull gene means from
#' @param isoform_distribution A dataframe with two columns. The 'isoforms_in_gene" column contains the number of isoforms per gene,
#' and the 'proportion' column contains a float between 0-1 that determines the proportion of genes to have that isoform number.
#' @param slope_jitter See the 'slope_jitter' argument in `sim_get_iso_proportions()`
#' @param show_warnings Shows/suppresses warnings regarding the isoform_distribution dataframe if manually provided.
#' @return A list of all generated data, with blank openings for data generated in subsequent functions.
#' @noRd

sim_prepare_sciso_data <- function(num_genes, cells_in_groups = c(200, 200),
                                   gamma_shape = 1.5, gamma_scale = 0.5,
                                   isoform_distribution = data.frame('isoforms_in_gene' = 1:10, 'proportion' = rep(0.1, 10)),
                                   slope_jitter = 0, show_warnings = T) {

    # Draw means from gamma distribution, which needs a shape (k) and a scale (theta) parameter
    # A distribution of k = 1.5, B = 0.5 centers the distribution around ~0.25, with 26% of genes with a mean >1
    gene_means <- rgamma(num_genes, shape = gamma_shape, scale = gamma_scale)
    gene_means_dataframe <- data.frame(list('gene_id' = paste0('gene_', 1:num_genes), 'gene_means' = gene_means))

    # Check that the proportions in isoform_distribution add up to 1.
    # If the sum is < 1, add proportion to the first row to fix and notify (so there'll be more genes with 1 isoform, which should be filtered out anyway).
    # Otherwise, error out.
    if (!(identical(sum(isoform_distribution$proportion), 1))) {
        if (sum(isoform_distribution$proportion) < 1) {
            isoform_distribution$proportion[1] <- isoform_distribution$proportion[1] + (1 - sum(isoform_distribution$proportion))
            if (show_warnings) {
                warning(paste0('Isoform proportions in the provided isoform distribution dataframe did not add up to 1.\n',
                               round((1 - sum(isoform_distribution$proportion)), 4), ' has been added to the first row of the distribtuon table.'))
            }
        } else {
            stop('Isoform proportions in the provided isoform distribution dataframe added up to > 1. Please fix the proportions column and rerun.')
        }
    }

    # Determine number of genes with various numbers of isoforms
    isoform_distribution$gene_numbers <- as.integer(num_genes * isoform_distribution$proportion)
    if (!(identical(as.integer(sum(isoform_distribution$gene_numbers)), as.integer(sum(num_genes))))) {
        isoform_distribution$gene_numbers[1] <- isoform_distribution$gene_numbers[1] + (sum(num_genes) - sum(isoform_distribution$gene_numbers))
    }

    # Pre-assign and create the isoform proportions table
    gene_names <- c()
    prev_end <- 0
    for (i in 1:nrow(isoform_distribution)) {
        isos_in_gene <- isoform_distribution$isoforms_in_gene[i]
        total_genes_to_assign_props <- isoform_distribution$gene_numbers[i]
        current_start <- prev_end + 1
        current_end <- prev_end + total_genes_to_assign_props
        tempgenes <- rep(paste0('gene_', current_start:current_end), each = isos_in_gene)
        gene_names <- c(gene_names, tempgenes)
        prev_end <- current_end
    }

    iso_names <- rep(NA, length(gene_names))
    previous_gene <- gene_names[1]
    j <- 0
    for (i in 1:length(iso_names)) {
        current_gene <- gene_names[i]
        if (current_gene == previous_gene) {
            j <- j + 1
            iso_names[i] <- paste0(current_gene, '_isoform_', j)
        } else {
            j <- 1
            previous_gene <- current_gene
            iso_names[i] <- paste0(current_gene, '_isoform_', j)
        }
    }

    isoform_props_table <- data.frame(list('gene_id' = gene_names, 'transcript_id' = iso_names, 'group_1_props' = rep(NA, length(iso_names))))

    # Loop over the different isoforms_in_gene genes, generate x of each and assign proportions.
    freqs <- table(isoform_props_table$gene_id)
    previous_gene <- 'a'
    j <- 0
    proportions <- NA
    for (i in 1:nrow(isoform_props_table)) {
        gene_id <- isoform_props_table$gene_id[i]
        isos_in_gene <- unname(freqs[gene_id])

        # If there's only 1 isoform in the gene, set the prop to 1 and go on.
        if (isos_in_gene == 1) {
            isoform_props_table$group_1_props[i] <- 1
            previous_gene <- gene_id
            next
        }

        # Otherwise, check to see if we've already generated proportions for this gene by comparing
        # the current gene_id to the previous row's gene_id. If so, keep taking from those proportions.
        # If not, generate new proportions and start taking from them. Importantly, any gene that has
        # 1 isoform in the gene should never reach this step because of the loop control in the previous block.
        if (gene_id == previous_gene) {
            j <- j + 1
            isoform_props_table$group_1_props[i] <- proportions[j]
        } else {
            proportions <- sim_get_iso_props(isos_in_gene, gene_id, slope_jitter)
            j <- 1
            isoform_props_table$group_1_props[i] <- proportions[j]
            previous_gene <- gene_id
        }
    }

    # Make the cell_labels here
    cell_labels <- sim_generate_cell_labels(sum(cells_in_groups))

    # Generate the isoform-level means here
    iso_means <- rep(NA, nrow(isoform_props_table))
    for (i in 1:nrow(isoform_props_table)) {
        gene_id <- isoform_props_table$gene_id[i]
        isoprop <- isoform_props_table$group_1_props[i]
        gene_avg <- gene_means_dataframe$gene_means[gene_means_dataframe$gene_id == gene_id]
        iso_means[i] <- gene_avg * isoprop
    }

    # Make a table of isoform-level proportions and means.
    isoform_means_table <- data.frame(list('gene_id' = gene_names, 'transcript_id' = iso_names, 'group_1_means' = iso_means))

    # Return an list of all generated data
    return_list <- list('gene_means_table' = gene_means_dataframe,
                        'cell_designations' = NA,
                        'counts_table' = NA,
                        'isoform_props_table' = isoform_props_table,
                        'isoform_means_table' = isoform_means_table,
                        'cells_per_group' = cells_in_groups,
                        'cell_labels' = cell_labels,
                        'gene_expression_outliers' = NA,
                        'cell_outliers' = NA,
                        'DGE_details' = NA,
                        'groups' = c('group_1'),
                        'other_details' = list('gamma_shape' = gamma_shape,
                                               'gamma_scale' = gamma_scale,
                                               'slope_jitter' = slope_jitter,
                                               'isoform_distribution' = isoform_distribution,
                                               'swapped_isos' = NA,
                                               'tool_version' = toolversion))

    return(return_list)
}
