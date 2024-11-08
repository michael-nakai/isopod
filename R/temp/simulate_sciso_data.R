library(tidyverse)
source("/oshlack_lab/michael.nakai/projects/simulate_sciso_data/generate_cell_labels.R")

# Update this whenever you CHANGE the simulation functions.
# If you're just adding separate functionality that doesn't affect existing
# functions, don't update this.
toolversion <- 1



# A helper function to get the assigned proportions for genes with x isoforms.
# Returns a vector of length x, where x is the number of isoforms in the gene.
# Jitter values are often VERY small values. The more isoforms you're calculating
# per gene, the smaller the jitter values should be. Generally, with isos_in_gene
# between 1-10, a jitter of 0.1 is OK. Anything more and the jitter should be scaled
# down to 0.01 - 0.05.
get_iso_props <- function(isos_in_gene, gene_id = '', slope_jitter = 0) {
  
  to_return <- rep(NA, isos_in_gene)
  
  m <- 1 / (-0.1617787 * isos_in_gene - 0.2111834)  # Slope calc
  b <- -log2(0.2965209 * isos_in_gene - 0.3441929)  # Intercept calc
  
  mtemp <- m + (slope_jitter * sample(c(1, -1), 1))  # Jitter addition or subtraction
  if (mtemp < 0) {  # If the slope flips to become positive, the jitter is IGNORED
    m <- mtemp
  }
  
  for (z in 1:isos_in_gene) {
    to_return[z] <- 2 ^ (m * z + b)  # General calc
  }
  
  # Often, the proportions will add up to < 1 here. The remaining proportion
  # is distributed between the two most abundant isoforms in a 3:1 ratio. This
  # avoids situations with >15 isoforms in the gene, where proportions don't
  # change much between each isoform, and reflects real data a bit better, since
  # you start seeing a better halving of counts between isos 1 and 2 than the 
  # formula makes naturally.
  # A known problem is that for genes with 3 or less isoforms in them, the total proportion
  # adds up to significantly over 1. This happens because the fitted curve for the slope calculation
  # is very accurate for values > 3, but doens't fit well for the data for genes with 2 or 3 isoforms.
  # If the proportion therefore adds up to greater than 1, the proportions for each isoform is scaled down
  # evenly. For example, a gene with two isoforms has the raw proportions of (1.0993, 0.3007). This 
  # is be scaled down to (0.785, 0.215). This also avoids problems when the jitter is set quite high, and
  # forces a steeper slope to be drawn, causing proportions to add up to greater than 1.
  if (sum(to_return) < 1) {
    to_add <- (1 - sum(to_return)) / 4
    to_return[1] <- to_return[1] + (to_add * 3)
    to_return[2] <- to_return[2] + to_add
  } else if ((sum(to_return) > 1)) {
    to_return <- to_return / sum(to_return)
  }
  
  return(to_return)
}


# Creates the start of the means table, as well as other required tables.
# Notably, this does NOT create the means for any group beyond the first.
# To create those means, either use the "random_swap_means()" function below, or fill in manually.
prepare_sciso_data <- function(num_genes, cells_in_groups = c(200, 200), 
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
      proportions <- get_iso_props(isos_in_gene, gene_id, slope_jitter)
      j <- 1
      isoform_props_table$group_1_props[i] <- proportions[j]
      previous_gene <- gene_id
    }
  }
  
  # Make the cell_labels here
  cell_labels <- generate_cell_labels(sum(cells_in_groups))
  
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


# Randomly swaps 10% of means between isoforms, and creates a new column in the isoform_means_table within the list.
# This will avoid genes with 1 transcript. It does NOT look at proportions and does NOT filter even if the proportion of 
# two isoforms is 0.5 - 0.5, since this situation is almost impossible with the updated isoform proportion generation.
# For genes with > 2 isoforms, the two most abundant isoforms are switched. All others are left untouched.
random_swap_means <- function(sciso_object, group_name, proportion_to_swap = 0.1, swap_from = 1, swap_to = 2) {
  
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


# Creates gene expression outliers by randomly selecting x% of all genes and applying a multiplier to those mean gene counts and mean isoform counts.
# Selects x% of genes without isoform switching, and x% of genes with isoform switching (balanced).
# If this function is invoked, a minimum of 1 gene with and without isoform switching will be selected, no matter how low the percentage. (unless there's literally nothing in the gorup)
# Note that with a default log-normal curve, almost all outliers will have HIGHER counts, since the multiplier will often be > 10.
make_gene_expression_outliers <- function(sciso_object, lognormal_location = 2.5, lognormal_scale = 0.3, proportion_for_outlier_genes = 0.01) {
  
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


# Simulate differential gene expression between the groups.
# First pull a DE factor from a log-normal distribution for x% of genes.
# Which group's counts were multiplied by the factor should also be recorded.
# Note that this can "stack" with the expression outlier multiplication, resulting in some genes having huge counts and huge DGE.
# Note (21/05/2024): lognormal_scale used to be 0.2, changed to 0.8 to widen the DGE expected inflator from 1.75-2.5 to 2-8 to reflect real data
make_differential_gene_expression <- function(sciso_object, lognormal_location = 0.8, lognormal_scale = 0.8, proportion_for_dge = 0.1) {
  
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


# Creates the counts table and designations from the sciso_object list.
# Uses the means table, cell labels vector, groups vector, and cells_per_group vector.
# Fills out the counts table and designations table.
# This should be one of the last steps, applied after all outlier and modifications are made (except for cell outliers).
# Mode can be 'poisson' or 'negative_binomial'
fill_out_counts <- function(sciso_object, runmode = 'poisson', negbinom_low = 0.2, negbinom_high = 0.7) {
  
  # Check for errors in the object provided
  if (!(length(sciso_object$groups) == length(sciso_object$cells_per_group))) {
    stop(paste0("The length of the object's groups is not equal to the length of object's cells_per_group vector. ",
                "The length of both must be equal to correctly construct a counts table."))
  }
  
  # Create the designations table first
  cell_group_vec <- c()
  for (i in 1:length(sciso_object$cells_per_group)) {
    cell_group_vec <- c(cell_group_vec, rep(sciso_object$groups[i], sciso_object$cells_per_group[i]))
  }
  designations <- data.frame(list('cell_id' = sciso_object$cell_labels,
                                  'cell_group' = cell_group_vec))
  sciso_object$cell_designations <- designations
  
  # Initialize the counts table by giving it a transcript_id and gene_id col
  counts_table <- sciso_object$isoform_props_table %>% select(transcript_id, gene_id)
  
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
  grouplist <- list()
  meancol_names <- colnames(sciso_object$isoform_means_table)[3:ncol(sciso_object$isoform_means_table)]
  for (i in 1:length(sciso_object$groups)) {
    groupname <- sciso_object$groups[i]
    meansgroupname <- meancol_names[i]
    temp <- list()
    groupmean_vec <- sciso_object$isoform_means_table[[meansgroupname]]
    cells_in_group <- sciso_object$cells_per_group[i]
    
    if (runmode == 'negative_binomial') {
      negbinom_size <- negbinom_recording$negative_binomial_size[i]
    }
    
    # Fill out temp with vectors corresponding to rows for that group
    # THIS DEPENDS ON MODE for how the counts are filled
    j <- 1
    for (mean_num in groupmean_vec) {
      if (runmode != 'negative_binomial') {
        rowvec <- rpois(cells_in_group, mean_num)
      } else {
        rowvec <- rnbinom(cells_in_group, mu = mean_num, size = negbinom_size)
      }
      tempname <- as.character(j)
      temp[[tempname]] <- rowvec
      j <- j + 1
    }
    grouplist[[groupname]] <- temp
  }
  
  # For each sublist of the sublists, combine them as rows into a table
  first_loop <- T
  for (groupname in names(grouplist)) {
    subtable <- as.data.frame(do.call(rbind, grouplist[[groupname]]))
    colnames(subtable) <- NULL
    rownames(subtable) <- NULL
    if (!(first_loop)) {
      bigtable <- cbind(bigtable, subtable)
      colnames(bigtable) <- NULL
      rownames(bigtable) <- NULL
    } else {
      bigtable <- subtable
      first_loop <- F
    }
  }
  
  # Add the necessary data to the counts table and assign it.
  colnames(bigtable) <- sciso_object$cell_labels
  counts_table <- cbind(counts_table, bigtable)
  sciso_object$counts_table <- counts_table
  
  if (runmode == 'negative_binomial') {
    sciso_object$other_details[['counts_generation_negative_binomial_details']] <- negbinom_recording
    sciso_object$other_details[['negative_binomial_coeff_low']] <- negbinom_low
    sciso_object$other_details[['negative_binomial_coeff_high']] <- negbinom_high
  }

  # Return all data as a list
  return(sciso_object)
}


# Creates counts outliers for x% of cells and records them into the sciso_object.
# Min 1 cell selected, NOT balanced between groups.
# This CAN cause DTU to appear if counts in a cell are sufficiently inflated.
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


# Calculates details about the table and records them into the sciso_object.
# Should run at the end.
# Should recalculate group proportion splits (and record them into a new table)
calculate_details <- function(sciso_object) {
  
  # Calculate total isoform counts
  counts_tab <- sciso_object$counts_table
  details_tab <- counts_tab %>% select(transcript_id, gene_id)
  details_tab$total_isoform_counts <- rowSums(counts_tab[, 3:ncol(counts_tab)])
  
  # Calculate total counts per group
  counts_tab_no_ids <- counts_tab[, 3:ncol(counts_tab)]
  groupsums <- t(rowsum(t(counts_tab_no_ids), sciso_object$cell_designations$cell_group[which(colnames(counts_tab_no_ids) %in% sciso_object$cell_designations$cell_id)]))
  groupsums <- as.data.frame(groupsums)
  colnames(groupsums) <- paste0('counts_in_', colnames(groupsums))
  details_tab <- cbind(details_tab, groupsums)
  
  # Calculate proportion per group (WITHIN the group, not BETWEEN the groups. So you can't reuse code from above)
  
  # Label each isoform with some info/flags, including:
  #   - Is it within a gene tagged for isoform switching?
  #   - What is the real calculated proportion split between the groups?
  
  # Is the gene an expression outlier? Basically just append the gene_expression_outliers table to the details table, and add a true/false column
  if (!(is.na(sciso_object$gene_expression_outliers))) {
    details_tab$gene_is_an_expression_outlier <- details_tab$gene_id %in% sciso_object$gene_expression_outliers$gene_id
    details_tab$temp_order <- 1:nrow(details_tab)
    details_tab <- merge(details_tab, sciso_object$gene_expression_outliers, by = 'gene_id', all.x = T)
    details_tab <- details_tab[order(details_tab$temp_order), ]
    details_tab$temp_order <- NULL
  }
  
  # Does the isoform have a cell with an expression outlier?
  # TODO: This currently only works for one outlier cell per isoform. Need to rewrite if you make it so that you can have multiple outlier cells per isoform.
  if (!(is.na(sciso_object$cell_outliers))) {
    details_tab$isoform_has_cell_expression_outlier <- details_tab$transcript_id %in% sciso_object$cell_outliers$transcript_id
    a <- sciso_object$cell_outliers
    a$gene_id <- NULL
    details_tab$temp_order <- 1:nrow(details_tab)
    details_tab <- merge(details_tab, a, by = 'transcript_id', all.x = T)
    details_tab <- details_tab[order(details_tab$temp_order), ]
    details_tab$temp_order <- NULL
  }
  
  # Is the gene tagged for isoform switching? (CODE THIS AFTER YOU FINISH IMPLEMENTING MULTIPLE ISOFORMS PER GENE)
  
  
}



