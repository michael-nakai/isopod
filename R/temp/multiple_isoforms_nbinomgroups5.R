source("/oshlack_lab/michael.nakai/projects/simulate_sciso_data/simulate_sciso_data.R")


# Multiple isoforms dataset. Same as before but with DGE and gene expression outliers
basic_name <- 'basic_poisson'
folder_name <- 'basic_poisson5'

# Options:
num_genes = 2000
cells_in_groups = c(200, 200)
gamma_shape = 1.5
gamma_scale = 0.5
isoform_distribution = data.frame('isoforms_in_gene' = 1:20, 'proportion' = 0.05)
show_warnings = T

proportion_to_swap = 0.3

lognormal_location = 2.5
lognormal_scale = 0.3
proportion_for_outlier_genes = 0.01

proportion_for_dge = 0.4

proportion_for_outlier_cells = 0.005


# MAIN
source("/oshlack_lab/michael.nakai/projects/simulate_sciso_data/simulate_sciso_data.R")

sciso_object <- prepare_sciso_data(num_genes, cells_in_groups, gamma_shape, gamma_scale, isoform_distribution)
sciso_object <- random_swap_means(sciso_object, 'group_2', proportion_to_swap)
sciso_object <- make_gene_expression_outliers(sciso_object, lognormal_location, lognormal_scale, proportion_for_outlier_genes)
sciso_object <- make_differential_gene_expression(sciso_object, proportion_for_dge = proportion_for_dge)
sciso_object <- fill_out_counts(sciso_object)
sciso_object <- make_cell_outliers(sciso_object, proportion_for_outlier_cells = proportion_for_outlier_cells)

newfold <- paste0("/oshlack_lab/michael.nakai/projects/simulate_sciso_data/test_output/", basic_name)
dir.create(newfold, showWarnings = F)
setwd(newfold)
saveRDS(sciso_object, paste0(folder_name, '.rds'))
