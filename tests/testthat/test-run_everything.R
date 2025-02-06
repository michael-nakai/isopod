RNGversion('4.2.0')
set.seed(1)

dtu_counts_table <- readRDS(system.file('tests', 'testdata', 'DTU_counts_table.rds', package = 'isopod'))
dtu_designations_table <- readRDS(system.file('tests', 'testdata', 'DTU_designations_table.rds', package = 'isopod'))
test_output_folder <- system.file('tests', 'testdata', 'test_output_folder', package = 'isopod')

# This option needs to be set to bypass the user input needed to pass the warning included in
isopod_test_mode <- getOption('isopod.test_mode')
options('isopod.test_mode' = TRUE)

# Because run_everything wraps filter_counts_table() and get_permutation_pvals(),
# any arguments passed to those functions are automatically checked by them.
# Therefore, argument checking for run_everything() is only required
# for wrapper function-specific arguments.
test_that("Output folder existing throws a warning before starting.", {
    
    expect_warning(run_everything(dtu_counts_table, dtu_designations_table,
                                  'transcript_id', 'gene_id',
                                  'cluster', 'cluster_1', test_output_folder,
                                  gene_count_threshold = 20, collapse_isoforms_with_counts_below = 0,
                                  permutations = 100, cores = 1,
                                  do_gene_level_comparisons = TRUE, disable_overwrite_warning = FALSE))
    
})

# Permutations are extremely low, since we're simply testing whether it runs correctly or not
test_that("run_everything() completes successfully when all arguments are set correctly", {
    
    expect_no_error(run_everything(dtu_counts_table, dtu_designations_table,
                                  'transcript_id', 'gene_id',
                                  'cluster', 'cluster_1', test_output_folder,
                                  gene_count_threshold = 20, collapse_isoforms_with_counts_below = 0,
                                  permutations = 10, cores = 1,
                                  do_gene_level_comparisons = TRUE, disable_overwrite_warning = TRUE))
    
})

options('isopod.test_mode' = original_test_mode)
