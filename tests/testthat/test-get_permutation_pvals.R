RNGversion('4.2.0')
set.seed(1)

counts_table <- readRDS(system.file('tests', 'testdata', 'test_counts_table.rds', package = 'isopod'))
designations_table <- readRDS(system.file('tests', 'testdata', 'test_designations_table.rds', package = 'isopod'))
dtu_counts_table <- readRDS(system.file('tests', 'testdata', 'DTU_counts_table.rds', package = 'isopod'))
dtu_designations_table <- readRDS(system.file('tests', 'testdata', 'DTU_designations_table.rds', package = 'isopod'))

test_that("Correct object types and results are returned for null dataset", {
    pval_output <- get_permutation_pvals(counts_table, designations_table, 'transcript_id', 'gene_id', 'clusters',
                                         analysis_group_1 = 'cluster_3', analysis_group_2 = 'cluster_4', permutations = 100,
                                         cutoff = 1)
    expect_s3_class(pval_output$permutation_pvalues, 'data.frame')
    expect_s3_class(pval_output$`first-loop_pvalues`, 'data.frame')
    expect_s3_class(pval_output$permutation_pvalues_gene, 'data.frame')
    expect_s3_class(pval_output$`first-loop_pvalues_gene`, 'data.frame')
    expect_s3_class(pval_output$`first-loop_contingency_tables`$cluster_3, 'data.frame')
    expect_s3_class(pval_output$`first-loop_contingency_tables_gene`$cluster_3, 'data.frame')
    expect_s3_class(pval_output$unadjusted_permutation_pvalues, 'data.frame')
    expect_s3_class(pval_output$unadjusted_permutation_pvalues_gene, 'data.frame')
    expect_s3_class(pval_output$`unadjusted_first-loop_pvalues`, 'data.frame')
    expect_s3_class(pval_output$`unadjusted_first-loop_pvalues_gene`, 'data.frame')
    expect_type(pval_output$pvalue_storage_list, 'logical')
    expect_s3_class(pval_output$odds_ratio_table, 'data.frame')
    expect_type(pval_output$transcripts_filtered_from_cutoff, 'list')
    expect_type(pval_output$genes_filtered_from_cutoff, 'list')
    
    # All permutation p-values should be ~1 here
    expect_equal(pval_output$permutation_pvalues$cluster_3, rep(1, nrow(pval_output$permutation_pvalues)))
    expect_equal(pval_output$permutation_pvalues_gene$cluster_3, rep(1, nrow(pval_output$permutation_pvalues_gene)))
})

test_that("Correct results are returned for simple DTU dataset", {
    # Because permutation analysis inherently is unstable with low permutations (which we do ONLY for testing, since
    # more permutations = more time), we have to do relative comparisons to make sure the permutation p-values are
    # correct.
    pval_output2 <- get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                          analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                          cutoff = 1)
    
    # Transcript 3 should have the largest p-value
    expect_true(all(pval_output2$permutation_pvalues$cluster_1[3] > pval_output2$permutation_pvalues$cluster_1[c(1,2,4,5)]))
    
    # Transcript 1 and 5 should have equal or smaller p-values than transcript 2 and 4 respectively
    expect_true(pval_output2$permutation_pvalues$cluster_1[1] <= pval_output2$permutation_pvalues$cluster_1[2])
    expect_true(pval_output2$permutation_pvalues$cluster_1[5] <= pval_output2$permutation_pvalues$cluster_1[4])
})

test_that("Function catches dataframes (transcript_counts and cell_labels) with no rows and/or no columns", {
    # No rows
    empty_df <- counts_table[0, ]
    expect_error(get_permutation_pvals(empty_df, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1),
                 regexp = "Your transcript_counts_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    expect_error(get_permutation_pvals(dtu_counts_table, empty_df, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1),
                 regexp = "Your cell_labels_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    
    # No columns
    empty_df <- counts_table[, 0]
    expect_error(get_permutation_pvals(empty_df, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1),
                 regexp = "Your transcript_counts_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    expect_error(get_permutation_pvals(dtu_counts_table, empty_df, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1),
                 regexp = "Your cell_labels_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    
    # No rows or columns
    empty_df <- counts_table[0, 0]
    expect_error(get_permutation_pvals(empty_df, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1),
                 regexp = "Your transcript_counts_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    expect_error(get_permutation_pvals(dtu_counts_table, empty_df, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1),
                 regexp = "Your cell_labels_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
})


test_that("Function catches wrong object types as dataframe", {
    # NULL obj
    empty_df <- NULL
    expect_error(get_permutation_pvals(empty_df, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1),
                 regexp = "Your transcript_counts_table is an empty dataframe or NULL object.\nPlease fill out your counts table and rerun this function.")
    expect_error(get_permutation_pvals(dtu_counts_table, empty_df, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1),
                 regexp = "Your cell_labels_table is an empty dataframe or NULL object.\nPlease fill out your counts table and rerun this function.")
    
    # Vectors
    expect_error(get_permutation_pvals(c(1, 2, 3, 4, 5), dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, c(1, 2, 3, 4, 5), 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    
    # Raw list obj
    expect_error(get_permutation_pvals(list('1' = 1, '2' = 2, '3' = 3, '4' = 4, '5' = 5), dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, list('1' = 1, '2' = 2, '3' = 3, '4' = 4, '5' = 5), 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
})


test_that("Function wrong object type passed as arguments that should be characters", {
    # transcript_id_colname
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, c(1, 2, 3, 4, 5), 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, c('hello', 'world'), 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, list('1' = 1, '2' = 2), 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    
    # gene_id_colname
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', c(1, 2, 3, 4, 5), 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', c('hello', 'world'), 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', list('1' = 1, '2' = 2), 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    
    # cell_labels_colname
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', c(1, 2, 3, 4, 5),
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', c('hello', 'world'),
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', list('1' = 1, '2' = 2),
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    
    # analysis_group_1
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = c(1, 2, 3, 4, 5), analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = c('hello', 'world'), analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = list('1' = 1, '2' = 2), analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    
    # analysis_group_2
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_2 = c(1, 2, 3, 4, 5), analysis_group_1 = 'cluster_1', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_2 = c('hello', 'world'), analysis_group_1 = 'cluster_1', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_2 = list('1' = 1, '2' = 2), analysis_group_1 = 'cluster_1', permutations = 100,
                                       cutoff = 1))
})


test_that("Function throws error when provided colnames are not present in the provided dataframes", {
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id_doesnt_exist', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id_doesnt_exist', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster_doesnt_exist',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1_not_present', analysis_group_2 = 'cluster_2', permutations = 100,
                                       cutoff = 1))
    expect_error(get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                       analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2_not_present', permutations = 100,
                                       cutoff = 1))
})
