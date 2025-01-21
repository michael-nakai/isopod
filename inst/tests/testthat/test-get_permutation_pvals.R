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
    expect_s3_class(pval_output$`first-loop_contingency_tables`$cluster_1, 'data.frame')
    expect_s3_class(pval_output$`first-loop_contingency_tables_gene`$cluster_1, 'data.frame')
    expect_s3_class(pval_output$unadjusted_permutation_pvalues, 'data.frame')
    expect_s3_class(pval_output$unadjusted_permutation_pvalues_gene, 'data.frame')
    expect_s3_class(pval_output$`unadjusted_first-loop_pvalues`, 'data.frame')
    expect_s3_class(pval_output$`unadjusted_first-loop_pvalues_gene`, 'data.frame')
    expect_type(pval_output$pvalue_storage_list, 'logical')
    expect_s3_class(pval_output$odds_ratio_table, 'data.frame')
    expect_type(pval_output$transcripts_filtered_from_cutoff, 'list')
    expect_s3_class(pval_output$genes_filtered_from_cutoff, 'list')
    
    # All permutation p-values should be ~1 here
    expect_equal(pval_output$permutation_pvalues$cluster_3, rep(1, nrow(pval_output$permutation_pvalues$cluster_3)))
    expect_equal(pval_output$permutation_pvalues_gene$cluster_3, rep(1, nrow(pval_output$permutation_pvalues_gene$cluster_3)))
})

test_that("Correct object types and results are returned for simple DTU dataset", {
    pval_output2 <- get_permutation_pvals(dtu_counts_table, dtu_designations_table, 'transcript_id', 'gene_id', 'cluster',
                                          analysis_group_1 = 'cluster_1', analysis_group_2 = 'cluster_2', permutations = 100,
                                          cutoff = 1)
    
})

test_that("Function catches dataframe with no rows and/or no columns", {
    empty_df <- counts_table[0, ]
    expect_error(filter_counts_table(empty_df, 'transcript_id', 'gene_id'),
                 regexp = "Your transcript_counts_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    empty_df <- counts_table[, 0]
    expect_error(filter_counts_table(empty_df, 'transcript_id', 'gene_id'),
                 regexp = "Your transcript_counts_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    empty_df <- counts_table[0, 0]
    expect_error(filter_counts_table(empty_df, 'transcript_id', 'gene_id'),
                 regexp = "Your transcript_counts_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
})


test_that("Function catches wrong object types as dataframe", {
    empty_df <- NULL
    expect_error(filter_counts_table(empty_df, 'transcript_id', 'gene_id'),
                 regexp = "Your transcript_counts_table is an empty dataframe or NULL object.\nPlease fill out your counts table and rerun this function.")
    expect_error(filter_counts_table(c(1, 2, 3, 4, 5), 'transcript_id', 'gene_id'))
    expect_error(filter_counts_table(list('1' = 1, '2' = 2, '3' = 3, '4' = 4, '5' = 5), 'transcript_id', 'gene_id'))
})


test_that("Function wrong object type passed as transcript_id_colname or gene_id_colname", {
    expect_error(filter_counts_table(counts_table, c(1, 2, 3, 4, 5), 'gene_id'))
    expect_error(filter_counts_table(counts_table, c('hello', 'world'), 'gene_id'))
    expect_error(filter_counts_table(counts_table, 'transcript_id', c(1, 2, 3, 4, 5)))
    expect_error(filter_counts_table(counts_table, 'transcript_id', c('hello', 'world')))
    expect_error(filter_counts_table(counts_table, list('1' = 1, '2' = 2), 'gene_id'))
    expect_error(filter_counts_table(counts_table, 'transcript_id', list('1' = 1, '2' = 2)))
})


test_that("Function throws error when transcript_ or gene_id_colname is not present in the dataframe", {
    expect_error(filter_counts_table(counts_table, 'transcript_id_column_doesnt_exist', 'gene_id'))
    expect_error(filter_counts_table(counts_table, 'transcript_id', 'gene_id_column_doesnt_exist'))
})
