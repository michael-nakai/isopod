RNGversion('4.2.0')
set.seed(1)

counts_table <- readRDS(system.file('tests', 'testdata', 'test_counts_table.rds', package = 'isopod'))
desig_table <- readRDS(system.file('tests', 'testdata', 'test_designations_table.rds', package = 'isopod'))
correct_result <- readRDS(system.file('tests', 'testdata', 'test_prop_table_correct.rds', package = 'isopod'))

test_that("Function works as expected on test data", {
    result <- generate_prop_tables(counts_table, desig_table, 'transcript_id', 'gene_id', 'clusters')
    expect_identical(result, correct_result)
})


test_that("Function catches dataframe with no rows and/or no columns", {
    empty_df <- counts_table[0, ]
    expect_error(generate_prop_tables(empty_df, desig_table, 'transcript_id', 'gene_id', 'clusters'),
                 regexp = "Your transcript_counts_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    expect_error(generate_prop_tables(counts_table, empty_df, 'transcript_id', 'gene_id', 'clusters'),
                 regexp = "Your cell_labels_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    
    empty_df <- counts_table[, 0]
    expect_error(generate_prop_tables(empty_df, desig_table, 'transcript_id', 'gene_id', 'clusters'),
                 regexp = "Your transcript_counts_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    expect_error(generate_prop_tables(counts_table, empty_df, 'transcript_id', 'gene_id', 'clusters'),
                 regexp = "Your cell_labels_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    
    empty_df <- counts_table[0, 0]
    expect_error(generate_prop_tables(empty_df, desig_table, 'transcript_id', 'gene_id', 'clusters'),
                 regexp = "Your transcript_counts_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
    expect_error(generate_prop_tables(counts_table, empty_df, 'transcript_id', 'gene_id', 'clusters'),
                 regexp = "Your cell_labels_table either has no rows or columns.\nPlease fill out your counts table and rerun this function.")
})


test_that("Function catches wrong object types as dataframe", {
    empty_df <- NULL
    expect_error(generate_prop_tables(empty_df, 'transcript_id', 'gene_id', 'clusters'),
                 regexp = "Your transcript_counts_table is an empty dataframe or NULL object.\nPlease fill out your counts table and rerun this function.")
    expect_error(generate_prop_tables(c(1, 2, 3, 4, 5), 'transcript_id', 'gene_id', 'clusters'))
    expect_error(generate_prop_tables(list('1' = 1, '2' = 2, '3' = 3, '4' = 4, '5' = 5), 'transcript_id', 'gene_id', 'clusters'))
})


test_that("Function wrong object type passed as transcript_id_colname, gene_id_colname, or cell_labels_colname", {
    expect_error(generate_prop_tables(counts_table, desig_table, c(1, 2, 3, 4, 5), 'gene_id', 'clusters'))
    expect_error(generate_prop_tables(counts_table, desig_table, 'transcript_id', c(1, 2, 3, 4, 5), 'clusters'))
    expect_error(generate_prop_tables(counts_table, desig_table, 'transcript_id', 'gene_id', c(1, 2, 3, 4, 5)))
    expect_error(generate_prop_tables(counts_table, desig_table, c('hello', 'world'), 'gene_id', 'clusters'))
    expect_error(generate_prop_tables(counts_table, desig_table, 'transcript_id', c('hello', 'world'), 'clusters'))
    expect_error(generate_prop_tables(counts_table, desig_table, 'transcript_id', 'gene_id', c('hello', 'world')))
    expect_error(generate_prop_tables(counts_table, desig_table, list('1' = 1, '2' = 2), 'gene_id', 'clusters'))
    expect_error(generate_prop_tables(counts_table, desig_table, 'transcript_id', list('1' = 1, '2' = 2), 'clusters'))
    expect_error(generate_prop_tables(counts_table, desig_table, 'transcript_id', 'gene_id', list('1' = 1, '2' = 2)))
})


test_that("Function throws error when transcript_ or gene_id_colname is not present in the dataframe, same for cell_label_colname", {
    expect_error(generate_prop_tables(counts_table, desig_table, 'transcript_id_column_doesnt_exist', 'gene_id', 'clusters'))
    expect_error(generate_prop_tables(counts_table, desig_table, 'transcript_id', 'gene_id_column_doesnt_exist', 'clusters'))
    expect_error(generate_prop_tables(counts_table, desig_table, 'transcript_id', 'gene_id', 'clusters_column_doesnt_exist'))
})

