RNGversion('4.2.0')
set.seed(1)

counts_table <- readRDS(system.file('tests', 'testdata', 'test_counts_table.rds', package = 'isopod'))

test_that("Filtering out genes with one transcript works", {
    filtered_table <- filter_counts_table(counts_table, 'transcript_id', 'gene_id')
    expect_equal(nrow(filtered_table), (nrow(counts_table) - as.integer(table(table(counts_table$gene_id))[1])))
    expect_equal(1000, nrow(filtered_table) + as.integer(table(table(counts_table$gene_id))[1]))
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
