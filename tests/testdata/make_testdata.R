# Setup
save_rds_files <- FALSE  # Keep this FALSE unless actually generating data. THIS WILL BREAK TESTS IF YOU RERUN IT.

# Helper func
# From https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r (thanks Museful!)
permutations <- function(n){
    if(n==1){
        return(matrix(1))
    } else {
        sp <- permutations(n-1)
        p <- nrow(sp)
        A <- matrix(nrow=n*p,ncol=n)
        for(i in 1:n){
            A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
        }
        return(A)
    }
}

# Creating a sample counts table (100 cells, 1000 transcripts, 300 genes)
x <- as.character(1:1000)
y <- as.character(sample.int(300, 1000, replace = TRUE))
a <- data.frame(x, y)
for (n in 1:100) {
    a[[as.character(n)]] <- as.integer(abs(rnorm(10, 1, 2)))
}
letters <- list('A', 'T', 'G', 'C', 'A', 'T', 'G', 'C')
cell_labels_mat <- matrix(letters[permutations(8)], ncol=8)
cell_labels <- apply(cell_labels_mat, 1, paste0, collapse = "")
cell_labels <- unique(cell_labels)[1:100]
colnames(a) <- c('transcript_id', 'gene_id', cell_labels)
if (save_rds_files) {
    saveRDS(a, 'tests/testdata/test_counts_table.rds')
}

# Create the designations dataframe (100 cells, 5 clusters)
x <- as.character(paste0('cluster_', sample.int(5, 100, replace = TRUE)))
a <- data.frame('clusters' = x, 'cell_barcodes' = cell_labels)
if (save_rds_files) {
    saveRDS(a, 'tests/testdata/test_designations_table.rds')
}
