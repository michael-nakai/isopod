#' Generates random cell labels consisting of 10 letters.
#' Can theoretically generate up to 10^10 cells
#' @param n The number of cell labels to generate
#' @noRd

sim_generate_cell_labels <- function(n) {
  while (T) {
    final_vec <- rep(NA, n)
    for (i in 1:n) {
      v <- c(sample(LETTERS, 10, replace = TRUE))
      v <- paste0(v,collapse = "")
      final_vec[i] <- v
    }

    if (length(unique(final_vec)) == length(final_vec)) {  # An additional check to make sure all labels are unique
      break
    }
  }
  return(final_vec)
}
