#' Generates isoform proportions within gene
#' @param isos_in_gene An integer of the number of isoforms in the gene
#' @param gene_id The name of the gene
#' @param slope_jitter The degree of random "jitter" to add/subtract to the proportion distributions
#' @noRd

sim_get_iso_props <- function(isos_in_gene, gene_id = '', slope_jitter = 0) {

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
