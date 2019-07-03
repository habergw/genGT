################################################################################
#
# File: simulate_screening.R
# Author: Gregory Haber (habergw@nih.gov)
# Description: Function to simulate screening process for given dataset.
# History: Updated 06/19/19
#
################################################################################
#' @title simulate_screening
#'
#' Simulates screening process for given dataset.
#'
#' @param x Vector of individual statuses (binary).
#' @param p Vector of individual prevalence estimates.
#' @param Se Assay sensitivity.
#' @param Sp Assay specificity.
#'
#' @return List given screening results, number of tests required, overall
#'   sample sensitivity, and overall sample specificity.
#' @export
simulate_screening <- function(x, p, Se, Sp, no_mc = TRUE) {
    q <- 1 - p
    q_ordered <- sort(q, decreasing = TRUE)
    q_order <- order(q, decreasing = TRUE)

    results <- sim_screen(x[q_order], q_ordered, Se, Sp, no_mc)
    
    list(results = data.frame(x = x,
                              screen_result = results$x_hat[order(q_order)],
                              p = p, groups = results$D[order(q_order)]),
         T = results$T,
         ET = results$ET,
         Se = results$Se,
         Sp = results$Sp)
}
