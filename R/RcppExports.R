# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title return_hdp_mc
#' constructs optimal hierarchical group testing design for ordered population
#' and uses Monte Carlo to estimate the expected numnber of tests, overall 
#' sensitivity, and overall specificity.
#'
#' @param q (Rcpp::NumericVector): ordered vector of individual prevalences.
#' @param Se (double): assay sensitivity.
#' @param Sp (double): assay specificity.
#' @param M  (int): number of Monte Carlo iterations
return_hdp_mc <- function(q, Se, Sp, M, no_mc_design) {
    .Call('_genGT_return_hdp_mc', PACKAGE = 'genGT', q, Se, Sp, M, no_mc_design)
}

#' @title return_hdp
#' constructs optimal hierarchical group testing design for ordered population
#' and returns expected number of tests
#'
#' @param q (Rcpp::NumericVector): ordered vector of individual prevalences.
#' @param Se (double): assay sensitivity.
#' @param Sp (double): assay specificity.
return_hdp <- function(q, Se, Sp, no_mc_design) {
    .Call('_genGT_return_hdp', PACKAGE = 'genGT', q, Se, Sp, no_mc_design)
}

sim_screen <- function(y, q, Se, Sp, no_mc) {
    .Call('_genGT_sim_screen', PACKAGE = 'genGT', y, q, Se, Sp, no_mc)
}

