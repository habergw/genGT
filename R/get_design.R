################################################################################
#
# File: get_design.R
# Author: Gregory Haber (habergw@nih.gov) and Yaakov Malinovsky
# Description: Definition of get_design() function for genGT package.
# History: Updated 07/01/19
#
################################################################################

#' @title get_design
#'
#' `get_design()` returns design information for optimal algorithm for 
#' generalized hierarchical group testing screening using dynamic programming.
#'
#' @param p Vector of prevalences for population. Optional, but either p or df 
#'          must be given.
#' @param df data frame containing id variables and prevalences for all members
#'           of population.
#' @param id_var Name of variable in df containing id information. Used to
#'               name id variable in result.
#' @param p_var Name of variable in df containing prevalence information. Used
#'              to name prevalence variable in result.
#' @param Se Assay sensitivity; assumed constant across group sizes.
#' @param Sp Assay specificity; assumed constant across group sizes.
#' @param monte_carlo A logical value indicating whether Monte Carlo methods
#'                    should be used to estimate overall sensitivity and overall
#'                    specificity..
#' @param M Number of Monte Carlo iterations.
#' @param no_mc_design A logical value indicating whether misclassifcation
#'                     values should be ignored when constructing the design.
#' @return Returns a list with several design elements. 
#'   1. A data frame containing, for each member of the population, an id
#'   value (this is the same as the orginal population id if provided),
#'   prevalence value (name is p_var, defaults to p), the individuals rank
#'   in list of ordered prevalences (order), initial group assignments (group).
#'   Note, the group labels are such that all individals with the same number
#'   should be tested togethor in the first stage of the algorithm.
#'
#'  2. The expected number of tests (ET) for completing the screening algorithm.
#'  
#'  3-4. The overall sensitivity (Se) and specificity (Sp). These are provided
#'       only if monte_carlo = TRUE.
#'
#'  5. A function (get_next_group) which indicated how a group should be divided
#'     in subsequent stages following a positive group test. This function takes
#'     two arguments, index_1 and index_2, representing population rank of the
#'     individuals with the smallest and largest prevalences, respectively, in
#'     the positive group. These ranks can be obtained directly from the order
#'     column of the data frame returned with the name initial_design (the first
#'     list component).
#' @export
get_design <- function(p = NULL, df = NULL, id_var = NULL, p_var = NULL, 
                       Se = 1, Sp = 1, monte_carlo = TRUE, M = 10000,
                       no_mc_design = T) {
    # Check and assign variables
    if (is.null(p) & is.null(df))
        stop("One of p or df must be provided.")

    if (is.null(p)) {
        p_var <- as.character(substitute(p_var))
        id_var <- as.character(substitute(id_var))

        if (length(id_var) == 0 | length(p_var) == 0)
            stop("If using df, both id_var and p_var must be specified.")

        p <- eval(parse(text = p_var), envir = df)
        id <- eval(parse(text = id_var), envir = df)
    } else {
        id <- NULL
    }

    # Only allow Monte Carlo if at least one of Se or Sp < 1
    if (Se == 1 & Sp == 1)
        monte_carlo <- FALSE

    # Order by prevalence
    q <- 1 - p
    q_ordered <- sort(q, decreasing = TRUE)
    indices <- rank(p, ties.method = "first")

    # Get HDP algorithm results
    if (monte_carlo) {
        hdp <- return_hdp_mc(q_ordered, Se, Sp, M, no_mc_design)
    } else {
        hdp <- return_hdp(q_ordered, Se, Sp, no_mc_design)
    }

    # If no misclassification is used in design, and testing error is present
    # use estimated number of tests
    if (no_mc_design & (Se != 1 | Sp != 1) & monte_carlo)
        hdp$ET <- hdp$ETmc

    # Assign numeric value to missing id 
    if(is.null(id)) {
        id <- seq_along(p)
    }

    # Assign defaults for id and p variable names if unprovided
    if (is.null(id_var))
        id_var <- "id"

    if (is.null(p_var))
        p_var <- "p"

    # Create data frame to hold results
    initial_design <- data.frame(id = id, p = p, order = indices, 
                                 group = hdp$D[indices])

    names(initial_design)[1:2] <- c(id_var, p_var)

    get_next_group <- create_next_stage_func(hdp$d, initial_design, id_var, 
                                             p_var)

    design <- list(initial_design = initial_design, ET = hdp$ET)

    if (monte_carlo) {
        design$Se <- hdp$Se
        design$Sp <- hdp$Sp
    }

    design$get_next_group <- get_next_group

    class(design) <- c("hdp")

    return(design)
}

#' @title print.hdp
#'
#' Generic function for printing out design object
#'
#' @export
print.hdp <- function(x) {
    print(x$initial_design)
}
