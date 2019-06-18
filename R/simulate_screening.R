################################################################################
#
# File: 
# Author: Gregory Haber (habergw@nih.gov)
# Description:
# History:
#
################################################################################
#' @export
simulate_screening <- function(x, p, Se, Sp) {
    q <- 1 - p
    q_ordered <- sort(q, decreasing = TRUE)
    q_order <- order(q, decreasing = TRUE)
    
    results <- sim_screen(x[q_order], q_ordered, Se, Sp)
    
    list(results = data.frame(x = x,
                              screen_result = results$x_hat[order(q_order)],
                              p = p),
         T = results$T,
         Se = results$Se,
         Sp = results$Sp)
}
