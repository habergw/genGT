################################################################################
#
# File: functions.R
# Author: Gregory Haber (habergw@nih.gov) and Yaakov Malinovsky
# Description: Functions for use in genGT package.
# History: Updated 05/30/19
#
################################################################################

# create_next_stage_func() creates a function which returns subsequent group
# following a positive test.
#
# This function can not be called directly.
#
# @param d: design matrix for defective set (as returned by return_hdp*).
# @param initial_design: data frame with initial design information.
# @param id_var: name of id variable in dataframe
# @param p_var: name of variable in dataframe showing prevalence
#
# @returns Function taking two arguments: index_1 and index_2. These represent
#          the ordered indices of the smallest and largest prevalences, 
#          respectively, of members in a group which has tested positive.
#          When called, this function returns a data.frame as in initial_design
#          including only individuals which should be tested togethor in the 
#          next stage.
create_next_stage_func <- function(d, initial_design, id_var, p_var) {
    force(d)
    force(initial_design)
    force(id_var)
    force(p_var)

    get_next_group <- function(index_1, index_2) {
        group_indices <- seq(index_1, index_1 + d[index_1, index_2] - 1)
        initial_design[initial_design$Order %in% group_indices, 
                       c(id_var, p_var, "Order")]
    }
}
