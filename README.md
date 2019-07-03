This package can be used to facilitate screening of large populations
using pooled tested when members of a population belong to different
risk sets for a given outcome. Specifically, it implements an optimal
ordered hierarchical group testing algorithm for this purpose. The
algorithm can incorporate assay missclassification with the assumption
that such errors are non-differential. In this document, we will explain
the basic usage of the main functions provided by the package.

Simulated data
==============

First, we simulate a dataset to be used for screening with 20
individuals. This will contain two variables

1.  **id**: unique identifier for each population member.
2.  **p**: prevalence of infection for each population member.

``` r
library(genGT)
set.seed(1)
data <- data.frame(id = 1:20, p = runif(20))

data
```

    ##    id          p
    ## 1   1 0.26550866
    ## 2   2 0.37212390
    ## 3   3 0.57285336
    ## 4   4 0.90820779
    ## 5   5 0.20168193
    ## 6   6 0.89838968
    ## 7   7 0.94467527
    ## 8   8 0.66079779
    ## 9   9 0.62911404
    ## 10 10 0.06178627
    ## 11 11 0.20597457
    ## 12 12 0.17655675
    ## 13 13 0.68702285
    ## 14 14 0.38410372
    ## 15 15 0.76984142
    ## 16 16 0.49769924
    ## 17 17 0.71761851
    ## 18 18 0.99190609
    ## 19 19 0.38003518
    ## 20 20 0.77744522

Initial study design
====================

To carry out screening, we must determine how the population members
will be grouped for the first stage of testing. To do this, we use the
**get\_design** function with our simulated data. Data can be passed to
get\_design in one of two ways: 1) as a data frame type object
containing prevalence estimates, an identifier variable, and possibly
other variables or 2) as a single vector of estimated prevalences. This
information is provided to the function using the following variables:

1.  **p**: (using method 2 above) single vector of prevalences.
2.  **df**: (using method 1 above) data.frame or similar object type
3.  **id\_var**: (using method 1 above) name of id variable in
    dataframe.
4.  **p\_var**: (using method 1 above) name of prevalence variable in
    data frame.

In addition to data information, the following variables can be passed
to **get\_design**:

1.  **Se**: Sensitivity of assay to be used for screening.
2.  **Sp**: Specificity of assay to be used for screening.
3.  **monte\_carlo** (boolean): Should Monte Carlo methods be used to
    estimate the overall testing sensitivity and specificity? (Defaults
    to TRUE.)
4.  **M**: If using Monte Carlo methods, how many iterations should be
    carried out? (Defaults to 10,000.)

For this example, we will assume Se = Sp = 0.95 and leave the other
arguments at their defaults.

We now call the function to get the initial design information.

``` r
design <- get_design(df = data, id_var = id, p_var = p, Se = 0.95, Sp = 0.95)
```

Printing the design object shows a data frame with the following
variables:

1.  **id**: This is the identifier variable passed to **id\_var**.
    (Note: the name will be the same as **id\_var** if provided.)
2.  **p**: Prevalence estimates. (Note: the name will be the same as
    **p\_var** if provided.)
3.  **order**: The rank of the individual in terms of prevalence
    estimate, from smallest to largest. This will be useful in
    determining subsequent stage groups.
4.  **group**: Initial group assignment for each individual (i.e., all
    those lablled 1 should be tested together, etc.)

Note: This same data frame can be accessed using
**design$initial\_design**.

``` r
design
```

    ##    id          p order group
    ## 1   1 0.26550866     5     2
    ## 2   2 0.37212390     6     2
    ## 3   3 0.57285336    10     5
    ## 4   4 0.90820779    18    13
    ## 5   5 0.20168193     3     1
    ## 6   6 0.89838968    17    12
    ## 7   7 0.94467527    19    14
    ## 8   8 0.66079779    12     7
    ## 9   9 0.62911404    11     6
    ## 10 10 0.06178627     1     1
    ## 11 11 0.20597457     4     1
    ## 12 12 0.17655675     2     1
    ## 13 13 0.68702285    13     8
    ## 14 14 0.38410372     8     3
    ## 15 15 0.76984142    15    10
    ## 16 16 0.49769924     9     4
    ## 17 17 0.71761851    14     9
    ## 18 18 0.99190609    20    15
    ## 19 19 0.38003518     7     3
    ## 20 20 0.77744522    16    11

Additional design information
-----------------------------

In addition to the design data frame described above, the
**get\_design** function returns additional design related variables as
follows:

1.  **ET**: Expected number of tests for carrying out screening.
2.  **Se**: Expected overall sensitivity (only returned if
    **monte\_carlo = TRUE**).
3.  **Sp**: Expected overall specificity (only returned if
    **monte\_carlo = TRUE**).

Subsequent stages
=================

Once the initial group sizes have been determined, each group can be
tested simultaneously. If a test is negative, all individuals comprising
that group are treated as negative. If positive and containing more than
one member, we must determine how to further divide the given group for
subsequent testing.

To facilitate this, the original design object is returned with a
function **get\_next\_group** which, when given the ranks of the
individuals in a group with the smallest and largest estimated
prevalences, will return a data frame listing the subgroup to be tested
in the following stage.

To illustrate this, suppose the first group, comprised of individuals
ranked 1-4, tested positive. We can then determine how this group should
be divided with the following function call.

``` r
design$get_next_group(1, 4)
```

    ##    id          p order
    ## 10 10 0.06178627     1
    ## 12 12 0.17655675     2

This output indicates that the individuals with prevalences ranked 1-2
should be retested in a new group test. The remaining initial group
members (those ranked 3 and 4) will be retested as well, with their
groupings determined by the results of the test on the members ranked 1
and 2, as described in the manuscript.
