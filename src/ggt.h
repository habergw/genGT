#ifndef GGT_H
#define GGT_H

#include <Rcpp.h>

// Structure to hold group information during simulations
struct group {
    int ind1, ind2;
};

// Structure to hold classification information during simulations
struct sim_vals {
    int est_0;
    int est_1;
    int true_0;
    int true_1;
    int T;
};

// Structure to hold overall simulation data
struct mc_data {
    double Se_overall, Sp_overall;
    double ET;
};

// hdpFunctions.cpp
Rcpp::List get_hdp(Rcpp::NumericVector q, double Se, double Sp);

Rcpp::NumericVector get_labels(Rcpp::NumericVector D);

// testFunctions.cpp
std::vector<int> get_test_val(int ind1, int ind2, std::vector<int> x, double Se, 
        double Sp, int *T);

std::vector<int> testH(int ind1, int ind2, std::vector<int> x, double Se, 
                       double Sp, int *T, Rcpp::NumericMatrix h);

std::vector<int> testh(int ind1, int ind2, std::vector<int> x, double Se, 
                       double Sp, int *T, Rcpp::NumericMatrix h);

std::list<group> initialize_groups(Rcpp::NumericVector D);

sim_vals sim_iter(Rcpp::NumericVector q, std::list<group> initial_groups, 
                  Rcpp::NumericMatrix h, double Se, double Sp);


mc_data mc_sims(Rcpp::NumericVector D, Rcpp::NumericMatrix h, 
                Rcpp::NumericVector q, double Se, double Sp, int M);

#endif
