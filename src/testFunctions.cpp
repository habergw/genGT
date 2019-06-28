#include <Rcpp.h>
#include "ggt.h" 

using namespace Rcpp;
using namespace std;

/**
 * get_test_val() carries out binomial test for a given group with indicated
 * misclassification values. Mostly useful for simulations.
 *
 * @param ind1 (int): population index of first group member.
 * @param ind2 (int): population index of last group member.
 * @param x (vector<int>): vector of outcomes from ordered population.
 * @param Se (double): assay sensitivity.
 * @oaram Sp (double): assay specificity.
 * @param T (*int): pointer to variable storing total number of tests.
 *
 * @return (vector<int>): test outcome (0 or 1).
 */
vector<int> get_test_val(int ind1, int ind2, vector<int> x, double Se, 
        double Sp, int *T)
{
    vector<int> res(1);
    double val = 0;

    for (int i = ind1; i < ind2 + 1; i++) {
        if (x[i] == 1) {
            val = 1;
            break;
        }
    }
    
    if (val == 1) {
        res[0] = R::runif(0, 1) < Se;
    } else {
        res[0] = R::runif(0, 1) > Sp;
    }
    *T += 1;

    return res;
}

/**
 * testH() recursively carries out tests for binomial groups.
 *
 * @param ind1 (int): population index of first group member.
 * @param ind2 (int): population index of last group member.
 * @param x (vector<int>): vector of outcomes from ordered population.
 * @param Se (double): assay sensitivity.
 * @oaram Sp (double): assay specificity.
 * @param T (*int): pointer to variable storing total number of tests.
 * @param h (Rcpp::NumericMatrix): next stage group sizes for defective groups.
 *
 * @return (vector<int>) test values for all individuals in initial group.
 */
vector<int> testH(int ind1, int ind2, vector<int> x, double Se, double Sp,
                    int *T, NumericMatrix h)
{
    if (ind1 == ind2) {
      // Rcout << "testH: " << ind1 << ", " << ind2 << endl;
        return get_test_val(ind1, ind2, x, Se, Sp, T);
    } 

    vector<int> test = get_test_val(ind1, ind2, x, Se, Sp, T);

    if (test[0] == 0) {
        return vector<int>(ind2 - ind1 + 1, 0);
    } else {
        return testh(ind1, ind2, x, Se, Sp, T, h);
    }
}

/**
 * testh() recursively carries out tests for defective groups.
 *
 * @param ind1 (int): population index of first group member.
 * @param ind2 (int): population index of last group member.
 * @param x (vector<int>): vector of outcomes from ordered population.
 * @param Se (double): assay sensitivity.
 * @oaram Sp (double): assay specificity.
 * @param T (*int): pointer to variable storing total number of tests.
 * @param h (Rcpp::NumericMatrix): next stage group sizes for defective groups.
 *
 * @return (vector<int>) test values for all individuals in initial group.
 */
vector<int> testh(int ind1, int ind2, vector<int> x,  double Se, double Sp,
                    int *T, NumericMatrix h)
{
    int ind_max = ind2;           // Maximum index for this function call
    ind2 = ind1 + h(ind1, ind2) - 1;
    vector<int> res, tmp2;
    // Rcout << ind1 << ", " << ind2 << endl;
    while(1) {
        vector<int> test = get_test_val(ind1, ind2, x, Se, Sp, T);
        if (test[0] == 0) {
            vector<int> tmp(ind2 - ind1 + 1, 0);
            res.insert(res.end(), tmp.begin(), tmp.end());
            ind1 = ind2 + 1;
            if (ind2 == ind_max) {
              return res;
            }
            if (ind1 == ind_max) {
                res.push_back(1);
                return res;
            }
           ind2 = ind1 + h(ind1, ind_max) - 1;

        } else {
            if (ind1 == ind2) {
                res.push_back(1);
                tmp2 = testH(ind2 + 1, ind_max, x, Se, Sp, T, h);
                res.insert(res.end(), tmp2.begin(), tmp2.end());
                return res;
            } else {
                tmp2 = testh(ind1, ind2, x, Se, Sp, T, h); 
                res.insert(res.end(), tmp2.begin(), tmp2.end());
                tmp2 = testH(ind2 + 1, ind_max, x, Se, Sp, T, h);
                res.insert(res.end(), tmp2.begin(), tmp2.end());
            }
            break;
        }
    }
    return res;
}

/**
 * initialize_groups() populates a list with elements containing indices for
 * initial group assignments.
 *
 * @param D (Rcpp::NumericVector): vector of initial group sizes output by
 *                                 get_hdp().
 *
 * @return (list<group>) list of structs of type group containing group indices.
 */
list<group> initialize_groups(NumericVector D) {
    int N = D.size();
    list<group> groups;
    group tmp;
    int i = 0, ind1, ind2;

    while(i < N) {
        ind1 = i;
        ind2 = i + D[i] - 1;
        // Rcout << ind1 << ", " << ind2 << endl;
        tmp.ind1 = ind1;
        tmp.ind2 = ind2;
        groups.push_back(tmp);

        i += D[i];
    }

    return groups;
}

/**
 * sim_iter() carries out a single simulation iteration for calculating
 * expected number of tests, overall sensitivity, and overall specificity.
 *
 * @param q (Rcpp::NumericVector): ordered vector of individual prevalences.
 * @param initial_groups (list<group>): list of group indices generated by 
 *                                      initialize_group().
 * @param h (Rcpp::NumericMatrix): element ij is expected number of tests for
 *                                 defective group comprised of ordered 
 *                                 population members i to j.
 * @param Se (double): assay sensitivity.
 * @param Sp (double): assay specificity.
 *
 * @return (sim_vals): struct containing number of correctly classified ones
 *                     and zeros, total number of ones and zeros, and total 
 *                     number of tests.
 */
sim_vals sim_iter(NumericVector q, list<group> initial_groups, NumericMatrix h, 
                  double Se, double Sp)
{
    vector<int> res;
    sim_vals values;
    values.T = 0;
    values.est_1 = 0;
    values.est_0 = 0;
    int s = initial_groups.size();
    int N = q.size();
    vector<int> x(N);
    vector<int> tmp;

    for (int i = 0; i < N; i++) {
        x[i] = R::runif(0, 1) > q(i);
    }
    
    auto it = initial_groups.begin();

    for (int i = 0; i < s; i++) {
        tmp = testH(it->ind1, it->ind2, x, Se, Sp, &values.T, h);
        res.insert(res.end(), tmp.begin(), tmp.end());
        advance(it, 1);
    }
    
    for (int i = 0; i < N; i++) {
        if (x[i] == 1 && res[i] == 1) {
            values.est_1 += 1;
        } else if (x[i] == 0 && res[i] == 0) {
            values.est_0 += 1;
        }
    }
    values.true_1 = accumulate(x.begin(), x.end(), 0);
    values.true_0 = N - values.true_1;
    return values;
}
    
/**
 * mc_sims() performs fixed number of simulation iterations.
 *
 * @param D (Rcpp::NumericVector): optimal partition sizes for binomial case.
 * @param h (Rcpp::NumericMatrix): element ij is expected number of tests for
 *                                 defective group comprised of ordered 
 *                                 population members i to j.
 * @param q (Rcpp::NumericVector): ordered vector of individual prevalences.
 * @param Se (double): assay sensitivity.
 * @param Sp (double): assay specificity.
 * @param M (int): number of Monte Carlo iterations.
 *
 * @return
 */
mc_data mc_sims(NumericVector D, NumericMatrix h, NumericVector q, 
                double Se, double Sp, int M)
{
    sim_vals tmp;
    mc_data res;
    
    list<group> initial_groups = initialize_groups(D);
    list<sim_vals> sim_list;
    double true_0 = 0, true_1 = 0, est_0 = 0, est_1 = 0;
    double T = 0;

    for (int i = 0; i < M; i++) {
        tmp = sim_iter(q, initial_groups, h, Se, Sp);
        sim_list.push_back(tmp);
    }

    auto sim_list_iter = sim_list.begin();

    for (int i = 0; i < M; i++) {
        true_0 += sim_list_iter->true_0;
        true_1 += sim_list_iter->true_1;
        est_0 += sim_list_iter->est_0;
        est_1 += sim_list_iter->est_1;
        T += sim_list_iter->T;
        advance(sim_list_iter, 1);
    }

    res.Se_overall = est_1 / true_1;
    res.Sp_overall = est_0 / true_0;
    res.ET = T / M;
    
    return res;
}
