#include <Rcpp.h>
#include "ggt.h"

using namespace Rcpp;
using namespace std;

//' @title return_hdp_mc
//' constructs optimal hierarchical group testing design for ordered population
//' and uses Monte Carlo to estimate the expected numnber of tests, overall 
//' sensitivity, and overall specificity.
//'
//' @param q (Rcpp::NumericVector): ordered vector of individual prevalences.
//' @param Se (double): assay sensitivity.
//' @param Sp (double): assay specificity.
//' @param M  (int): number of Monte Carlo iterations
// [[Rcpp::export]]
List return_hdp_mc(NumericVector q, double Se, double Sp, int M) 
{
    List hdp = get_hdp(q, Se, Sp);

    NumericVector D = hdp[3];
    NumericMatrix d = hdp[1];
    NumericVector H = hdp[2];
    NumericVector h = hdp[0];

    mc_data vals = mc_sims(D, d, q, Se, Sp, M);

    D = get_labels(D);

    return List::create(_["D"] = D,
                        _["ET"] = H(0),
                        _["ETmc"] = vals.ET,
                        _["Se"] = vals.Se_overall,
                        _["Sp"] = vals.Sp_overall,
                        _["h"] = h,
                        _["d"] = d);
}

//' @title return_hdp
//' constructs optimal hierarchical group testing design for ordered population
//' and returns expected number of tests
//'
//' @param q (Rcpp::NumericVector): ordered vector of individual prevalences.
//' @param Se (double): assay sensitivity.
//' @param Sp (double): assay specificity.
// [[Rcpp::export]]
List return_hdp(NumericVector q, double Se, double Sp) 
{
    List hdp = get_hdp(q, Se, Sp);

    NumericVector D = hdp[3];
    NumericMatrix d = hdp[1];
    NumericVector H = hdp[2];
    NumericVector h = hdp[0];

    D = get_labels(D);
    
    return List::create(_["D"] = D, 
                        _["ET"] = H(0),
                        _["h"] = h,
                        _["d"] = d);
}

// [[Rcpp::export]]
List sim_screen(NumericVector y, NumericVector q, double Se, double Sp,
                bool no_mc)
{
  List hdp;
  if (no_mc) {  
    hdp = get_hdp(q, 1, 1);
  } else {
    hdp = get_hdp(q, Se, Sp);
  }
  
    NumericVector D = hdp[3];
    NumericVector H = hdp[2];
    NumericMatrix h = hdp[1];

    vector<int> x = as<vector<int> >(y);
    list<group> initial_groups = initialize_groups(D);
    vector<int> res;
    sim_vals values;
    values.T = 0;
    values.est_1 = 0;
    values.est_0 = 0;
    int s = initial_groups.size();
    int N = q.size();
    vector<int> tmp;
    double ESe, ESp;
    
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
 
    ESe = (double) values.est_1 / values.true_1;
    ESp = (double) values.est_0 / values.true_0;

    D = get_labels(D);

    return List::create(_["x_hat"] = res,
                        _["D"] = D,
                        _["Se"] = ESe,
                        _["Sp"] = ESp,
                        _["T"] = values.T,
                        _["ET"] = H(0));
}
