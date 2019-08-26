#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List model_grm_prob_rawC(NumericVector t, NumericVector a, NumericMatrix b, double D=1.702) {
  int n_p = t.size(), n_i = b.nrow(), n_c = b.ncol();
  a = a.size() == 1 ? rep_len(a, n_i) : a;
  if(a.size() != n_i)
    stop("Invalid dimension in item paraemters");
  List out(n_i);
  NumericVector ones(n_p, 1.0);
  NumericVector zeros(n_p, 0.0);
  
  for(int j = 0; j < n_i; j++) {
    NumericMatrix x(n_p, n_c + 2);
    x(_, 0) = ones;
    for(int k = 0; k < n_c; k++) {
      x(_, k + 1) = 1 / (1 + exp(-D * a[j] * (t - b(j, k))));
    }
    x(_, sum(!is_na(b(j, _))) + 1) = zeros;
    out[j] = x;
  }
  return out;
}


// [[Rcpp::export]]
List model_grm_probC(NumericVector t, NumericVector a, NumericMatrix b, double D=1.702) {
  int n_p = t.size(), n_i = b.nrow(), n_c = b.ncol();
  a = a.size() == 1 ? rep_len(a, n_i) : a;
  if(a.size() != n_i)
    stop("Invalid dimension in item paraemters");
  List out(n_i);
  NumericVector ones(n_p, 1.0);
  NumericVector zeros(n_p, 0.0);
  
  for(int j = 0; j < n_i; j++) {
    NumericMatrix x1(n_p, n_c + 1);
    NumericMatrix x0(n_p, n_c + 1);
    x1(_, 0) = ones;
    for(int k = 0; k < n_c; k++) {
      NumericVector p = 1 / (1 + exp(-D * a[j] * (t - b(j, k))));
      x1(_, k + 1) = p;
      x0(_, k) = p;
    }
    x0(_, sum(!is_na(b(j, _)))) = zeros;
    out[j] = x1 - x0;
  }
  return out;
}


    
