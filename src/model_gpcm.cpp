#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List model_gpcm_probC(NumericVector t, NumericVector a, NumericVector b, NumericMatrix d, double D=1.702) {
  int n_p = t.size(), n_i = b.size(), n_c = d.ncol();
  a = a.size() == 1 ? rep_len(a, n_i) : a;
  if(a.size() != n_i || b.size() != n_i || d.nrow() != n_i || d.ncol() != n_c)
    stop("Invalid dimension in item parameters");
  
  List out(n_i);
  for(int j = 0; j < n_i; j++) {
    NumericMatrix x(n_p, n_c);
    for(int i = 0; i < n_p; i++) {
      for(int k = 0; k < n_c; k++) {
        x(i, k) = D * a[j] * (t[i] - b[j] + d(j, k));
        x(i, k) += k == 0 ? 0.0 : x(i, k - 1);
      }
      x(i, _) = exp(x(i, _));
      x(i, _) = x(i, _) / sum(na_omit(x(i, _)));
    }
    out[j] = x;
  }
  return out;
}




