#include <Rcpp.h>
using namespace Rcpp;


//' @rdname model_3pl
//' @param t ability parameters, 1d vector
//' @param a discrimination parameters, 1d vector
//' @param b difficulty parameters, 1d vector
//' @param c guessing parameters, 1d vector
//' @param D the scaling constant, default=1.702
//' @return \code{model_3pl_prob} returns the resulting probabilities in a matrix
//' @examples
//' with(model_3pl_gendata(10, 5), model_3pl_prob(t, a, b, c))
//' @export
// [[Rcpp::export]]
NumericMatrix model_3pl_prob(NumericVector t, NumericVector a, NumericVector b, NumericVector c, double D=1.702) {
  int n_p = t.size(), n_i=b.size();
  a = a.size() == 1 ? rep_len(a, n_i) : a;
  c = c.size() == 1 ? rep_len(c, n_i) : c;
  if(a.size() != n_i || c.size() != n_i)
    stop("Invalid dimension in item parameters");
  
  NumericMatrix out(n_p, n_i);
  for(int j = 0; j < n_i; j++) {
    out(_, j) = c[j] + (1 - c[j]) / (1 + exp(-D * a[j] * (t - b[j])));
  }
  
  return out;
}


//' @rdname model_3pl
//' @return \code{model_3pl_info} returns the resulting information in a matrix
//' @examples
//' with(model_3pl_gendata(10, 5), model_3pl_info(t, a, b, c))
//' @export
// [[Rcpp::export]]
NumericMatrix model_3pl_info(NumericVector t, NumericVector a, NumericVector b, NumericVector c, double D=1.702) {
  int n_p = t.size(), n_i=b.size();
  a = a.size() == 1 ? rep_len(a, n_i) : a;
  c = c.size() == 1 ? rep_len(c, n_i) : c;
  NumericMatrix p = model_3pl_prob(t, a, b, c, D);
  NumericMatrix out(n_p, n_i);
  
  for(int j = 0; j < n_i; j++) {
    out(_, j) = pow(D * a[j] * (p(_, j) - c[j]) / (1 - c[j]), 2) * (1 - p(_, j)) / p(_, j);
  }
  
  return out;
}


// [[Rcpp::export]]
List model_3pl_dvC(NumericMatrix u, List quad, NumericVector a, NumericVector b, NumericVector c, double D=1.702) {
  NumericVector quad_t = quad["t"], quad_w = quad["w"];
  NumericMatrix p = model_3pl_prob(quad_t, a, b, c, D);
  int n_p = u.nrow(), n_i = u.ncol(), n_q = quad_t.size();
  
  List p0(n_q);
  NumericMatrix p1(n_p, n_q);
  NumericVector p2(n_p);
  for(int k = 0; k < n_q; k++) {
    NumericMatrix lh(n_p, n_i);
    for(int i = 0; i < n_p; i++) {
      p1(i, k) = 1.0;
      for(int j = 0; j < n_i; j++) {
        lh(i, j) = pow(p(k, j), u(i, j)) * pow(1 - p(k, j), 1 - u(i, j));
        p1(i, k) *= NumericVector::is_na(lh(i, j)) ? 1.0 : lh(i, j);
      }
    }
    p0[k] = lh;
    p2 += p1(_, k) * quad_w[k];
  }
  
  return List::create(p0, p1, p2);          
}