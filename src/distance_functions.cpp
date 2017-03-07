#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix euclidean_distance_rcpp (const Rcpp::NumericMatrix & x, const Rcpp::NumericMatrix & y){
  unsigned int nrow = x.nrow();
  unsigned int ncol = y.nrow();
  unsigned int i = 0, j = 0;
  double d;
  Rcpp::NumericMatrix out(nrow, ncol);
  rownames(out) = rownames(x);
  colnames(out) = rownames(y);

  for (i = 0; i < nrow - 1; i++) {
    Rcpp::NumericVector xi = x.row(i);
    for (j = i + 1; j < ncol; j++) {
      d = sqrt(sum(pow(xi - y.row(j), 2.0)));
      out(j, i) = d;
      out(i, j) = d;
    }
  }

  return out;
}
