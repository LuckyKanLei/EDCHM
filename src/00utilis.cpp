#include <Rcpp.h>

using namespace Rcpp;  

NumericVector vecpow(NumericVector a, NumericVector b) {
  NumericVector out(a.size());
  std::transform(a.begin(), a.end(),
                 b.begin(), out.begin(), ::pow);
  return out;
}
