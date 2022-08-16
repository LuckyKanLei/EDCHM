#include <Rcpp.h>

using namespace Rcpp;  

NumericVector vecpow(NumericVector base, NumericVector exp) {
  NumericVector out(base.size());
  std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), ::pow);
  return out;
}


NumericVector vecpow10(NumericVector exp) {
  
  NumericVector base10 (exp.size(), 10.0);
  
  return vecpow(base10, exp);
}
