#include "00utilis.h"

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

double sum_product(NumericVector lhs,
                   NumericVector rhs)
{
  int n = lhs.size();
  double s = 0.0;
  for (int i= 0; i < n; i++) {
    s += lhs(i) * rhs(i);
  }
  return s;
}
void resetVector(Rcpp::NumericVector& x) {
  // Fill the vector with zeros
  std::fill(x.begin(), x.end(), 0.0);
}

