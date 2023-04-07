#ifndef __UTILITIES__
#define __UTILITIES__

#include <Rcpp.h>
using namespace Rcpp;

NumericVector vecpow(NumericVector base, NumericVector exp);
NumericVector vecpow10(NumericVector exp);
double sum_product(NumericVector lhs, NumericVector rhs);
void resetVector(Rcpp::NumericVector& x);
#endif // __UTILITIES__