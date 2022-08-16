#ifndef __UTILITIES__
#define __UTILITIES__

#include <Rcpp.h>
using namespace Rcpp;

NumericVector vecpow(NumericVector base, NumericVector exp);
NumericVector vecpow10(NumericVector exp);

#endif // __UTILITIES__