#include <Rcpp.h>
#include "00utilis.h"
using namespace Rcpp;
// [[Rcpp::interfaces(r, cpp)]]



//' **baseflow**
//' @name baseflow
//' @param ground_water_mm (mm/m2) water volum in `groundLy`
//' @param ground_capacity_mm (mm/m2) water storage capacity in `soilLy` or interceptof `landLy`
//' @param param_baseflow_sup_k,param_baseflow_sup_gamma parameters for [baseflow_SupplyPow()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_SupplyPow(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector param_baseflow_sup_k,
    NumericVector param_baseflow_sup_gamma
)
{
  NumericVector baseflow_, k_;
  
  k_ = param_baseflow_sup_k * vecpow((ground_water_mm / ground_capacity_mm), param_baseflow_sup_gamma);
  baseflow_ = k_ * ground_water_mm;
  
  return ifelse(baseflow_ > ground_water_mm, ground_water_mm, baseflow_) ;
}

//' @rdname baseflow
//' @param param_baseflow_sur_k parameters for [baseflow_SupplyRatio()]
//' @return baseflow (mm/m2) 
//' @export
// [[Rcpp::export]]
NumericVector baseflow_SupplyRatio(
    NumericVector ground_water_mm,
    NumericVector param_baseflow_sur_k
)
{
  
  return param_baseflow_sur_k * ground_water_mm;
  
}
