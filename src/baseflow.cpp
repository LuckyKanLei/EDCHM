#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]



//' **baseflow**
//' @name baseflow
//' @inheritParams all_vari
//' @return baseflow (mm/m2) 
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
//' @export
// [[Rcpp::export]]
NumericVector baseflow_SupplyRatio(
    NumericVector ground_water_mm,
    NumericVector param_baseflow_sur_k
)
{
  
  return param_baseflow_sur_k * ground_water_mm;
  
}


//' @rdname baseflow
//' @export
// [[Rcpp::export]]
NumericVector baseflow_GR4J(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm
)
{
  NumericVector baseflow_, k_;
  
  k_ = 1 - pow((1 + pow(ground_water_mm / ground_capacity_mm, 4)), -0.25);
  baseflow_ = k_ * ground_water_mm;
  
  return ifelse(baseflow_ > ground_water_mm, ground_water_mm, baseflow_) ;
}

//' @rdname baseflow
//' @param param_baseflow_grf_gamma parameters for [baseflow_GR4Jfix()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_GR4Jfix(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector param_baseflow_grf_gamma
)
{
  NumericVector baseflow_, k_;
  
  k_ = 1 - vecpow((1 + vecpow(ground_water_mm / ground_capacity_mm, param_baseflow_grf_gamma)), -1.0 / param_baseflow_grf_gamma);
  baseflow_ = k_ * ground_water_mm;
  
  return ifelse(baseflow_ > ground_water_mm, ground_water_mm, baseflow_) ;
}
