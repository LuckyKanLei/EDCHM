#include <Rcpp.h>
#include "00utilis.h"
using namespace Rcpp;
// [[Rcpp::interfaces(r, cpp)]]



//' **lateral flux**
//' @name lateral
//' @inheritParams all_vari
//' @param param_lateral_gr4_k parameters
//' @return lateral_mm (mm/m2)
//' @export
// [[Rcpp::export]]
NumericVector lateral_GR4J(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector param_lateral_gr4_k
) 
{
  NumericVector ground_lateral_mm, k_;
  NumericVector ground_diff_mm = (ground_capacity_mm - ground_water_mm);
  k_ = param_lateral_gr4_k * pow((ground_water_mm / ground_capacity_mm), 3.5);
  ground_lateral_mm = k_ * ground_water_mm;
  ground_lateral_mm = ifelse(ground_lateral_mm > ground_diff_mm, ground_diff_mm, ground_lateral_mm) ;
  return ifelse(ground_lateral_mm > - ground_water_mm, ground_lateral_mm, - ground_water_mm) ;
}

//' @rdname lateral
//' @param param_lateral_sur_k parameters for [lateral_SupplyRatio()]
//' @export
// [[Rcpp::export]]
NumericVector lateral_SupplyRatio(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector param_lateral_sur_k
)
{
  
  NumericVector ground_diff_mm = (ground_capacity_mm - ground_water_mm);
  NumericVector ground_lateral_mm =  param_lateral_sur_k * ground_water_mm;
  ground_lateral_mm = ifelse(ground_lateral_mm > ground_diff_mm, ground_diff_mm, ground_lateral_mm) ;
  return ifelse(ground_lateral_mm > - ground_water_mm, ground_lateral_mm, - ground_water_mm) ;
}

//' @rdname lateral
//' @param param_lateral_sup_k,param_lateral_sup_gamma parameters for [lateral_SupplyPow()]
//' @export
// [[Rcpp::export]]
NumericVector lateral_SupplyPow(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector param_lateral_sup_k,
    NumericVector param_lateral_sup_gamma
)
{
  NumericVector ground_lateral_mm, k_;
  NumericVector ground_diff_mm = (ground_capacity_mm - ground_water_mm);
  
  k_ = param_lateral_sup_k * vecpow((ground_water_mm / ground_capacity_mm), param_lateral_sup_gamma);
  ground_lateral_mm = k_ * ground_water_mm;
  ground_lateral_mm = ifelse(ground_lateral_mm > ground_diff_mm, ground_diff_mm, ground_lateral_mm) ;
  return ifelse(ground_lateral_mm > - ground_water_mm, ground_lateral_mm, - ground_water_mm) ;
}
