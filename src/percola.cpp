#include <Rcpp.h>
#include "00utilis.h"
using namespace Rcpp;
// [[Rcpp::interfaces(r, cpp)]]


NumericVector vecpow(NumericVector base, NumericVector exp);

//' **percolation**
//' @name percola
//' @param soil_water_mm (mm/m2) water volum in `soilLy`
//' @param soil_capacity_mm (mm/m2) average soil Capacity (maximal storage capacity)
//' @param param_percola_gr4_k parameters
//' @return percola_mm (mm/m2)
//' @export
// [[Rcpp::export]]
NumericVector percola_GR4J(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector param_percola_gr4_k
) 
{
  return soil_water_mm * (1 - pow((1 + pow(param_percola_gr4_k * soil_water_mm / soil_capacity_mm, 4)), -0.25));
}

//' @rdname percola
//' @param param_percola_sur_k parameters for [percola_SupplyRatio()]
//' @export
// [[Rcpp::export]]
NumericVector percola_SupplyRatio(
    NumericVector soil_water_mm,
    NumericVector param_percola_sur_k
)
{
  
  return param_percola_sur_k * soil_water_mm;
}

//' @rdname percola
//' @param param_percola_sup_k,param_percola_sup_gamma parameters for [percola_SupplyPow()]
//' @export
// [[Rcpp::export]]
NumericVector percola_SupplyPow(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector param_percola_sup_k,
    NumericVector param_percola_sup_gamma
)
{
  NumericVector soil_percola_mm, k_;
  
  k_ = param_percola_sup_k * vecpow((soil_water_mm / soil_capacity_mm), param_percola_sup_gamma);
  soil_percola_mm = k_ * soil_water_mm;
  return ifelse(soil_percola_mm > soil_water_mm, soil_water_mm, soil_percola_mm) ;
}
