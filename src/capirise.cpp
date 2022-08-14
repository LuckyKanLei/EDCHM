#include <Rcpp.h>
using namespace Rcpp;

NumericVector vecpow(NumericVector a, NumericVector b);

//' **capilarise**
//' @name capirise
//' @param ground_water_mm (mm/m2) water volum in `groundLy`
//' @param soil_water_mm (mm/m2) water volum in `soilLy`
//' @param soil_capacity_mm (mm/m2) average soil Capacity (maximal storage capacity)
//' @param param_capirise_sur_k parameters for[capirise_SupplyRatio()]
//' @return  capilarise (mm/m2)
// [[Rcpp::export]]
NumericVector capirise_SupplyRatio(
    NumericVector ground_water_mm,
    NumericVector soil_water_mm ,
    NumericVector soil_capacity_mm, 
    NumericVector param_capirise_sur_k
)
{
  NumericVector soil_diff_mm, capirise_mm, limit_mm;
  
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  
  capirise_mm = param_capirise_sur_k * ground_water_mm;
  
  limit_mm = ifelse(soil_diff_mm > ground_water_mm, ground_water_mm, soil_diff_mm) ;
  return ifelse(capirise_mm > limit_mm, limit_mm, capirise_mm) ;
}

//' @rdname capirise
//' @param param_capirise_sup_k,param_capirise_sup_gamma parameters for [capirise_SupplyPow()]
// [[Rcpp::export]]
NumericVector capirise_SupplyPow(
    NumericVector ground_water_mm,
    NumericVector soil_water_mm ,
    NumericVector soil_capacity_mm, 
    NumericVector param_capirise_sup_k,
    NumericVector param_capirise_sup_gamma
)
{
  NumericVector soil_diff_mm, capirise_mm, k_, limit_mm;
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  
  k_ = param_capirise_sup_k * vecpow((ground_water_mm / soil_capacity_mm), param_capirise_sup_gamma);
  capirise_mm = k_ * ground_water_mm;
  
  limit_mm = ifelse(soil_diff_mm > ground_water_mm, ground_water_mm, soil_diff_mm) ;
  return ifelse(capirise_mm > limit_mm, limit_mm, capirise_mm) ;
}

//' @rdname capirise
//' @param param_capirise_acr_k parameters [capirise_AcceptRatio()]
// [[Rcpp::export]]
NumericVector capirise_AcceptRatio(
    NumericVector ground_water_mm, 
    NumericVector soil_water_mm ,
    NumericVector soil_capacity_mm, 
    NumericVector param_capirise_acr_k
)
{
  NumericVector soil_diff_mm, capirise_mm, limit_mm;
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  
  capirise_mm = soil_diff_mm * param_capirise_acr_k;
  
  limit_mm = ifelse(soil_diff_mm > ground_water_mm, ground_water_mm, soil_diff_mm) ;
  return ifelse(capirise_mm > limit_mm, limit_mm, capirise_mm) ;
}



//' @rdname capirise
//' @param param_capirise_acp_k,param_capirise_acp_gamma parameters for [capirise_AcceptPow()]
// [[Rcpp::export]]
NumericVector capirise_AcceptPow(
    NumericVector ground_water_mm, 
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector param_capirise_acp_k,
    NumericVector param_capirise_acp_gamma
)
{
  NumericVector capirise_mm, k_, soil_diff_mm, limit_mm;
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  
  k_ = param_capirise_acp_k * vecpow((soil_water_mm / soil_capacity_mm), param_capirise_acp_gamma);
  capirise_mm = k_ * soil_diff_mm;
  
  limit_mm = ifelse(soil_diff_mm > ground_water_mm, ground_water_mm, soil_diff_mm) ;
  return ifelse(capirise_mm > limit_mm, limit_mm, capirise_mm) ;
}
