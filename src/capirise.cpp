#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]



//' **capilarise**
//' @name capirise
//' @inheritParams all_vari
//' @param param_capirise_acr_k <0.01, 1> coefficient parameter [capirise_AcceptRatio()]
//' @export
// [[Rcpp::export]]
NumericVector capirise_AcceptRatio(
    NumericVector ground_water_mm, 
    NumericVector soil_water_mm ,
    NumericVector soil_capacity_mm, 
    NumericVector soil_fieldCapacityPerc_1,
    NumericVector param_capirise_acr_k
)
{
  NumericVector soil_diff_mm, capirise_mm, limit_mm;
  soil_diff_mm = soil_capacity_mm * (1 - soil_fieldCapacityPerc_1) - soil_water_mm;
  soil_diff_mm = ifelse(soil_diff_mm < 0, 0, soil_diff_mm);
  capirise_mm = soil_diff_mm * param_capirise_acr_k;
  
  limit_mm = ifelse(soil_diff_mm > ground_water_mm, ground_water_mm, soil_diff_mm) ;
  return ifelse(capirise_mm > limit_mm, limit_mm, capirise_mm) ;
}



//' @rdname capirise
//' @param param_capirise_acp_k <0.01, 1> coefficient parameter for [capirise_AcceptPow()]
//' @param param_capirise_acp_gamma <0.01, 1> exponential parameter for [capirise_AcceptPow()]
//' @export
// [[Rcpp::export]]
NumericVector capirise_AcceptPow(
    NumericVector ground_water_mm, 
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector soil_fieldCapacityPerc_1,
    NumericVector param_capirise_acp_k,
    NumericVector param_capirise_acp_gamma
)
{
  NumericVector capirise_mm, k_, soil_diff_mm, limit_mm;
  soil_diff_mm = soil_capacity_mm * (1 - soil_fieldCapacityPerc_1) - soil_water_mm;
  soil_diff_mm = ifelse(soil_diff_mm < 0, 0, soil_diff_mm);
  
  k_ = param_capirise_acp_k * vecpow((soil_diff_mm / (soil_capacity_mm * (1 - soil_fieldCapacityPerc_1))), param_capirise_acp_gamma);
  capirise_mm = k_ * soil_diff_mm;
  capirise_mm = ifelse(capirise_mm < 0, 0, capirise_mm);
  
  limit_mm = ifelse(soil_diff_mm > ground_water_mm, ground_water_mm, soil_diff_mm) ;
  return ifelse(capirise_mm > limit_mm, limit_mm, capirise_mm) ;
}
