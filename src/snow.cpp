#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]


//' **snow**
//' @name snow
//' @inheritParams all_vari
//' @param param_snow_kus_fE <0.0005, 0.003> (mm/m2/MJ) snow melt temperature parameter for [snowMelt_Factor()]
//' @param param_snow_kus_fT <0.05, 1> (mm/m2/h/Cel) potential melt volum per Cel per hour parameter for [snowMelt_Factor()]
//' @return snow_melt_mm (mm/m2) melted snow
//' @export
// [[Rcpp::export]]
NumericVector snowMelt_Kustas(
    NumericVector snow_ice_mm,
    NumericVector atmos_temperature_Cel,
    NumericVector atmos_netRadiat_MJ,
    NumericVector time_step_h,
    NumericVector param_snow_kus_fE,
    NumericVector param_snow_kus_fT
)
{
  NumericVector snow_melt_mm = ifelse(atmos_temperature_Cel < 0, 0, atmos_temperature_Cel) * param_snow_kus_fT * time_step_h + param_snow_kus_fE * atmos_netRadiat_MJ;
  return ifelse(snow_melt_mm > snow_ice_mm, snow_ice_mm, snow_melt_mm) ;
  
}

//' @rdname snow
//' @param param_snow_fac_Tmelt <0, 3> (Cel) snow melt temperature parameter for [snowMelt_Factor()]
//' @param param_snow_fac_f <0.05, 2> (mm/m2/h/Cel) potential melt volum per Cel per hour parameter for [snowMelt_Factor()]
//' @export
// [[Rcpp::export]]
NumericVector snowMelt_Factor(
    NumericVector snow_ice_mm,
    NumericVector atmos_temperature_Cel,
    NumericVector time_step_h,
    NumericVector param_snow_fac_f,
    NumericVector param_snow_fac_Tmelt
)
{
  NumericVector diff_T, snow_melt_mm;
  diff_T = atmos_temperature_Cel - param_snow_fac_Tmelt;
  diff_T = ifelse(diff_T > 0, diff_T, 0);
  
  snow_melt_mm = param_snow_fac_f * time_step_h * diff_T;
  return ifelse(snow_melt_mm > snow_ice_mm, snow_ice_mm, snow_melt_mm) ;
}

