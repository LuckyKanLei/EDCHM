#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]


//' caculate **snowfall**
//' @name atmosSnow
//' @inheritParams all_vari
//' @param param_atmos_thr_Ts (Cel) parameters for [atmosSnow_ThresholdT()]
//' @return atmos_snow_mm (mm/m2/TS) snowfall volum
//' @export
// [[Rcpp::export]]
NumericVector atmosSnow_ThresholdT(
    NumericVector atmos_precipitation_mm, 
    NumericVector atmos_temperature_Cel, 
    NumericVector param_atmos_thr_Ts
)
{
  return ifelse(atmos_temperature_Cel > param_atmos_thr_Ts, 0, atmos_precipitation_mm);
}

//' @rdname atmosSnow
//' @param param_atmos_ubc_A0FORM (0, 3 Cel) parameters for [atmosSnow_UBC()]
//' @export
// [[Rcpp::export]]
NumericVector atmosSnow_UBC(
    NumericVector atmos_precipitation_mm, 
    NumericVector atmos_temperature_Cel, 
    NumericVector param_atmos_ubc_A0FORM
)
{
  NumericVector atmos_snow_mm;
  atmos_snow_mm = atmos_temperature_Cel / param_atmos_ubc_A0FORM * atmos_precipitation_mm;
  atmos_snow_mm = ifelse(atmos_temperature_Cel < 0, atmos_precipitation_mm, atmos_snow_mm);
  return ifelse(atmos_temperature_Cel > param_atmos_ubc_A0FORM, 0, atmos_snow_mm);
}

