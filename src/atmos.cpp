#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]


//' caculate **snowfall**
//' @name atmosSnow
//' @inheritParams all_vari
//' @return atmos_snow_mm (mm/m2/TS) snowfall volum
//' @param param_atmos_thr_Ts <-1, 3> (Cel) threshold air temperature that snow, parameter for [atmosSnow_ThresholdT()]
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
//' @param param_atmos_ubc_A0FORM <0.01, 3> (Cel) threshold air temperature that snow, it can not equal or small than 0, parameter for [atmosSnow_UBC()]
//' @export
// [[Rcpp::export]]
NumericVector atmosSnow_UBC(
    NumericVector atmos_precipitation_mm, 
    NumericVector atmos_temperature_Cel, 
    NumericVector param_atmos_ubc_A0FORM
)
{
  NumericVector atmos_snow_mm;
  atmos_snow_mm = (1 - atmos_temperature_Cel / param_atmos_ubc_A0FORM) * atmos_precipitation_mm;
  atmos_snow_mm = ifelse(atmos_temperature_Cel <= 0, atmos_precipitation_mm, atmos_snow_mm);
  return ifelse(atmos_temperature_Cel > param_atmos_ubc_A0FORM, 0, atmos_snow_mm);
}

