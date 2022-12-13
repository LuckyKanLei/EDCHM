#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]


//' caculate **snowfall**
//' @name atmosSnow
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' Under the concept of the conceptional HM, the snowfall is always calculated by 
//' the temperature \mjseqn{T} and 
//' the precipitation availability, the portion of snowfall is always decided by the air tempature.
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{P_s = f_{atmosSnow}(D_{atms})}
//' 
//' 
//' to:
//' 
//' \mjsdeqn{P_s = f_{atmosSnow}(P, T) = k^*P}
//' \mjsdeqn{0 \leq k^* \leq 1}
//' where
//'   - \mjseqn{P} is `atmos_precpitation_mm`
//'   - \mjseqn{T} is `atmos_teperature_Cel`
//' - \mjseqn{k^*} is estimated portion
//' 
//' Then the different `atmosSnow` methods will estimate the portion \mjseqn{k^*}.
//' 
//' 
//' The output density distribution from 2 methods:
//' 
//' \if{html}{\figure{mdl_atmosSnow.svg}}
//' \if{latex}{\figure{mdl_atmosSnow.pdf}{options: width=140mm}}
//' @references
//' \insertAllCited{}
//' @return atmos_snow_mm (mm/m2/TS) snowfall volume
//' @details
//' # **_ThresholdT**: 
//' 
//' \if{html}{\figure{mdl_atmosSnow_thr.svg}}
//' \if{latex}{\figure{mdl_atmosSnow_thr.pdf}{options: width=140mm}}
//' 
//' Only a temperature is as the threshold defined, so estimate the portion \mjseqn{k^*} as: 
//' \mjsdeqn{k^{*}=1, \quad T \leq T_s}
//' where
//'   - \mjseqn{T_s} is `param_atmos_thr_Ts`
//' 
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
//' @details
//' # **_UBC** \insertCite{UBC_Quick_1977}{EDCHM}: 
//' 
//' \if{html}{\figure{mdl_atmosSnow_ubc.svg}}
//' \if{latex}{\figure{mdl_atmosSnow_ubc.pdf}{options: width=140mm}}
//' 
//' estimate the portion \mjseqn{k^*}{} as:
//' \mjsdeqn{k^* = 1- \frac{T}{T_0}}
//' \mjsdeqn{k^* \geq 0}
//' where
//'   - \mjseqn{T_0} is `param_atmos_ubc_A0FORM`
//' 
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
