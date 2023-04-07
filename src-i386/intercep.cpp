#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]


//' **interception** water from land go into the soil.
//' @name intercep
//' @inheritParams all_vari
//' @description 
//' \loadmathjax
//' 
//' In hydrological modeling, interception refers to the process by which water from precipitation is temporarily retained on the surfaces of vegetation, such as leaves and branches, before being returned to the atmosphere through evaporation or drip.
//' 
//' Under the concept of the conceptional HM, the interception will simply be calculated with the maximal interception of the land.
//' And the interception water will also not go to the land, but will be evaporated.
//' The maximal Interception of the canopy is maybe difficult to estimate 
//' but the process is really simple and there is also not so many method to describe it. 
//' 
//' @details
//' # **_Full** : 
//' 
//'
//' \if{html}{\figure{mdl_intercep_ful.svg}}
//' \if{latex}{\figure{mdl_intercep_ful.pdf}{options: width=140mm}}
//' 
//' consider only the radiation and temperature as the main factors. 
//' \mjsdeqn{F_{itcp} = C_{icpt} - W_{icpt}}
//' where
//'   - \mjseqn{F_{icp}} is `intercept_water_mm`
//'   - \mjseqn{C_{icpt}} is `land_intercepCapaciy_mm`
//'   - \mjseqn{W_{icpt}} is `land_intercepWater_mm`
//' @return intercept_water_mm (mm/m2) intercepted water in this timestep
//' @export
// [[Rcpp::export]]
NumericVector intercep_Full(
    NumericVector atmos_precipitation_mm,
    NumericVector land_interceptWater_mm,
    NumericVector land_interceptCapacity_mm
)
{
  NumericVector water_diff_mm = land_interceptCapacity_mm - land_interceptWater_mm;
  return ifelse(water_diff_mm > atmos_precipitation_mm, atmos_precipitation_mm, water_diff_mm) ;
}
