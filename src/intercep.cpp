#include <Rcpp.h>
#include "00utilis.h"
using namespace Rcpp;

//' **interception** water from land go into the soil.
//' @name intercep
//' @param atmos_rain_mm (mm/m2) preciptation in rain form
//' @param land_interceptWater_mm (mm/m2) initial water volum that can be intercepted
//' @param land_interceptCapacity_mm (mm/m2) average intercept Capacity (maximal storage capacity)
//' @return intercept_water_mm (mm/m2) intercepted water in this timestep
//' @export
// [[Rcpp::export]]
NumericVector intercep_Full(
    NumericVector atmos_rain_mm,
    NumericVector land_interceptWater_mm,
    NumericVector land_interceptCapacity_mm
)
{
  NumericVector water_diff_mm = land_interceptCapacity_mm - land_interceptWater_mm;
  return ifelse(water_diff_mm > atmos_rain_mm, atmos_rain_mm, water_diff_mm) ;
}
