#include <Rcpp.h>
#include "00utilis.h"
using namespace Rcpp;
// [[Rcpp::interfaces(r, cpp)]]


//' functions for **atmosphere**
//' @name atmos
//' @param time_dayOfYear_ (1, 366) the number of the day in the year between 1 (1 January) and 365 or 366 (31 December)
//' @param atmos_temperature_Cel (Cel) the average air temperature in the time phase
//' @param atmos_temperatureMax_Cel (Cel) the maximal air temperature in the time phase
//' @param atmos_temperatureMin_Cel (Cel) the minimal air temperature in the time phase
//' @param atmos_relativeHumidity_1 (0, 1) relative humidity
//' @param atmos_solarRadiat_MJ (MJ/m2/TS) the solar radiation that actually reaches the earths surface
//' @param atmos_windSpeed_m_s (m/s) measured wind speed at z m above ground surface
//' @param atmos_windMeasureHeight_m (m) height of measurement above ground surface
//' @param land_latitude_Degree (degree) average latitude
//' @param land_elevation_m (m) average elevation
//' @return atmos_netRadiat_MJ (MJ/m2/TS) the balance between the energy absorbed, reflected and emitted by the earths surface or the difference between the incoming net shortwave (Rns) and the net outgoing longwave (Rnl) radiation
//' @export
// [[Rcpp::export]]
NumericVector atmos_NettoRadiat(
    NumericVector time_dayOfYear_,
    NumericVector atmos_temperature_Cel, 
    NumericVector atmos_temperatureMax_Cel, 
    NumericVector atmos_temperatureMin_Cel, 
    NumericVector atmos_relativeHumidity_1,
    NumericVector atmos_solarRadiat_MJ,
    NumericVector land_latitude_Degree,
    NumericVector land_elevation_m
)
{
  NumericVector e_s, e_a, phi_, d_r, delta_, omega_s, R_a, R_so, R_ns, R_nl;
  // start const double
  const double alpha_ = 0.23; // Eq.38
  const double sigma_ = 4.903e-09; // Stefan-Boltzmann constant Eq.39
  // end const double
  
  //// e_s,
  e_s = 0.3054 * exp(17.27 * atmos_temperatureMax_Cel / (atmos_temperatureMax_Cel + 237.3)) + 0.3054 * exp(17.27 * atmos_temperatureMin_Cel / (atmos_temperatureMin_Cel + 237.3)); // Eq.11, 12
  
  //// e_a,
  e_a = e_s * atmos_relativeHumidity_1; // Eq.19
  
  //// R_n,
  phi_ = M_PI / 180 * land_latitude_Degree; // Eq.22
  d_r = 1 + 0.033 * cos(2 * M_PI / 365 * time_dayOfYear_); // Eq.23
  delta_ = 0.409 * sin(2 * M_PI / 365 * time_dayOfYear_ - 1.39); // Eq.24
  omega_s = acos(-tan(phi_) * tan(delta_)); // Eq.25
  //// (24 (60)/M_PI) G_sc = 37.58603, G_sc <- 0.0820 MJ / (m^2 d)
  R_a = 37.58603 * d_r * (omega_s * sin(phi_) * sin(delta_) + cos(phi_) * cos(delta_) * sin(omega_s)); // Eq.21
  
  R_so = (0.75 + 2e-5 * land_elevation_m) * R_a; // Eq.37
  R_ns = (1 - alpha_) * atmos_solarRadiat_MJ; // Eq.38
  R_nl = sigma_ * (pow((atmos_temperatureMax_Cel + 273.16), 4) + pow((atmos_temperatureMin_Cel + 273.16), 4)) / 2 * (0.34 - 0.14 * sqrt(e_a)) * (1.35 * atmos_solarRadiat_MJ / R_so - 0.35); // Eq.39
  return R_ns - R_nl; // Eq.40
}

//' @rdname  atmos
//' @return atmos_netRadiat_MJ	(MJ/m2/TS) the balance between the energy absorbed, reflected and emitted by the earth's surface or the difference between the incoming net shortwave (Rns) and the net outgoing longwave (Rnl) radiation
//' @export
// [[Rcpp::export]]
NumericVector atmos_SaturatVaporPress(
    NumericVector atmos_temperature_Cel 
)
{
  return 6.1078 * exp(17.27 * atmos_temperature_Cel / (atmos_temperature_Cel + 237.3)); // Eq.40
}

//' @rdname  atmos
//' @return atmos_vaporPress_hPa	(hPa) actual vapour pressure, can be calculated by [atmos_VaporPress()]
//' @export
// [[Rcpp::export]]
NumericVector atmos_VaporPress(
    NumericVector atmos_temperature_Cel,
    NumericVector atmos_relativeHumidity_1 
)
{
  return 6.1078 * exp(17.27 * atmos_temperature_Cel / (atmos_temperature_Cel + 237.3)) * atmos_relativeHumidity_1; // Eq.40
}

//' @rdname  atmos
//' @return atmos_WindSpeed2m_m_s (m/s) wind speed at 2 m above ground surface
//' @export
// [[Rcpp::export]]
NumericVector atmos_WindSpeed2m(
    NumericVector atmos_windSpeed_m_s,
    NumericVector atmos_windMeasureHeight_m
)
{
  return atmos_windMeasureHeight_m * 4.87 / (log(67.8 * atmos_windSpeed_m_s - 5.42)); // Eq.47
}


