#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]



//' **percolation**
//' @name percola
//' @inheritParams all_vari
//' @param param_percola_grf_k <0.01, 1> coefficient parameter for [percola_GR4Jfix()]
//' @return percola_mm (mm/m2)
//' @export
// [[Rcpp::export]]
NumericVector percola_GR4Jfix(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector param_percola_grf_k
) 
{
  return soil_water_mm * (1 - pow((1 + pow(param_percola_grf_k * soil_water_mm / soil_capacity_mm, 4)), -0.25));
}


//' @rdname percola
//' @export
// [[Rcpp::export]]
NumericVector percola_GR4J(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm
) 
{
  return soil_water_mm * (1 - pow((1 + pow(4.0/9.0 * soil_water_mm / soil_capacity_mm, 4)), -0.25));
}

//' @rdname percola
//' @param param_percola_sur_k <0.01, 1> coefficient parameter for [percola_SupplyRatio()]
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
//' @param param_percola_sup_k <0.01, 1> coefficient parameter for [percola_SupplyPow()]
//' @param param_percola_sup_gamma <0, 7> parameter for [percola_SupplyPow()]
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

//' @rdname percola
//' @param param_percola_map_gamma <0.1, 5> exponential parameter for [percola_MaxPow()]
//' @export
// [[Rcpp::export]]
NumericVector percola_MaxPow(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector soil_potentialPercola_mm,
    NumericVector param_percola_map_gamma
)
{
  NumericVector percola_;
  
  percola_ = soil_potentialPercola_mm * vecpow(soil_water_mm / soil_capacity_mm, param_percola_map_gamma);
  
  return ifelse(percola_ > soil_water_mm, soil_water_mm, percola_) ;
}

//' @rdname percola
//' @param param_percola_thp_thresh <0.1, 0.9> coefficient parameter for [percola_ThreshPow()]
//' @param param_percola_thp_gamma <0.1, 5> exponential parameter for [percola_ThreshPow()]
//' @export
// [[Rcpp::export]]
NumericVector percola_ThreshPow(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector soil_potentialPercola_mm,
    NumericVector param_percola_thp_thresh,
    NumericVector param_percola_thp_gamma
)
{
  NumericVector percola_, percola_temp;
  percola_temp = (soil_water_mm / soil_capacity_mm - param_percola_thp_thresh);
  percola_temp = ifelse(percola_temp < 0, 0, percola_temp);
  percola_ = soil_potentialPercola_mm * vecpow(percola_temp / (1 - param_percola_thp_thresh), param_percola_thp_gamma);
  percola_ = ifelse(percola_ > soil_potentialPercola_mm, soil_potentialPercola_mm, percola_);
  return ifelse(percola_ > soil_water_mm, soil_water_mm, percola_) ;
}


//' @rdname percola
//' @param param_percola_arn_thresh <0.1, 0.9> coefficient parameter for [percola_ThreshPow()]
//' @param param_percola_arn_k <0.1, 1> exponential parameter for [percola_ThreshPow()]
//' @export
// [[Rcpp::export]]
NumericVector percola_Arno(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector soil_potentialPercola_mm,
    NumericVector param_percola_arn_thresh,
    NumericVector param_percola_arn_k
)
{
  NumericVector percola_, percola_1, percola_2, Ws_Wc;
  Ws_Wc = soil_capacity_mm * param_percola_arn_thresh;
  percola_1 = param_percola_arn_k * soil_potentialPercola_mm / (soil_capacity_mm) * soil_water_mm;
  percola_2 = param_percola_arn_k * soil_potentialPercola_mm / (soil_capacity_mm) * soil_water_mm + soil_potentialPercola_mm * (1 - param_percola_arn_k) * pow((soil_water_mm - Ws_Wc) / (soil_capacity_mm - Ws_Wc),2);
  percola_ = ifelse(soil_water_mm < Ws_Wc, percola_1, percola_2);
  percola_ = ifelse(soil_potentialPercola_mm > Ws_Wc, soil_water_mm, percola_);
  percola_ = ifelse(percola_ > soil_potentialPercola_mm, soil_potentialPercola_mm, percola_);
  return ifelse(percola_ > soil_water_mm, soil_water_mm, percola_) ;
}


//' @rdname percola
//' @export
// [[Rcpp::export]]
NumericVector percola_BevenWood(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector soil_fieldCapacityPerc_1,
    NumericVector soil_potentialPercola_mm
)
{
  NumericVector soil_percola_mm, soil_percolaAvilibale_mm, soil_diff_mm, k_;
  soil_percolaAvilibale_mm = soil_water_mm - soil_capacity_mm * (1-soil_fieldCapacityPerc_1);
  soil_percolaAvilibale_mm = ifelse(soil_percolaAvilibale_mm < 0, 0, soil_percolaAvilibale_mm);
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  k_ = ifelse(soil_diff_mm < 1, soil_water_mm, soil_water_mm / soil_diff_mm);
  soil_percola_mm = k_ * soil_potentialPercola_mm;
  soil_percola_mm = ifelse(soil_water_mm > soil_percolaAvilibale_mm, soil_percola_mm, 0.0);
  return ifelse(soil_percola_mm > soil_percolaAvilibale_mm, soil_percolaAvilibale_mm, soil_percola_mm) ;
}

