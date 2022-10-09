#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]



//' **baseflow**
//' @name baseflow
//' @inheritParams all_vari
//' @return ground_baseflow_mm (mm/m2/TS) 
//' @export
// [[Rcpp::export]]
NumericVector baseflow_GR4J(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm
)
{
  NumericVector baseflow_, k_;
  
  k_ = 1 - pow((1 + pow(ground_water_mm / ground_capacity_mm, 4)), -0.25);
  baseflow_ = k_ * ground_water_mm;
  
  return ifelse(baseflow_ > ground_water_mm, ground_water_mm, baseflow_) ;
}

//' @rdname baseflow
//' @param param_baseflow_grf_gamma <2, 7> exponential parameter for [baseflow_GR4Jfix()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_GR4Jfix(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector param_baseflow_grf_gamma
)
{
  NumericVector baseflow_, k_;
  
  k_ = 1 - vecpow((1 + vecpow(ground_water_mm / ground_capacity_mm, param_baseflow_grf_gamma)), -1.0 / param_baseflow_grf_gamma);
  baseflow_ = k_ * ground_water_mm;
  
  return ifelse(baseflow_ > ground_water_mm, ground_water_mm, baseflow_) ;
}



//' @rdname baseflow
//' @param param_baseflow_sup_k <0.01, 1> coefficient parameter for [baseflow_SupplyPow()]
//' @param param_baseflow_sup_gamma <0, 1> exponential parameter for [baseflow_SupplyPow()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_SupplyPow(
    NumericVector ground_water_mm,
    NumericVector param_baseflow_sup_k,
    NumericVector param_baseflow_sup_gamma
)
{
  NumericVector baseflow_;
  
  baseflow_ = param_baseflow_sup_k * vecpow(ceil(ground_water_mm), param_baseflow_sup_gamma);

  return ifelse(baseflow_ > ground_water_mm, ground_water_mm, baseflow_) ;
}

//' @rdname baseflow
//' @param param_baseflow_sur_k <0.01, 1> coefficient parameter for [baseflow_SupplyRatio()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_SupplyRatio(
    NumericVector ground_water_mm,
    NumericVector param_baseflow_sur_k
)
{
  
  return param_baseflow_sur_k * ground_water_mm;
  
}

//' @rdname baseflow
//' @param param_baseflow_map_gamma <0.1, 5> exponential parameter for [baseflow_MaxPow()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_MaxPow(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector ground_potentialBaseflow_mm,
    NumericVector param_baseflow_map_gamma
)
{
  NumericVector baseflow_;
  
  baseflow_ = ground_potentialBaseflow_mm * vecpow(ground_water_mm / ground_capacity_mm, param_baseflow_map_gamma);
  
  return ifelse(baseflow_ > ground_water_mm, ground_water_mm, baseflow_) ;
}

//' @rdname baseflow
//' @param param_baseflow_thp_thresh <0.1, 0.9> coefficient parameter for [baseflow_ThreshPow()]
//' @param param_baseflow_thp_gamma <0.1, 5> exponential parameter for [baseflow_ThreshPow()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_ThreshPow(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector ground_potentialBaseflow_mm,
    NumericVector param_baseflow_thp_thresh,
    NumericVector param_baseflow_thp_gamma
)
{
  NumericVector baseflow_, baseflow_temp;
  baseflow_temp = (ground_water_mm / ground_capacity_mm - param_baseflow_thp_thresh);
  baseflow_temp = ifelse(baseflow_temp < 0, 0, baseflow_temp);
  
  baseflow_ = ground_potentialBaseflow_mm * vecpow(baseflow_temp / (1 - param_baseflow_thp_thresh), param_baseflow_thp_gamma);
  baseflow_ = ifelse(baseflow_ > ground_potentialBaseflow_mm, ground_potentialBaseflow_mm, baseflow_);
  
  return ifelse(baseflow_ > ground_water_mm, ground_water_mm, baseflow_) ;
}


//' @rdname baseflow
//' @param param_baseflow_arn_thresh <0.1, 0.9> coefficient parameter for [baseflow_ThreshPow()]
//' @param param_baseflow_arn_k <0.1, 1> exponential parameter for [baseflow_ThreshPow()]
//' @export
// [[Rcpp::export]]
NumericVector baseflow_Arno(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector ground_potentialBaseflow_mm,
    NumericVector param_baseflow_arn_thresh,
    NumericVector param_baseflow_arn_k
)
{
  NumericVector baseflow_, baseflow_1, baseflow_2, Ws_Wc;
  Ws_Wc = ground_capacity_mm * param_baseflow_arn_thresh;
  baseflow_1 = param_baseflow_arn_k * ground_potentialBaseflow_mm / (ground_capacity_mm) * ground_water_mm;
  baseflow_2 = param_baseflow_arn_k * ground_potentialBaseflow_mm / (ground_capacity_mm) * ground_water_mm + ground_potentialBaseflow_mm * (1 - param_baseflow_arn_k) * pow((ground_water_mm - Ws_Wc) / (ground_capacity_mm - Ws_Wc),2);
  baseflow_ = ifelse(ground_water_mm < Ws_Wc, baseflow_1, baseflow_2);
  baseflow_ = ifelse(ground_potentialBaseflow_mm > Ws_Wc, ground_water_mm, baseflow_);
  baseflow_ = ifelse(baseflow_ > ground_potentialBaseflow_mm, ground_potentialBaseflow_mm, baseflow_);
  return ifelse(baseflow_ > ground_water_mm, ground_water_mm, baseflow_) ;
}
