#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]



//' **lateral flux**
//' @name lateral
//' @inheritParams all_vari
//' @return lateral_mm (mm/m2)
//' @export
// [[Rcpp::export]]
NumericVector lateral_GR4J(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector ground_potentialLateral_mm
) 
{
  NumericVector ground_lateral_mm;
  NumericVector ground_diff_mm = (ground_capacity_mm - ground_water_mm);
  ground_lateral_mm = ground_potentialLateral_mm * pow((ground_water_mm / ground_capacity_mm), 3.5);
  ground_lateral_mm = ifelse(ground_lateral_mm > ground_diff_mm, ground_diff_mm, ground_lateral_mm) ;
  return ifelse(ground_lateral_mm > - ground_water_mm, ground_lateral_mm, - ground_water_mm) ;
}


//' @rdname lateral
//' @param param_lateral_grf_gamma <0.01, 5> parameter for [lateral_GR4Jfix()]
//' @export
// [[Rcpp::export]]
NumericVector lateral_GR4Jfix(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector ground_potentialLateral_mm,
    NumericVector param_lateral_grf_gamma
) 
{
  NumericVector ground_lateral_mm;
  NumericVector ground_diff_mm = (ground_capacity_mm - ground_water_mm);
  ground_lateral_mm = ground_potentialLateral_mm * vecpow((ground_water_mm / ground_capacity_mm), param_lateral_grf_gamma);
  ground_lateral_mm = ifelse(ground_lateral_mm > ground_diff_mm, ground_diff_mm, ground_lateral_mm) ;
  return ifelse(ground_lateral_mm > - ground_water_mm, ground_lateral_mm, - ground_water_mm) ;
}

//' @rdname lateral
//' @param param_lateral_sur_k <-2, 1> coefficient parameter for [lateral_SupplyRatio()]
//' @export
// [[Rcpp::export]]
NumericVector lateral_SupplyRatio(
    NumericVector ground_water_mm,
    NumericVector param_lateral_sur_k
)
{
  
  NumericVector ground_lateral_mm =  param_lateral_sur_k * ground_water_mm;
  return ifelse(ground_lateral_mm > - ground_water_mm, ground_lateral_mm, - ground_water_mm) ;
}

//' @rdname lateral
//' @param param_lateral_sup_k <-1, 1> coefficient parameter for [lateral_SupplyPow()]
//' @param param_lateral_sup_gamma <0.01, 5> parameters for [lateral_SupplyPow()]
//' @export
// [[Rcpp::export]]
NumericVector lateral_SupplyPow(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector param_lateral_sup_k,
    NumericVector param_lateral_sup_gamma
)
{
  NumericVector ground_lateral_mm, k_;
  NumericVector ground_diff_mm = (ground_capacity_mm - ground_water_mm);
  
  k_ = param_lateral_sup_k * vecpow((ground_water_mm / ground_capacity_mm), param_lateral_sup_gamma);
  ground_lateral_mm = k_ * ground_water_mm;
  ground_lateral_mm = ifelse(ground_lateral_mm > ground_diff_mm, ground_diff_mm, ground_lateral_mm) ;
  return ifelse(ground_lateral_mm > - ground_water_mm, ground_lateral_mm, - ground_water_mm) ;
}

//' @rdname lateral
//' @param param_lateral_map_gamma <0.1, 5> exponential parameter for [lateral_MaxPow()]
//' @export
// [[Rcpp::export]]
NumericVector lateral_MaxPow(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector ground_potentialLateral_mm,
    NumericVector param_lateral_map_gamma
)
{
  NumericVector ground_lateral_mm;
  NumericVector ground_diff_mm = (ground_capacity_mm - ground_water_mm);
  
  ground_lateral_mm = ground_potentialLateral_mm * vecpow(ground_water_mm / ground_capacity_mm, param_lateral_map_gamma);
  
  ground_lateral_mm = ifelse(ground_lateral_mm > ground_diff_mm, ground_diff_mm, ground_lateral_mm) ;
  return ifelse(ground_lateral_mm > - ground_water_mm, ground_lateral_mm, - ground_water_mm) ;
}

//' @rdname lateral
//' @param param_lateral_thp_thresh <0.1, 0.9> coefficient parameter for [lateral_ThreshPow()]
//' @param param_lateral_thp_gamma <0.1, 5> exponential parameter for [lateral_ThreshPow()]
//' @export
// [[Rcpp::export]]
NumericVector lateral_ThreshPow(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector ground_potentialLateral_mm,
    NumericVector param_lateral_thp_thresh,
    NumericVector param_lateral_thp_gamma
)
{
  NumericVector ground_lateral_mm, lateral_temp;
  NumericVector ground_diff_mm = (ground_capacity_mm - ground_water_mm);
  lateral_temp = (ground_water_mm / ground_capacity_mm - param_lateral_thp_thresh);
  lateral_temp = ifelse(lateral_temp < 0, 0, lateral_temp);
  
  ground_lateral_mm = ground_potentialLateral_mm * vecpow(lateral_temp / (1 - param_lateral_thp_thresh), param_lateral_thp_gamma);

  ground_lateral_mm = ifelse(ground_lateral_mm > ground_diff_mm, ground_diff_mm, ground_lateral_mm) ;
  return ifelse(ground_lateral_mm > - ground_water_mm, ground_lateral_mm, - ground_water_mm) ;
}


//' @rdname lateral
//' @param param_lateral_arn_thresh <0.1, 0.9> coefficient parameter for [lateral_ThreshPow()]
//' @param param_lateral_arn_k <0.1, 1> exponential parameter for [lateral_ThreshPow()]
//' @export
// [[Rcpp::export]]
NumericVector lateral_Arno(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector ground_potentialLateral_mm,
    NumericVector param_lateral_arn_thresh,
    NumericVector param_lateral_arn_k
)
{
  NumericVector ground_lateral_mm, lateral_1, lateral_2, Ws_Wc;
  NumericVector ground_diff_mm = (ground_capacity_mm - ground_water_mm);
  Ws_Wc = ground_capacity_mm * param_lateral_arn_thresh;
  lateral_1 = param_lateral_arn_k * ground_potentialLateral_mm / (ground_capacity_mm) * ground_water_mm;
  lateral_2 = param_lateral_arn_k * ground_potentialLateral_mm / (ground_capacity_mm) * ground_water_mm + ground_potentialLateral_mm * (1 - param_lateral_arn_k) * pow((ground_water_mm - Ws_Wc) / (ground_capacity_mm - Ws_Wc),2);
  ground_lateral_mm = ifelse(ground_water_mm < Ws_Wc, lateral_1, lateral_2);
  ground_lateral_mm = ifelse(ground_potentialLateral_mm > Ws_Wc, ground_water_mm, ground_lateral_mm);

  ground_lateral_mm = ifelse(ground_lateral_mm > ground_diff_mm, ground_diff_mm, ground_lateral_mm) ;
  return ifelse(ground_lateral_mm > - ground_water_mm, ground_lateral_mm, - ground_water_mm) ;
}
