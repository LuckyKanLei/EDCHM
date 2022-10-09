#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]


//' **infiltration**
//' @name infilt
//' @inheritParams all_vari
//' @export
// [[Rcpp::export]]
NumericVector infilt_GR4J(
    NumericVector land_water_mm,
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm
) 
{
  NumericVector soil_diff_mm, tanh_pn_x1, s_x1, infilt_water_mm, limit_mm;
  
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  limit_mm = ifelse(soil_diff_mm > land_water_mm, land_water_mm, soil_diff_mm);
  
  tanh_pn_x1 = tanh(land_water_mm / soil_capacity_mm);
  s_x1 = soil_water_mm / soil_capacity_mm;
  infilt_water_mm = soil_capacity_mm * (1 - (s_x1) * (s_x1)) * tanh_pn_x1 / (1 + s_x1 * tanh_pn_x1); //// Eq.3
  
  return ifelse(infilt_water_mm > limit_mm, limit_mm, infilt_water_mm);
}

//' @rdname infilt
//' @param param_infilt_ubc_P0AGEN <0.1, 4> coefficient parameter for [infilt_UBC()]
//' @export
// [[Rcpp::export]]
NumericVector infilt_UBC(
    NumericVector land_water_mm, 
    NumericVector land_impermeableFrac_1, 
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm, 
    NumericVector param_infilt_ubc_P0AGEN
)
{
  NumericVector soil_diff_mm, k_, infilt_water_mm, limit_mm;
  
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  limit_mm = ifelse(soil_diff_mm > land_water_mm, land_water_mm, soil_diff_mm);
  
  k_ = (1 - land_impermeableFrac_1 * vecpow10(- soil_diff_mm / (soil_capacity_mm * param_infilt_ubc_P0AGEN)));
  infilt_water_mm = land_water_mm * k_;
  
  return ifelse(infilt_water_mm > limit_mm, limit_mm, infilt_water_mm);
}

//' @rdname infilt
//' @param param_infilt_sur_k <0.01, 1> coefficient parameter for [infilt_SupplyRatio()]
//' @return infilt_mm (mm/m2) 
//' @export
// [[Rcpp::export]]
NumericVector infilt_SupplyRatio(
    NumericVector land_water_mm,
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm, 
    NumericVector param_infilt_sur_k
)
{
  NumericVector soil_diff_mm, infilt_water_mm, limit_mm;
  
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  limit_mm = ifelse(soil_diff_mm > land_water_mm, land_water_mm, soil_diff_mm);
  
  infilt_water_mm = param_infilt_sur_k * land_water_mm;
  
  return ifelse(infilt_water_mm > limit_mm, limit_mm, infilt_water_mm);
}

//' @rdname infilt
//' @param param_infilt_sup_k <0.01, 1> coefficient parameter for [infilt_SupplyPow()]
//' @param param_infilt_sup_gamma <0, 1> parameters for [infilt_SupplyPow()]
//' @export
// [[Rcpp::export]]
NumericVector infilt_SupplyPow(
    NumericVector land_water_mm,
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm, 
    NumericVector param_infilt_sup_k,
    NumericVector param_infilt_sup_gamma
)
{
  NumericVector soil_diff_mm, infilt_water_mm, k_, limit_mm;
  
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  limit_mm = ifelse(soil_diff_mm > land_water_mm, land_water_mm, soil_diff_mm);
  
  infilt_water_mm = param_infilt_sup_k * vecpow(ceil(land_water_mm), param_infilt_sup_gamma);
  infilt_water_mm = ifelse(infilt_water_mm > land_water_mm, land_water_mm, infilt_water_mm);
  return ifelse(infilt_water_mm > limit_mm, limit_mm, infilt_water_mm);
}

//' @rdname infilt
//' @param param_infilt_acr_k <0.01, 1> coefficient parameter for [infilt_AcceptRatio()]
//' @export
// [[Rcpp::export]]
NumericVector infilt_AcceptRatio(
    NumericVector land_water_mm, 
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm, 
    NumericVector param_infilt_acr_k
)
{
  NumericVector soil_diff_mm, infilt_water_mm, limit_mm;
  
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  limit_mm = ifelse(soil_diff_mm > land_water_mm, land_water_mm, soil_diff_mm);
  
  infilt_water_mm = soil_diff_mm * param_infilt_acr_k;
  
  return ifelse(infilt_water_mm > limit_mm, limit_mm, infilt_water_mm);
}



//' @rdname infilt
//' @param param_infilt_acp_k <0.01, 1> coefficient parameter for [infilt_AcceptPow()]
//' @param param_infilt_acp_gamma <0.001, 5> parameters for [infilt_AcceptPow()]
//' @export
// [[Rcpp::export]]
NumericVector infilt_AcceptPow(
    NumericVector land_water_mm, 
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector param_infilt_acp_k,
    NumericVector param_infilt_acp_gamma
)
{
  NumericVector infilt_water_mm, k_, soil_diff_mm, limit_mm;
  
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  limit_mm = ifelse(soil_diff_mm > land_water_mm, land_water_mm, soil_diff_mm);
  
  k_ = param_infilt_acp_k * vecpow((soil_diff_mm / soil_capacity_mm), param_infilt_acp_gamma);
  infilt_water_mm = k_ * soil_diff_mm;
  return ifelse(infilt_water_mm > land_water_mm, land_water_mm, infilt_water_mm);
}

//' @rdname infilt
//' @param param_infilt_hbv_beta <0.001, 5> parameters for [infilt_HBV()]
//' @export
// [[Rcpp::export]]
NumericVector infilt_HBV(
    NumericVector land_water_mm, 
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm, 
    NumericVector param_infilt_hbv_beta 
)
{
  NumericVector soil_diff_mm, infilt_water_mm, k_, limit_mm;
  
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  limit_mm = ifelse(soil_diff_mm > land_water_mm, land_water_mm, soil_diff_mm);
  
  k_ = (1 - vecpow(soil_water_mm / soil_capacity_mm, param_infilt_hbv_beta));
  
  infilt_water_mm = land_water_mm * k_;
  
  return ifelse(infilt_water_mm > limit_mm, limit_mm, infilt_water_mm);
}



//' @rdname infilt
//' @param param_infilt_xaj_B <0.01, 3> parameters for [infilt_XAJ()]
//' @export
// [[Rcpp::export]]
NumericVector infilt_XAJ(
    NumericVector land_water_mm, 
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm, 
    NumericVector param_infilt_xaj_B
)
{
  NumericVector soil_diff_mm, MM_, k_, infilt_water_mm, limit_mm, AU_, AU_L_MM, MM_AU, B_1, B_B_1, B_p_1;
  
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  limit_mm = ifelse(soil_diff_mm > land_water_mm, land_water_mm, soil_diff_mm) ;
  
  MM_ = soil_capacity_mm * (param_infilt_xaj_B + 1);
  
  B_p_1 = (param_infilt_xaj_B + 1);
  B_B_1 = param_infilt_xaj_B / B_p_1;
  B_1 = 1 / param_infilt_xaj_B;
  
  
  AU_ = MM_ * (1 - vecpow(1 - soil_water_mm * B_p_1 / MM_, B_1));
  
  AU_L_MM = (MM_ - AU_ - land_water_mm) / MM_;
  AU_L_MM = ifelse(AU_L_MM < 0, 0, AU_L_MM) ;
  MM_AU = (MM_ - AU_) / MM_;
  
  infilt_water_mm = - MM_ * (vecpow(AU_L_MM, B_p_1) - vecpow(MM_AU, B_p_1)) / B_p_1;
  
  return ifelse(infilt_water_mm > limit_mm, limit_mm, infilt_water_mm) ;
}

//' @rdname infilt
//' @param param_infilt_vic_B <0.01, 3> parameters for [infilt_VIC()]
//' @export
// [[Rcpp::export]]
NumericVector infilt_VIC(
    NumericVector land_water_mm, 
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm, 
    NumericVector param_infilt_vic_B
)
{
  NumericVector infilt_water_mm, limit_mm, soil_diff_mm, i_0, i_m, B_1, B_p_1, A_s;
  
  i_m = soil_capacity_mm * (param_infilt_vic_B + 1);
  
  B_p_1 = (param_infilt_vic_B + 1);
  B_1 = 1 / B_p_1;
  
  i_0 = i_m * (1 - vecpow(1 - soil_water_mm / soil_capacity_mm, B_1));
  
  
  soil_diff_mm = soil_capacity_mm - soil_water_mm;
  infilt_water_mm = ifelse((i_0 + land_water_mm) > i_m, soil_diff_mm, 
                           soil_diff_mm - soil_capacity_mm * vecpow((1 - (i_0 + land_water_mm) / i_m), B_p_1)) ;
  // A_s = 1 - vecpow((1 - i_0 / i_m), param_infilt_vic_B);
  // infilt_water_mm = i_0 * (1 - A_s);
  
  limit_mm = ifelse(soil_diff_mm > land_water_mm, land_water_mm, soil_diff_mm) ;
  return ifelse(infilt_water_mm > limit_mm, limit_mm, infilt_water_mm);
}
