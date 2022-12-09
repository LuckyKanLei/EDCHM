#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]


//' **infiltration**
//' @name infilt
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' Under the concept of the conceptional HM, the flux of infiltration always be calculated by the pounded water on the land \mjseqn{W_{land}}, 
//' which can be precipitation, precipitation after interception or precipitation with sonow melt and so on. 
//' The second point is the water acceptability of the soil layer (\mjseqn{C_{soil} - W_{soil}}).
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{F_{iflt} = f_{infilt}(D_{land}, D_{soil})}
//' 
//' 
//' to:
//' 
//' \mjsdeqn{F_{iflt} = f_{infilt}(W_{land}, W_{soil}, C_{soil}, ...)}
//' 
//' some methods will tread the infiltartion as the part of th pounded water so there is also:
//' 
//' \mjsdeqn{F_{iflt} = k^* W_{land}}
//' 
//' 
//' where
//' - \mjseqn{F_{iflt}} is `infilt_mm`
//' - \mjseqn{W_{land}} is `land_water_mm`
//' - \mjseqn{W_{soil}} is `soil_water_mm`
//' - \mjseqn{C_{soil}} is `soil_capacity_mm`
//' - \mjseqn{k^*} is estimated ratio.
//' 
//' The output density distribution from 9 methods:
//'
//' \if{html}{\figure{mdl_infilt.svg}}
//' \if{latex}{\figure{mdl_infilt.pdf}{options: width=140mm}}
//' @references
//' \insertAllCited{}
//' @return flux of infiltration from land surface to soil layer
//' 
//' @details
//' # **_GR4J** \insertCite{GR4J_Perrin_2003}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_infilt_gr4.svg}}
//' \if{latex}{\figure{mdl_infilt_gr4.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{iflt}=\frac{C_{soil}\left(1-\left(\frac{W_{soil}}{C_{soil}}\right)^{2}\right) \tanh \left(\frac{W_{land}}{C_{soil}}\right)}{1+\frac{W_{soil}}{C_{soil}} \tanh \left(\frac{W_{land}}{C_{soil}}\right)}}
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
//' @details
//' # **_UBC** \insertCite{UBC_Quick_1977}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_infilt_ubc.svg}}
//' \if{latex}{\figure{mdl_infilt_ubc.pdf}{options: width=140mm}}
//' 
//' estimate the ratio \mjseqn{k^*} as:
//' \mjsdeqn{k^* = p_{imper} 10^{\frac{W_{soil}-C_{soil}}{p_{AGEN}}}}
//' where
//'   - \mjseqn{p_{imper}} is `land_impermeableFrac_1`
//'   - \mjseqn{p_{AGEN}} is `param_infilt_ubc_P0AGEN`
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
//' @details
//' # **_SupplyRatio**: 
//'
//' \if{html}{\figure{mdl_infilt_sur.svg}}
//' \if{latex}{\figure{mdl_infilt_sur.pdf}{options: width=140mm}}
//' 
//' is a very simple method, which estimate only the pounded water:
//' \mjsdeqn{k^* = k}
//' where
//'   - \mjseqn{k} is `param_infilt_sur_k`
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
//' @details
//' # **_AcceptRatio**: 
//'
//' \if{html}{\figure{mdl_infilt_acr.svg}}
//' \if{latex}{\figure{mdl_infilt_acr.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{iflt} = k (C_{soil} - W_{soil})}
//' where
//'   - \mjseqn{k} is `param_infilt_acr_k`
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
//' @details
//' # **_SupplyPow**: 
//'
//' \if{html}{\figure{mdl_infilt_sup.svg}}
//' \if{latex}{\figure{mdl_infilt_sup.pdf}{options: width=140mm}}
//' 
//' is a very simple method, which estimate only the pounded water:
//' \mjsdeqn{F_{iflt} = kW_{land}^{\gamma}}
//' where
//'   - \mjseqn{k} is `param_infilt_sup_k`
//'   - \mjseqn{\gamma} is `param_infilt_sup_gamma`
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
//' @details
//' # **_AcceptPow**: 
//'
//' \if{html}{\figure{mdl_infilt_acp.svg}}
//' \if{latex}{\figure{mdl_infilt_acp.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{iflt} = k \left(\frac{C_{soil} - W_{soil}}{C_{soil}} \right)^{\gamma}}
//' where
//'   - \mjseqn{k} is `param_infilt_acp_k`
//'   - \mjseqn{\gamma} is `param_infilt_acp_gamma`
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
//' @details
//' # **_HBV** \insertCite{HBV_Lindstrom_1997}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_infilt_hbv.svg}}
//' \if{latex}{\figure{mdl_infilt_hbv.pdf}{options: width=140mm}}
//' 
//' estimate the ratio \mjseqn{k^*} as:
//' \mjsdeqn{k^* = 1-\left(\frac{W_{soil}}{C_{soil}}\right)^{\beta}}
//' where
//'   - \mjseqn{\beta} is `param_infilt_hbv_beta`
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
//' @details
//' # **_XAJ** \insertCite{XAJ_Zhao_1992}{EDCHM}:
//'
//' \if{html}{\figure{mdl_infilt_xaj.svg}}
//' \if{latex}{\figure{mdl_infilt_xaj.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{iflt} = MM  \frac{\left( \frac{MM - AU}{MM} \right)^{B+1} - \left( \frac{MM - AU - W_{land}}{MM} \right)^{B+1}}{B+1}}
//' \mjsdeqn{AU = MM - \left( \frac{(1 - W_{soil})(B+1)}{MM} \right)^{1 / B - 1}  }
//' \mjsdeqn{MM = C_{soil}(B+1)  }
//' where
//'   - \mjseqn{B} is `param_infilt_xaj_B`
//' 
//' ![](xaj_infilt.png)
//' 
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
//' @details
//' # **_VIC** \insertCite{VIC_Wood_1992}{EDCHM}:
//'
//' \if{html}{\figure{mdl_infilt_vic.svg}}
//' \if{latex}{\figure{mdl_infilt_vic.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{infilt} = \int_{i_{0}}^{i_{0}+P} A(i) {\rm d} i}
//' \mjsdeqn{i = C_{soil}(B+1) \left[ 1 - (1-A)^{1/B} \right]}
//' where
//'   - \mjseqn{B} is `param_infilt_vic_B`
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
