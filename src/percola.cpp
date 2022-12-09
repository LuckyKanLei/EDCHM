#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]



//' **percolation**
//' @name percola
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' Under the concept of the conceptional HM, the flux of capillary rise always 
//' be calculated by the water in the soil layer \mjseqn{W_{soil}},
//' it can also be tread as the part of the \mjseqn{W_{soil}}.
//' So we can give the function from:
//' 
//' \mjsdeqn{F_{prcl} = f_{percola}(D_{grnd}, D_{soil})}
//' 
//' 
//' to:
//' 
//' \mjsdeqn{F_{prcl} = f_{percola}(W_{soil}, C_{soil}, W_{grnd}, C_{grnd}, ...) = k^* W_{soil}}
//' \mjsdeqn{F_{prcl} \leq W_{soil}}
//' \mjsdeqn{F_{prcl} \leq C_{grnd} - W_{grnd}}
//' 
//' 
//' where
//' - \mjseqn{F_{prcl}} is `soil_percola_mm`
//' - \mjseqn{W_{soil}} is `water_soil_mm`
//' - \mjseqn{C_{soil}} is `capacity_soil_mm`
//' - \mjseqn{W_{grnd}} is `ground_water_mm`
//' - \mjseqn{C_{grnd}} is `capacity_water_mm`
//' - \mjseqn{k^*} is estimated ratio
//' 
//' The output density distribution from 8 methods:
//'
//' \if{html}{\figure{mdl_percola.svg}}
//' \if{latex}{\figure{mdl_percola.pdf}{options: width=140mm}}
//' @references
//' \insertAllCited{}
//' @return percola_mm (mm/m2)
//' @details
//' # **_GR4J** \insertCite{GR4J_Perrin_2003}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_percola_gr4.svg}}
//' \if{latex}{\figure{mdl_percola_gr4.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{k^* = 1 - \left[ 1 + \left(\frac{4}{9} \frac{W_{soil}}{C_{soil}} \right)^4 \right]^{-1/4}}
//' where
//'   - \mjseqn{k^*} is estimated ratio
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
//' @details
//' # **_GR4Jfix** \insertCite{GR4J_Perrin_2003}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_percola_grf.svg}}
//' \if{latex}{\figure{mdl_percola_grf.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{k^* = 1 - \left[ 1 + \left(k \frac{W_{soil}}{C_{soil}} \right)^4 \right]^{-1/4}}
//' where
//'   - \mjseqn{k} is `param_percola_grf_k`
//' @param param_percola_grf_k <0.01, 1> coefficient parameter for [percola_GR4Jfix()]
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
//' @details
//' # **_MaxPow**: 
//'
//' \if{html}{\figure{mdl_percola_map.svg}}
//' \if{latex}{\figure{mdl_percola_map.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{prcl} = M_{prcl} \left(\frac{W_{soil}}{C_{soil}} \right)^\gamma}
//' where
//'   - \mjseqn{M_{prcl}} is `soil_potentialPercola_mm`
//'   - \mjseqn{\gamma} is `param_baseflow_map_gamma`
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
//' @details
//' # **_ThreshPow** 
//'
//' \if{html}{\figure{mdl_percola_thp.svg}}
//' \if{latex}{\figure{mdl_percola_thp.pdf}{options: width=140mm}}
//' 
//' based on the `_MaxPow` and add the one threshold \mjseqn{\phi_b}: 
//' \mjsdeqn{F_{prcl} = 0, \quad \frac{W_{soil}}{C_{soil}} < \phi_b}
//' \mjsdeqn{F_{prcl} = M_{prcl} \left(\frac{\frac{W_{soil}}{C_{soil}} - \phi_b}{1-\phi_b} \right)^\gamma, \quad \frac{W_{soil}}{C_{soil}} \geq \phi_b}
//' where
//'   - \mjseqn{\phi_b} is `param_percola_thp_thresh`
//'   - \mjseqn{\gamma} is `param_percola_thp_gamma`
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
//' @details
//' # **_Arno** \insertCite{baseflow_Arno_1991,VIC2_Liang_1994}{EDCHM} 
//'
//' \if{html}{\figure{mdl_percola_arn.svg}}
//' \if{latex}{\figure{mdl_percola_arn.pdf}{options: width=140mm}}
//' 
//' has also in two cases divided by a threshold water content \mjseqn{\phi_b}:
//' (*This method is actually not the original method, but an analogy with `baseflow_Arno`*) 
//' \mjsdeqn{F_{prcl} = k M_{prcl} \frac{W_{soil}}{C_{soil}}, \quad \frac{W_{soil}}{C_{soil}} < \phi_b}
//' \mjsdeqn{F_{prcl} = k M_{prcl} \frac{W_{soil}}{C_{soil}} + (1-k) M_{prcl} \left(\frac{W_{soil} - W_s}{C_{soil} - W_s} \right)^2, \quad \frac{W_{soil}}{C_{soil}} \geq \phi_b}
//' \mjsdeqn{W_s = k C_{soil}}
//' where
//'   - \mjseqn{\phi_b} is `param_percola_arn_thresh`
//'   - \mjseqn{k} is `param_percola_arn_k`
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
//' @details
//' # **_BevenWood** \insertCite{percola_BevenWood_1983,TOPMODEL_Beven_1995}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_percola_bew.svg}}
//' \if{latex}{\figure{mdl_percola_bew.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{k =  \frac{W_{soil}}{C_{soil} - W_{soil}} \quad {\rm and} \quad k \leq 1}
//' \mjsdeqn{F_{prcl} = k M_{prcl}}
//' where
//'   - \mjseqn{k_{fc}} is `soil_fieldCapacityPerc_1`
//'   - \mjseqn{\gamma} is `param_percola_sup_gamma`
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



//' @rdname percola
//' @details
//' # **_SupplyPow** \insertCite{GR4J_Perrin_2003}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_percola_sup.svg}}
//' \if{latex}{\figure{mdl_percola_sup.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{k^* = k \left(\frac{W_{soil}}{C_{soil}} \right)^\gamma}
//' where
//'   - \mjseqn{k} is `param_percola_sup_k`
//'   - \mjseqn{\gamma} is `param_percola_sup_gamma`
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
//' @details
//' # **_SupplyRatio** \insertCite{GR4J_Perrin_2003}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_percola_sur.svg}}
//' \if{latex}{\figure{mdl_percola_sur.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{k^* = k}
//' where
//'   - \mjseqn{k} is `param_percola_sur_k`
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
