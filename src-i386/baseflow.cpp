#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]



//' **baseflow**
//' @name baseflow
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' 
//' In hydrological modeling, baseflow refers to the flow of water in rivers and streams that is sustained by the release of water from the groundwater.
//' Or baseflow refers to the flow of water from an aquifer or deeper soil horizon to surface water, typically due to a head gradient between fully saturated soil and stream  \insertCite{Raven_Manual_35}{EDCHM}. 
//' It may be considered the sum of the contribution of deep groundwater exchange with a river and delayed storage  \insertCite{Raven_Manual_35}{EDCHM}.
//' 
//' It is always calculated (only) by the water in the ground layer \mjseqn{W_{grnd}}, which can also be treated as part of \mjseqn{W_{grnd}}. 
//' However, the impact of other RUs (response units) on the route to the river will be ignored.
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{F_{base} = f_{baseflow}(D_{grnd})}
//' 
//' 
//' to:
//' 
//' \mjsdeqn{F_{base} = f_{baseflow}(W_{grnd}, C_{grnd}, M_{base}, ...)}
//' \mjsdeqn{F_{base} = k^* W_{grnd} \quad {\rm or} \quad F_{base} = k^* M_{base}}
//' \mjsdeqn{0 \leq k^* \leq 1}
//' 
//' 
//' where
//' - \mjseqn{W_{grnd}} is `ground_water_mm`
//' - \mjseqn{M_{base}} is `ground_potentialBaseflow_mm`
//' - \mjseqn{C_{grnd}} is `ground_capacity_mm`, but not all the methods need the \mjseqn{C_{grnd}}
//' - \mjseqn{k^*} is estimated ratio
//' 
//' The output density distribution from 7 methods:
//'
//' \if{html}{\figure{mdl_baseflow.svg}}
//' \if{latex}{\figure{mdl_baseflow.pdf}{options: width=140mm}}
//' @references
//' \insertAllCited{}
//' @return ground_baseflow_mm (mm/m2/TS) 
//' @details
//' # **_GR4J** \insertCite{GR4J_Perrin_2003}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_baseflow_gr4.svg}}
//' \if{latex}{\figure{mdl_baseflow_gr4.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{base} = k^* W_{grnd}}
//' \mjsdeqn{k^* = 1 - \left[ 1 + \left(\frac{W_{grnd}}{C_{grnd}} \right)^4 \right]^{-1/4}}
//' where
//'   - \mjseqn{k^*} is estimated ratio
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
//' @details
//' # **_GR4Jfix** \insertCite{GR4J_Perrin_2003}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_baseflow_grf.svg}}
//' \if{latex}{\figure{mdl_baseflow_grf.pdf}{options: width=140mm}}
//' 
//' This method based on `_GR4J` use a new parameter to replace the numer 4:
//' \mjsdeqn{F_{base} = k^* W_{grnd}}
//' \mjsdeqn{k^* = 1 - \left[ 1 + \left(\frac{W_{grnd}}{C_{grnd}} \right)^\gamma \right]^{-1/\gamma}}
//' where
//'   - \mjseqn{\gamma} is `param_baseflow_grf_gamma`
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
//' @details
//' # **_SupplyRatio**: 
//'
//' \if{html}{\figure{mdl_baseflow_sur.svg}}
//' \if{latex}{\figure{mdl_baseflow_sur.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{base} = k W_{grnd}}
//' where
//'   - \mjseqn{k} is `param_baseflow_sur_k`
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
//' @details
//' # **_SupplyPow**: 
//'
//' \if{html}{\figure{mdl_baseflow_sup.svg}}
//' \if{latex}{\figure{mdl_baseflow_sup.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{base} = k(W_{grnd})^\gamma}
//' where
//'   - \mjseqn{k} is `param_baseflow_sup_k`
//'   - \mjseqn{\gamma} is `param_baseflow_sup_gamma`
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
//' @details
//' # **_MaxPow**: 
//'
//' \if{html}{\figure{mdl_baseflow_map.svg}}
//' \if{latex}{\figure{mdl_baseflow_map.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{base} = M_{base} \left(\frac{W_{grnd}}{C_{grnd}} \right)^\gamma}
//' where
//'   - \mjseqn{M_{base}} is `ground_potentialBaseflow_mm`
//'   - \mjseqn{\gamma} is `param_baseflow_map_gamma`
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
//' @details
//' # **_ThreshPow** 
//'
//' \if{html}{\figure{mdl_baseflow_thp.svg}}
//' \if{latex}{\figure{mdl_baseflow_thp.pdf}{options: width=140mm}}
//' 
//' This method based on the `_MaxPow` and add the one threshold \mjseqn{\phi_b}: 
//' \mjsdeqn{F_{base} = 0, \quad \frac{W_{grnd}}{C_{grnd}} < \phi_b}
//' \mjsdeqn{F_{base} = M_{base} \left(\frac{\frac{W_{grnd}}{C_{grnd}} - \phi_b}{1-\phi_b} \right)^\gamma, \quad \frac{W_{grnd}}{C_{grnd}} \geq \phi_b}
//' where
//'   - \mjseqn{\phi_b} is `param_baseflow_thp_thresh`
//'   - \mjseqn{\gamma} is `param_baseflow_thp_gamma`
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
//' @details
//' # **_Arno** \insertCite{baseflow_Arno_1991,VIC2_Liang_1994}{EDCHM}:
//'
//' \if{html}{\figure{mdl_baseflow_arn.svg}}
//' \if{latex}{\figure{mdl_baseflow_arn.pdf}{options: width=140mm}}
//' 
//' This method has also in two cases divided by a threshold water content \mjseqn{\phi_b}: 
//' \mjsdeqn{F_{base} = k M_{base} \frac{W_{grnd}}{C_{grnd}}, \quad \frac{W_{grnd}}{C_{grnd}} < \phi_b}
//' \mjsdeqn{F_{base} = k M_{base} \frac{W_{grnd}}{C_{grnd}} + (1-k) M_{base} \left(\frac{W_{grnd} - W_s}{C_{grnd} - W_s} \right)^2, \quad \frac{W_{grnd}}{C_{grnd}} \geq \phi_b}
//' \mjsdeqn{W_s = k C_{grnd}}
//' where
//'   - \mjseqn{\phi_b} is `param_baseflow_arn_thresh`
//'   - \mjseqn{k} is `param_baseflow_arn_k`
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
