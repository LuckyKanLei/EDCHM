#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]



//' **capilarise**
//' @name capirise
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' 
//' In hydrological modeling, capillary rise refers to the process by which water is drawn upward from groundwater (table) through the soil due to the force of capillary action.
//' In conceptual watershed models, the capillary rise term often refers to a process that moves water from lower to higher soil water stores, 
//' which may also implicitly include lateral groundwater flow processes in a sloping domain  \insertCite{Raven_Manual_35}{EDCHM}.
//' 
//' It can be calculated by the water in the ground layer \mjseqn{W_{grnd}}, which can also be treated as part of \mjseqn{W_{grnd}}. 
//' There are not many methods to describe this process. Most HMs ignore this process, 
//' perhaps because it is not significant in most situations, or because the process of percolation can deal with this process at the same time.
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{F_{capi} = f_{capirise}(D_{grnd}, D_{soil})}
//' 
//' 
//' to:
//' 
//' \mjsdeqn{F_{capi} = f_{capirise}(W_{grnd}, W_{soil}, C_{soil}, ...)}
//' \mjsdeqn{F_{capi} \leq W_{grnd}}
//' \mjsdeqn{F_{capi} \leq C_{soil} - W_{soil}}
//' 
//' 
//' where
//' - \mjseqn{F_{capi}} is `ground_capirise_mm`
//' - \mjseqn{W_{grnd}} is `ground_water_mm`
//' - \mjseqn{W_{soil}} is `water_soil_mm`
//' - \mjseqn{C_{soil}} is `capacity_soil_mm`
//' 
//' The output density distribution from 4 methods:
//'
//' \if{html}{\figure{mdl_capirise.svg}}
//' \if{latex}{\figure{mdl_capirise.pdf}{options: width=140mm}}
//' @references
//' \insertAllCited{}
//' @return ground_capirise_mm (mm/m2/TS) capillary rise
//' 
//' @details
//' # **_HBV** \insertCite{HBV_Lindstrom_1997}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_capirise_hbv.svg}}
//' \if{latex}{\figure{mdl_capirise_hbv.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{capi} = M_{capi} \left( 1 - \frac{W_{soil}}{C_{soil}} \right)}
//' where
//'   - \mjseqn{M_{capi}} is `soil_potentialCapirise_mm`
//' @export
// [[Rcpp::export]]
NumericVector capirise_HBV(
    NumericVector ground_water_mm, 
    NumericVector soil_water_mm ,
    NumericVector soil_capacity_mm, 
    NumericVector soil_potentialCapirise_mm
)
{
  NumericVector soil_diff_mm, capirise_mm, limit_mm;
  soil_diff_mm = soil_capacity_mm  - soil_water_mm;
  soil_diff_mm = ifelse(soil_diff_mm < 0, 0, soil_diff_mm);
  capirise_mm = soil_potentialCapirise_mm * (soil_diff_mm / soil_capacity_mm);
  
  limit_mm = ifelse(soil_diff_mm > ground_water_mm, ground_water_mm, soil_diff_mm) ;
  return ifelse(capirise_mm > limit_mm, limit_mm, capirise_mm) ;
}



//' @rdname capirise
//' @details
//' # **_HBVfix** \insertCite{HBV_Lindstrom_1997}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_capirise_hbf.svg}}
//' \if{latex}{\figure{mdl_capirise_hbf.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{capi} = M_{capi} \left( 1 - \frac{W_{soil}}{k_{fc}C_{soil}} \right), \quad W_{soil} < k_{fc}C_{soil}}
//' where
//'   - \mjseqn{k_{fc}} is `soil_fieldCapacityPerc_1`
//' @export
// [[Rcpp::export]]
NumericVector capirise_HBVfix(
    NumericVector ground_water_mm, 
    NumericVector soil_water_mm ,
    NumericVector soil_capacity_mm, 
    NumericVector soil_fieldCapacityPerc_1,
    NumericVector soil_potentialCapirise_mm
)
{
  NumericVector soil_diff_mm, capirise_mm, limit_mm;
  soil_diff_mm = soil_capacity_mm * (1 - soil_fieldCapacityPerc_1) - soil_water_mm;
  soil_diff_mm = ifelse(soil_diff_mm < 0, 0, soil_diff_mm);
  capirise_mm = soil_potentialCapirise_mm * (soil_diff_mm / soil_capacity_mm);
  
  limit_mm = ifelse(soil_diff_mm > ground_water_mm, ground_water_mm, soil_diff_mm) ;
  return ifelse(capirise_mm > limit_mm, limit_mm, capirise_mm) ;
}

//' @rdname capirise
//' @details
//' # **_AcceptRatio**: 
//'
//' \if{html}{\figure{mdl_capirise_acr.svg}}
//' \if{latex}{\figure{mdl_capirise_acr.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{capi} = k \left( W_{soil} - k_{fc}C_{soil} \right), \quad W_{soil} < k_{fc}C_{soil}}
//' where
//'   - \mjseqn{k} is `param_capirise_acr_k`
//'   - \mjseqn{k_{fc}} is `soil_fieldCapacityPerc_1`
//' @param param_capirise_acr_k <0.01, 1> coefficient parameter [capirise_AcceptRatio()]
//' @export
// [[Rcpp::export]]
NumericVector capirise_AcceptRatio(
    NumericVector ground_water_mm, 
    NumericVector soil_water_mm ,
    NumericVector soil_capacity_mm, 
    NumericVector soil_fieldCapacityPerc_1,
    NumericVector param_capirise_acr_k
)
{
  NumericVector soil_diff_mm, capirise_mm, limit_mm;
  soil_diff_mm = soil_capacity_mm * (1 - soil_fieldCapacityPerc_1) - soil_water_mm;
  soil_diff_mm = ifelse(soil_diff_mm < 0, 0, soil_diff_mm);
  capirise_mm = soil_diff_mm * param_capirise_acr_k;
  
  limit_mm = ifelse(soil_diff_mm > ground_water_mm, ground_water_mm, soil_diff_mm) ;
  return ifelse(capirise_mm > limit_mm, limit_mm, capirise_mm) ;
}



//' @rdname capirise
//' @details
//' # **_AcceptRatio**: 
//'
//' \if{html}{\figure{mdl_capirise_acp.svg}}
//' \if{latex}{\figure{mdl_capirise_acp.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{capi} = k \left( W_{soil} - k_{fc}C_{soil} \right)^\gamma, \quad W_{soil} < k_{fc}C_{soil}}
//' where
//'   - \mjseqn{k} is `param_capirise_acp_k`
//'   - \mjseqn{\gamma} is `param_capirise_acp_gamma`
//' @param param_capirise_acp_k <0.01, 1> coefficient parameter for [capirise_AcceptPow()]
//' @param param_capirise_acp_gamma <0.01, 1> exponential parameter for [capirise_AcceptPow()]
//' @export
// [[Rcpp::export]]
NumericVector capirise_AcceptPow(
    NumericVector ground_water_mm, 
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector soil_fieldCapacityPerc_1,
    NumericVector param_capirise_acp_k,
    NumericVector param_capirise_acp_gamma
)
{
  NumericVector capirise_mm, k_, soil_diff_mm, limit_mm;
  soil_diff_mm = soil_capacity_mm * (1 - soil_fieldCapacityPerc_1) - soil_water_mm;
  soil_diff_mm = ifelse(soil_diff_mm < 0, 0, soil_diff_mm);
  
  k_ = param_capirise_acp_k * vecpow((soil_diff_mm / (soil_capacity_mm * (1 - soil_fieldCapacityPerc_1))), param_capirise_acp_gamma);
  capirise_mm = k_ * soil_diff_mm;
  capirise_mm = ifelse(capirise_mm < 0, 0, capirise_mm);
  
  limit_mm = ifelse(soil_diff_mm > ground_water_mm, ground_water_mm, soil_diff_mm) ;
  return ifelse(capirise_mm > limit_mm, limit_mm, capirise_mm) ;
}



