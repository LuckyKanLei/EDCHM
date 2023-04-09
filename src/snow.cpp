#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]


//' **snow**
//' @name snowMelt
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' 
//' Under the concept of the conceptual HM, the melt of snowpack is always calculated by 
//' the energy availability (the state-variable temperature \mjseqn{T} or flux-variable (nett-) radiation \mjseqn{Rn}) 
//' and the solid water (snow or ice) availability \mjseqn{W_{snow}} of the snowpack. 
//' 
//' Some more complex processes, such as refrozen and residual water, will be ignored. 
//' To simplify the model, the layer snowLy will store only the solid water and will melt it as much as possible when the energy is sufficient.
//' 
//' So we can give the function from:
//' 
//' \mjsdeqn{F_{melt} = f_{snowMelt}(D_{atms}, D_{snow})}
//' 
//' 
//' to:
//' \mjsdeqn{F_{melt}  = f_{snowMelt}(T, ...)}
//' \mjsdeqn{F_{melt} \leq W_{snow} }
//' 
//' where
//'   - \mjseqn{F_{melt}} is `snow_melt_mm`
//'   - \mjseqn{W_{snow}} is `snow_ice_mm`
//'   - \mjseqn{T} is average temperature
//' 
//' Then the different `snowMelt` methods will estimate the maximal snow melt \mjseqn{M_{max}}.
//' 
//' The output density distribution from 2 methods:
//'
//' \if{html}{\figure{mdl_snowMelt.svg}}
//' \if{latex}{\figure{mdl_snowMelt.pdf}{options: width=140mm}}
//' @references
//' \insertAllCited{}
//' @return snow_melt_mm (mm/m2) melted snow
//' 
//' @details
//' # **_Kustas** \insertCite{snow_kustas_1994}{EDCHM}: 
//' 
//'
//' \if{html}{\figure{mdl_snowMelt_kus.svg}}
//' \if{latex}{\figure{mdl_snowMelt_kus.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{melt}  = m_T T + m_E R_n}
//' but due to the temperature is one energy-state-variable, 
//' in order to adjust to subday scale we need to add a new time interval \mjseqn{t_h} from 1 to 24 hour
//' \mjsdeqn{F_{melt}  = m_T T t_h + m_E R_n}
//' where
//'   - \mjseqn{m_T} is `param_snow_kus_fT`
//'   - \mjseqn{m_E} is `param_snow_kus_fE`
//'   - \mjseqn{R_n} is daily net radiation
//' 
//' @param param_snow_kus_fE <0.0005, 0.003> (mm/m2/MJ) snow melt temperature parameter for [snowMelt_Factor()]
//' @param param_snow_kus_fT <0.05, 1> (mm/m2/h/Cel) potential melt volum per Cel per hour parameter for [snowMelt_Factor()]
//' @export
// [[Rcpp::export]]
NumericVector snowMelt_Kustas(
    NumericVector snow_ice_mm,
    NumericVector atmos_temperature_Cel,
    NumericVector atmos_netRadiat_MJ,
    NumericVector param_snow_kus_fE,
    NumericVector param_snow_kus_fT
)
{
  // NumericVector snow_melt_mm = ifelse(atmos_temperature_Cel < 0, 0, atmos_temperature_Cel) * param_snow_kus_fT * time_step_h + param_snow_kus_fE * atmos_netRadiat_MJ;
  NumericVector snow_melt_mm = ifelse(atmos_temperature_Cel < 0, 0, atmos_temperature_Cel) * param_snow_kus_fT * 24 + param_snow_kus_fE * atmos_netRadiat_MJ;
  return ifelse(snow_melt_mm > snow_ice_mm, snow_ice_mm, snow_melt_mm) ;
  
}

//' @rdname snowMelt
//' @details
//' # **_Factor** \insertCite{phyHydro_dingman_2014}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_snowMelt_fac.svg}}
//' \if{latex}{\figure{mdl_snowMelt_fac.pdf}{options: width=140mm}}
//' 
//' \mjsdeqn{F_{melt}  = m_T (T-T_b), \quad T > T_b}
//' where
//'   - \mjseqn{m_T} is `param_snow_fac_f`
//'   - \mjseqn{T_b} is `param_snow_fac_Tmelt`
//' 
//' @param param_snow_fac_Tmelt <0, 3> (Cel) snow melt temperature parameter for [snowMelt_Factor()]
//' @param param_snow_fac_f <0.05, 2> (mm/m2/h/Cel) potential melt volum per Cel per hour parameter for [snowMelt_Factor()]
//' @export
// [[Rcpp::export]]
NumericVector snowMelt_Factor(
    NumericVector snow_ice_mm,
    NumericVector atmos_temperature_Cel,
    NumericVector param_snow_fac_f,
    NumericVector param_snow_fac_Tmelt
)
{
  NumericVector diff_T, snow_melt_mm;
  diff_T = atmos_temperature_Cel - param_snow_fac_Tmelt;
  diff_T = ifelse(diff_T > 0, diff_T, 0);
  
  snow_melt_mm = param_snow_fac_f * 24 * diff_T;
  // snow_melt_mm = param_snow_fac_f * time_step_h * diff_T;
  return ifelse(snow_melt_mm > snow_ice_mm, snow_ice_mm, snow_melt_mm) ;
}

