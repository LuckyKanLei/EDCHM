#include "EDCHM_snow.h"
// [[Rcpp::interfaces(r, cpp)]]


//' @name modells
//' @details
//' # **EDCHM_snow_full**: 
//' With full Output
//' @export
// [[Rcpp::export]]
List EDCHM_snow_full(
int n_time, 
int n_spat,
int time_step_h,
NumericMatrix atmos_potentialEvatrans_mm, 
NumericMatrix atmos_precipitation_mm, 
NumericMatrix atmos_temperature_Cel, 
NumericVector ground_capacity_mm, 
NumericVector ground_water_mm, 
NumericVector land_impermeableFrac_1, 
NumericVector snow_ice_mm, 
NumericVector soil_capacity_mm, 
NumericVector soil_potentialPercola_mm, 
NumericVector soil_water_mm, 
NumericVector confluenLand_responseTime_TS, 
NumericVector confluenGround_responseTime_TS, 
NumericVector param_atmos_thr_Ts, 
NumericVector param_baseflow_grf_gamma, 
NumericVector param_confluenLand_kel_k, 
NumericVector param_evatrans_ubc_gamma, 
NumericVector param_infilt_ubc_P0AGEN, 
NumericVector param_percola_arn_k, 
NumericVector param_percola_arn_thresh, 
NumericVector param_snow_fac_f, 
NumericVector param_snow_fac_Tmelt
)
{

NumericVector land_water_mm, atmos_snow_mm, soil_evatrans_mm, snow_melt_mm, soil_infilt_mm, soil_percolation_mm, confluenLand_iuh_1, confluenGround_iuh_1;
NumericMatrix land_runoff_mm(n_time, n_spat), ground_baseflow_mm(n_time, n_spat), confluen_streamflow_mm(n_time, n_spat);
NumericMatrix out_evatrans(n_time, n_spat), out_soilwater(n_time, n_spat), out_groundwater(n_time, n_spat), out_snowice(n_time, n_spat), out_snowmelt(n_time, n_spat);
for (int i= 0; i < n_time; i++) {

atmos_snow_mm = atmosSnow_ThresholdT(atmos_precipitation_mm(i, _), atmos_temperature_Cel(i, _), param_atmos_thr_Ts);
atmos_precipitation_mm(i, _) = atmos_precipitation_mm(i, _) - atmos_snow_mm;

soil_evatrans_mm = evatransActual_UBC(atmos_potentialEvatrans_mm(i, _), soil_water_mm, soil_capacity_mm, param_evatrans_ubc_gamma);
soil_water_mm += - soil_evatrans_mm;
land_water_mm = atmos_precipitation_mm(i, _);

snow_melt_mm = snowMelt_Factor(time_step_h, snow_ice_mm, atmos_temperature_Cel(i, _), param_snow_fac_f, param_snow_fac_Tmelt);
land_water_mm += snow_melt_mm;
snow_ice_mm += -snow_melt_mm;
snow_ice_mm += atmos_snow_mm;

soil_infilt_mm = infilt_UBC(land_water_mm, land_impermeableFrac_1, soil_water_mm, soil_capacity_mm, param_infilt_ubc_P0AGEN);
soil_water_mm += soil_infilt_mm;
land_runoff_mm(i, _) = land_water_mm - soil_infilt_mm;

soil_percolation_mm = percola_Arno(soil_water_mm, soil_capacity_mm, soil_potentialPercola_mm, param_percola_arn_thresh, param_percola_arn_k);
ground_water_mm += soil_percolation_mm;
soil_water_mm += - soil_percolation_mm;

NumericVector baseflow_temp = ifelse(ground_water_mm < ground_capacity_mm, 0, ground_water_mm - ground_capacity_mm);

ground_water_mm = ifelse(ground_water_mm < ground_capacity_mm,ground_water_mm, ground_capacity_mm);
ground_baseflow_mm(i, _) = baseflow_GR4Jfix(ground_water_mm, ground_capacity_mm, param_baseflow_grf_gamma);
ground_water_mm += - ground_baseflow_mm(i, _);
ground_baseflow_mm(i, _) = ground_baseflow_mm(i, _) + baseflow_temp;

    out_evatrans(i, _) = soil_evatrans_mm;
    out_soilwater(i, _) = soil_water_mm;
    out_groundwater(i, _) = ground_water_mm;
    out_snowice(i, _) = snow_ice_mm;
    out_snowmelt(i, _) = snow_melt_mm;
    
}
for (int j= 0; j < n_spat; j++) {
confluenLand_iuh_1 = confluenIUH_Kelly(confluenLand_responseTime_TS(j), param_confluenLand_kel_k(j));
confluenGround_iuh_1 = confluenIUH_GR4J1(confluenGround_responseTime_TS(j));

confluen_streamflow_mm(_, j) = confluen_IUH2S(land_runoff_mm(_, j), ground_baseflow_mm(_, j), confluenLand_iuh_1, confluenGround_iuh_1);
  }
  return List::create(
    _["evatrans_mm"] = out_evatrans,
    _["soilwater_mm"] = out_soilwater,
    _["groundwater_mm"] = out_groundwater,
    _["snowice_mm"] = out_snowice,
    _["snowmelt_mm"] = out_snowmelt,
    _["runoff_mm"] = land_runoff_mm,
    _["baseflow_mm"] = ground_baseflow_mm,
    _["streamflow_mm"] = confluen_streamflow_mm
  );
}
