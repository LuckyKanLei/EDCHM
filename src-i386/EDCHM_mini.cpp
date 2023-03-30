#include "EDCHM_mini.h"
// [[Rcpp::interfaces(r, cpp)]]


//' modells build with EDCHM modulas
//' @name modells
//' @description some example models with EDCHM modulas
//' @inheritParams all_param
//' @inheritParams all_vari
//' @references
//' \insertAllCited{}
//' @return 
//' - `EDCHM_xxxx`: stream flow in mm/TS
//' - `EDCHM_xxxx_full`: list of variablen
//' @details
//' # **EDCHM_mini**: 
//' A model based on mini-structure with only six process:
//' - evatransActual (evatransSoil)
//' - infilt
//' - percola
//' - baseflow
//' - confluenIUH (confluenLand)
//' - confluenIUH (confluenGround)
//' @param confluenLand_responseTime_TS,confluenGround_responseTime_TS,param_confluenLand_kel_k, parameters in modells
//' @param n_time,n_spat number of time step and spatial unit
//' @export
// [[Rcpp::export]]
NumericMatrix EDCHM_mini(
int n_time, 
int n_spat,
NumericMatrix atmos_potentialEvatrans_mm, 
NumericMatrix atmos_precipitation_mm, 
NumericVector ground_capacity_mm, 
NumericVector ground_water_mm, 
NumericVector land_impermeableFrac_1, 
NumericVector soil_capacity_mm, 
NumericVector soil_potentialPercola_mm, 
NumericVector soil_water_mm, 
NumericVector confluenLand_responseTime_TS, 
NumericVector confluenGround_responseTime_TS, 
NumericVector param_baseflow_grf_gamma, 
NumericVector param_confluenLand_kel_k, 
NumericVector param_evatrans_ubc_gamma, 
NumericVector param_infilt_ubc_P0AGEN, 
NumericVector param_percola_arn_k, 
NumericVector param_percola_arn_thresh
)
{

NumericVector land_water_mm, soil_evatrans_mm, soil_infilt_mm, soil_percolation_mm, confluenLand_iuh_1, confluenGround_iuh_1;
NumericMatrix land_runoff_mm(n_time, n_spat), ground_baseflow_mm(n_time, n_spat), confluen_streamflow_mm(n_time, n_spat);
for (int i= 0; i < n_time; i++) {

soil_evatrans_mm = evatransActual_UBC(atmos_potentialEvatrans_mm(i, _), soil_water_mm, soil_capacity_mm, param_evatrans_ubc_gamma);
soil_water_mm += - soil_evatrans_mm;
land_water_mm = atmos_precipitation_mm(i, _);

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

}
for (int j= 0; j < n_spat; j++) {
confluenLand_iuh_1 = confluenIUH_Kelly(confluenLand_responseTime_TS(j), param_confluenLand_kel_k(j));
confluenGround_iuh_1 = confluenIUH_GR4J1(confluenGround_responseTime_TS(j));

confluen_streamflow_mm(_, j) = confluen_IUH2S(land_runoff_mm(_, j), ground_baseflow_mm(_, j), confluenLand_iuh_1, confluenGround_iuh_1);
}
return confluen_streamflow_mm;
}
