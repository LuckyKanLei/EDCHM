// Defines a header file containing function for EDCHM_snow/
#ifndef EDCHM_SNOW_H
#define EDCHM_SNOW_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector atmosSnow_ThresholdT(
    NumericVector atmos_precipitation_mm, 
    NumericVector atmos_temperature_Cel, 
    NumericVector param_atmos_thr_Ts
);
NumericVector evatransActual_UBC(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_evatrans_ubc_gamma
);
NumericVector snowMelt_Factor(
    NumericVector snow_ice_mm,
    NumericVector atmos_temperature_Cel,
    NumericVector param_snow_fac_f,
    NumericVector param_snow_fac_Tmelt
);

NumericVector infilt_UBC(
    NumericVector land_water_mm, 
    NumericVector land_impermeableFrac_1, 
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm, 
    NumericVector param_infilt_ubc_P0AGEN
);
NumericVector percola_Arno(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm,
    NumericVector soil_potentialPercola_mm,
    NumericVector param_percola_arn_thresh,
    NumericVector param_percola_arn_k
);
NumericVector baseflow_GR4Jfix(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector param_baseflow_grf_gamma
);
NumericVector confluen_IUH2S(
    NumericVector land_runoff_mm,
    NumericVector ground_baseflow_mm, 
    NumericVector confluen_iuhLand_1,
    NumericVector confluen_iuhGround_1
);

NumericVector confluenIUH_Kelly(
    double confluen_responseTime_TS,
    double param_confluen_kel_k
);

NumericVector confluenIUH_GR4J1(
    double confluen_responseTime_TS
);

NumericVector vecpow(NumericVector base, NumericVector exp);
NumericVector vecpow10(NumericVector exp);

#endif