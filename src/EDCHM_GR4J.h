// Defines a header file containing function for EDCHM_GR4J/
#ifndef EDCHM_GR4J_H
#define EDCHM_GR4J_H

#include <Rcpp.h>
using namespace Rcpp;

double sum_product(NumericVector lhs, NumericVector rhs);
void resetVector(Rcpp::NumericVector& x);

NumericVector infilt_GR4J(
    NumericVector land_water_mm,
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm
);

NumericVector evatransActual_GR4J(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm
);

NumericVector percola_GR4J(
    NumericVector soil_water_mm,
    NumericVector soil_capacity_mm
);
NumericVector baseflow_GR4J(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm
);
NumericVector lateral_GR4J(
    NumericVector ground_water_mm,
    NumericVector ground_capacity_mm,
    NumericVector ground_potentialLateral_mm
);
NumericVector confluen_IUH2S(
    NumericVector land_runoff_mm,
    NumericVector ground_baseflow_mm, 
    NumericVector confluen_iuhLand_1,
    NumericVector confluen_iuhGround_1
);

NumericVector confluenIUH_GR4J2(
    double confluen_responseTime_TS
);

NumericVector confluenIUH_GR4J1(
    double confluen_responseTime_TS
);

NumericVector vecpow(NumericVector base, NumericVector exp);
NumericVector vecpow10(NumericVector exp);

#endif