#include <Rcpp.h>
#include "00utilis.h"
using namespace Rcpp;
// [[Rcpp::interfaces(r, cpp)]]


//' **confluence**
//' @description Routing methods with 
//' - `IUH`: IUH (Instant Unit Hydrograph) with one watercourse, 
//' - `IUH2S`; IUH with tweo watersource, those have the different IUH-vector, 
//' @inheritParams all_vari
//' @name confluen
//' @return confluenced water (mm/m2)
//' @export
// [[Rcpp::export]]
NumericVector confluen_IUH(
    NumericVector confluen_inputWater_mm, 
    NumericVector confluen_iuh_1
)
{
  
  int n_iuh = confluen_iuh_1.size(), n_time = confluen_inputWater_mm.size();
  NumericVector confluen_outputWater_mm (n_time);
  for (int i = 0; i < n_iuh; i++) {
    for (int j = 0; j <= i; j++) {
      confluen_outputWater_mm[i] += confluen_inputWater_mm[i-j] * confluen_iuh_1[j];
    }
  }
  for (int i = n_iuh; i < n_time; i++) {
    for (int j = 0; j < n_iuh; j++) {
      confluen_outputWater_mm[i] += confluen_inputWater_mm[i-j] * confluen_iuh_1[j];
    }
  }
  
  return confluen_outputWater_mm;
  
}

//' @rdname confluen
//' @export
// [[Rcpp::export]]
NumericVector confluen_IUH2S(
    NumericVector land_runoff_mm,
    NumericVector ground_baseflow_mm, 
    NumericVector confluen_iuhLand_1,
    NumericVector confluen_iuhGround_1
)
{
  NumericVector confluen_runoff_mm (land_runoff_mm.size()), confluen_baseflow_mm (ground_baseflow_mm.size());
  confluen_runoff_mm = confluen_IUH(
    land_runoff_mm, 
    confluen_iuhLand_1
  );
  confluen_baseflow_mm = confluen_IUH(
    ground_baseflow_mm, 
    confluen_iuhGround_1
  );
  

  return confluen_runoff_mm + confluen_baseflow_mm;
  
}
