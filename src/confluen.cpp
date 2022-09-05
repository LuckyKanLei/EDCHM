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


//' create **IUH** (Instant Unit Graphy)
//' @name confluenIUH
//' @param confluen_resposeTime_TS (TS) respose time in every routeline
//' @return IUH (list of num vector) 
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_GR4J1(
    double confluen_resposeTime_TS
)
{
  double t_max = ceil(confluen_resposeTime_TS);
  IntegerVector seq_t = seq(1, t_max);
  NumericVector seq_t2 = as<NumericVector>(seq_t);
  NumericVector SH_1 = pow(( seq_t2/ confluen_resposeTime_TS), 2.5);
  SH_1(t_max - 1) = 1;
  SH_1[Range(1, t_max - 1)] = diff(SH_1);
  return SH_1;
}



//' @rdname confluenIUH
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_GR4J2(
    double confluen_resposeTime_TS
)
{
  double t_max_1 = ceil(confluen_resposeTime_TS);
  double t_max_2 = ceil(2 * confluen_resposeTime_TS);
  IntegerVector seq_t1 = seq(1, t_max_1 - 1);
  NumericVector seq_t12 = as<NumericVector>(seq_t1);
  IntegerVector seq_t2 = seq(t_max_1, (t_max_2 - 1));
  NumericVector seq_t22 = as<NumericVector>(seq_t2);
  
  NumericVector SH_2_1 = .5 * pow((seq_t12 / confluen_resposeTime_TS),2.5);
  NumericVector SH_2_2 = 1 - .5 * pow((2 - seq_t22 / confluen_resposeTime_TS),2.5);
  NumericVector SH_2(t_max_2, 1);
  SH_2[Range(0, t_max_1 - 2)] = SH_2_1;
  SH_2[Range(t_max_1 - 1, t_max_2 - 2)] = SH_2_2;
  SH_2[Range(1, t_max_2 - 1)] = diff(SH_2);
  
  return SH_2;
}
