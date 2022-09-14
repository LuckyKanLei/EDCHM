#include "00utilis.h"
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
//' @inheritParams all_vari
//' @return IUH (list of num vector) 
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_GR4J1(
    double confluen_responseTime_TS
)
{
  double t_max = ceil(confluen_responseTime_TS);
  IntegerVector seq_t = seq(1, t_max);
  NumericVector seq_t2 = as<NumericVector>(seq_t);
  NumericVector SH_1 = pow(( seq_t2/ confluen_responseTime_TS), 2.5);
  SH_1(t_max - 1) = 1;
  SH_1[Range(1, t_max - 1)] = diff(SH_1);
  return SH_1;
}



//' @rdname confluenIUH
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_GR4J2(
    double confluen_responseTime_TS
)
{
  double t_max_1 = ceil(confluen_responseTime_TS);
  double t_max_2 = ceil(2 * confluen_responseTime_TS);
  IntegerVector seq_t1 = seq(1, t_max_1 - 1);
  NumericVector seq_t12 = as<NumericVector>(seq_t1);
  IntegerVector seq_t2 = seq(t_max_1, (t_max_2 - 1));
  NumericVector seq_t22 = as<NumericVector>(seq_t2);
  
  NumericVector SH_2_1 = .5 * pow((seq_t12 / confluen_responseTime_TS),2.5);
  NumericVector SH_2_2 = 1 - .5 * pow((2 - seq_t22 / confluen_responseTime_TS),2.5);
  NumericVector SH_2(t_max_2, 1);
  SH_2[Range(0, t_max_1 - 2)] = SH_2_1;
  SH_2[Range(t_max_1 - 1, t_max_2 - 2)] = SH_2_2;
  SH_2[Range(1, t_max_2 - 1)] = diff(SH_2);
  
  return SH_2;
}

//' @rdname confluenIUH
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_Clark(
    double confluen_responseTime_TS
)
{
  double t_max = ceil(- confluen_responseTime_TS * log(confluen_responseTime_TS * 0.005));
  IntegerVector seq_t = seq(1, 20 * t_max);
  NumericVector seq_t2 = as<NumericVector>(seq_t) / 20.0;
  NumericVector iuh_ = 1 / confluen_responseTime_TS * exp(- seq_t2 / confluen_responseTime_TS);
  NumericMatrix mat_iuh = NumericMatrix(20, t_max, iuh_.begin());
  NumericVector vct_iuh = colMeans(mat_iuh);
  return vct_iuh / sum(vct_iuh);
}

//' @rdname confluenIUH
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_Kelly(
    double confluen_responseTime_TS,
    double confluen_concentratTime_TS
)
{
  double num_temp_tc2 = (confluen_concentratTime_TS * confluen_concentratTime_TS);
  double num_temp_12_34 = 4 * confluen_responseTime_TS  / num_temp_tc2 * 
    (1 - 2 * exp(confluen_concentratTime_TS / confluen_responseTime_TS * 0.5));
  double num_temp_12_35 = 4 * confluen_responseTime_TS  / num_temp_tc2 * 
    (1 - 2 * exp(confluen_concentratTime_TS / confluen_responseTime_TS * 0.5) + exp(confluen_concentratTime_TS / confluen_responseTime_TS));
  double t_max = ceil(std::max(confluen_concentratTime_TS, - confluen_responseTime_TS * log(0.002 / num_temp_12_35)));
  NumericVector iuh_1, iuh_2, iuh_3, iuh_, vct_iuh, temp_etK;
  IntegerVector seq_t = seq(1, 20 * t_max);
  NumericVector seq_t2 = as<NumericVector>(seq_t) / 20.0;
  temp_etK = exp(- seq_t2 / confluen_responseTime_TS);
  iuh_1 = 4 / num_temp_tc2 * (seq_t2 + confluen_responseTime_TS * (temp_etK - 1));
  iuh_2 = num_temp_12_34 * temp_etK - 4 / num_temp_tc2 * (seq_t2 - confluen_responseTime_TS - confluen_concentratTime_TS);
  iuh_3 = num_temp_12_35 * temp_etK;
  iuh_ = ifelse(seq_t2 > confluen_concentratTime_TS * 0.5, iuh_2, iuh_1);
  iuh_ = ifelse(seq_t2 > confluen_concentratTime_TS, iuh_3, iuh_);
  NumericMatrix mat_iuh = NumericMatrix(20, t_max, iuh_.begin());
  vct_iuh = colMeans(mat_iuh);
  return vct_iuh / sum(vct_iuh);
}


//' @rdname confluenIUH
//' @param param_confluen_nas_n parameters for[confluenIUH_Nash()]
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_Nash(
    double confluen_responseTime_TS,
    double param_confluen_nas_n
)
{
  double t_max = ceil(std::max(4.0, param_confluen_nas_n) * 3 * confluen_responseTime_TS);
  NumericVector iuh_, vct_iuh;
  IntegerVector seq_t = seq(1, 20 * t_max);
  NumericVector seq_t2 = as<NumericVector>(seq_t) / 20.0;
  iuh_ = pow(seq_t2 / confluen_responseTime_TS, param_confluen_nas_n - 1) * exp(- seq_t2 / confluen_responseTime_TS) / 
    confluen_responseTime_TS / tgamma(param_confluen_nas_n);
  NumericMatrix mat_iuh = NumericMatrix(20, t_max, iuh_.begin());
  vct_iuh = colMeans(mat_iuh);
  return vct_iuh / sum(vct_iuh);
}

//' @rdname confluenIUH
//' @param param_confluen_nak_b,param_confluen_nak_n parameters for[confluenIUH_NashKumar()]
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_NashKumar(
    double param_confluen_nak_b,
    double param_confluen_nak_n
)
{
  NumericVector iuh_, vct_iuh;
  double confluen_responseTime_TS = param_confluen_nak_b * 1 / (param_confluen_nak_n - 1);
  double t_max = ceil(std::max(4.0, param_confluen_nak_n) * 3 * confluen_responseTime_TS);
  IntegerVector seq_t = seq(1, 20 * t_max);
  NumericVector seq_t2 = as<NumericVector>(seq_t) / 20.0;
  iuh_ = pow(seq_t2 / confluen_responseTime_TS, param_confluen_nak_n - 1) * exp(- seq_t2 / confluen_responseTime_TS) / 
    confluen_responseTime_TS / tgamma(param_confluen_nak_n);
  NumericMatrix mat_iuh = NumericMatrix(20, t_max, iuh_.begin());
  vct_iuh = colMeans(mat_iuh);
  return vct_iuh / sum(vct_iuh);
}
