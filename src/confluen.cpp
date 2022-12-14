#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]


//' **confluence**
//' @description 
//' \loadmathjax
//' 
//' In hydrological modeling, routing (named as [confluen] in EDCHM) refers to the process of simulating the movement of water through a river network or other drainage system. 
//' It allows the model to predict the flow of water in rivers and streams. 
//' In hydrological models, routing is typically performed using mathematical algorithms that account for the physical properties of the river network, 
//' such as its geometry, roughness, and discharge capacity. 
//' The parameters that govern routing, such as flow velocity and channel roughness, 
//' can have a significant impact on the accuracy of the model.
//' 
//' `confluence` is a calculation function that causes water to flow into the gauge point.
//' - `IUH`: IUH (Instant Unit Hydrograph) with one watercourse, 
//' - `IUH2S`: IUH with two water sources, each with a different IUH vector, 
//' - `IUH3S`: IUH with three water sources, each with a different IUH vector.
//' 
//' Under the concept of the conceptual HM, the water flux to the water flow will be calculated using the confluence process. 
//' This process does not calculate the water balance, but rather the time-varying nature of the water flow. 
//' The "Instant Unit Hydrograph" method is the most effective way to deal with time-varying flows. 
//' In the first stage, only [confluenIUH] will be supported.
//' 
//' So we can give the function:
//' 
//' \mjsdeqn{Q = f_{confluen}(F, u)}
//' 
//' 
//' 
//' where
//' - \mjseqn{Q} is stream flow, but still in mm/TS not m3/TS or m3/S
//' - \mjseqn{F} is flux that will into river conflen, e.g.`land_runoff_mm`, `soil_interflow_mm` or `ground_baseflow_mm`
//' - \mjseqn{u} is Instant Unit Hydrograph series
//' 
//' @references
//' \insertAllCited{}
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

//' @rdname confluen
//' @export
// [[Rcpp::export]]
NumericVector confluen_IUH3S(
    NumericVector land_runoff_mm,
    NumericVector soil_interflow_mm, 
    NumericVector ground_baseflow_mm, 
    NumericVector confluen_iuhLand_1,
    NumericVector confluen_iuhSoil_1,
    NumericVector confluen_iuhGround_1
)
{
  NumericVector confluen_runoff_mm (land_runoff_mm.size()), confluen_interflow_mm (soil_interflow_mm.size()), confluen_baseflow_mm (ground_baseflow_mm.size());
  confluen_runoff_mm = confluen_IUH(
    land_runoff_mm, 
    confluen_iuhLand_1
  );
  confluen_interflow_mm = confluen_IUH(
    soil_interflow_mm, 
    confluen_iuhSoil_1
  );
  confluen_baseflow_mm = confluen_IUH(
    ground_baseflow_mm, 
    confluen_iuhGround_1
  );
  
  
  return confluen_runoff_mm + confluen_interflow_mm + confluen_baseflow_mm;
  
}


//' create **IUH** (Instant Unit Hydrograph)
//' @name confluenIUH
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' 
//' The process `confluenIUH` return a series of portions, that means how many flux water will
//' in those moment into the river.
//' The sum of this series will always in 1.
//' 
//' So we can give the function:
//' 
//' \mjsdeqn{u = f_{confluenIUH}(t_r, ...)}
//' 
//' 
//' 
//' where
//' - \mjseqn{u} is series of portions
//' - \mjseqn{t_r} is  `confluen_responseTime_TS`
//' 
//' @references
//' \insertAllCited{}
//' @return IUH (list of num vector) 
//' @details
//' # **_GR4J1** \insertCite{GR4J_Perrin_2003}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_iuh_gr1.svg}}
//' \if{latex}{\figure{mdl_iuh_gr1.pdf}{options: width=100mm}}
//' 
//' \mjsdeqn{u(i) = S(i) - S(i-1)}
//' \mjsdeqn{S(i) = \left( \frac{i}{t_r} \right)^{2.5}, \quad 0 \leq i \leq t_r}
//' where
//'   - \mjseqn{u} is IUH series
//'   - \mjseqn{i} is index
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
//' @details
//' # **_GR4J2** \insertCite{GR4J_Perrin_2003}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_iuh_gr2.svg}}
//' \if{latex}{\figure{mdl_iuh_gr2.pdf}{options: width=100mm}}
//' 
//' \mjsdeqn{u(i) = S(i) - S(i-1)}
//' \mjsdeqn{S(i) = 0.5\left( \frac{i}{t_r} \right)^{2.5}, \quad 0 \leq i \leq t_r}
//' \mjsdeqn{S(i) = 1 - 0.5\left(2 - \frac{i}{t_r} \right)^{2.5}, \quad t_r < i < 2t_r}
//' \mjsdeqn{S(i) = 0, \quad i = 2t_r}
//' where
//'   - \mjseqn{u} is IUH series
//'   - \mjseqn{i} is index
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
//' @details
//' # **_Kelly** \insertCite{iuh_Kelly_1955}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_iuh_kel.svg}}
//' \if{latex}{\figure{mdl_iuh_kel.pdf}{options: width=100mm}}
//' 
//' \mjsdeqn{u(i) = \frac{4}{t_r^2} \left( i + k \left( e^{-i/k} \right) \right), \quad i \leq t_r / 2 }
//' \mjsdeqn{u(i) = - \frac{4}{t_r^2}(i - k - t_r) + \frac{4ke^{-i/k}}{t_r^2} (1 - 2 e^{t_r/(2k)}), \quad t_r / 2 < i \leq t_r }
//' \mjsdeqn{u(i) =  \frac{4ke^{-i/k}}{t_r^2} (1 - 2 e^{t_r/(2k)} +  e^{t_r/k}), \quad i > t_r }
//' where
//'   - \mjseqn{k} is `param_confluen_kel_k`
//' @param param_confluen_kel_k <1, 4> parameter for[confluenIUH_Kelly()]
//' @export
// [[Rcpp::export]]
NumericVector confluenIUH_Kelly(
    double confluen_responseTime_TS,
    double param_confluen_kel_k
)
{
  double confluen_concentratTime_TS = confluen_responseTime_TS * param_confluen_kel_k;
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
//' @details
//' # **_Nash** \insertCite{iuh_Nash_1957}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_iuh_nas.svg}}
//' \if{latex}{\figure{mdl_iuh_nas.pdf}{options: width=100mm}}
//' 
//' \mjsdeqn{u(i) = \frac{1}{t_r\Gamma(n)} \left(\frac{4}{t_r^2}\right)^{n -1}e^{-i/t_r}}
//' where
//'   - \mjseqn{n} is `param_confluen_nas_n`
//' @param param_confluen_nas_n <1, 8> parameter for[confluenIUH_Nash()]
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
//' @details
//' # **_Clark** \insertCite{iuh_Clark_1945}{EDCHM}: 
//'
//' \if{html}{\figure{mdl_iuh_cla.svg}}
//' \if{latex}{\figure{mdl_iuh_cla.pdf}{options: width=100mm}}
//' 
//' \mjsdeqn{u(i) = \frac{1}{t_r} e^{-i/t_r} }
//' where
//'   - \mjseqn{t_r} is `confluen_responseTime_TS`
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
