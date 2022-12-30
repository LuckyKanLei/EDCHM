#include "EDCHM_GR4J.h"
// [[Rcpp::interfaces(r, cpp)]]

double sum_product(NumericVector lhs,
                   NumericVector rhs)
{
  int n = lhs.size();
  double s = 0.0;
  for (int i= 0; i < n; i++) {
    s += lhs(i) * rhs(i);
  }
  return s;
}
void resetVector(Rcpp::NumericVector& x) {
  // Fill the vector with zeros
  std::fill(x.begin(), x.end(), 0.0);
}

//' @name modells
//' @param S_,R_ storage water S and R
//' @param X_1,X_2,X_3,X_4 parameters in GR4J
//' @details
//' # **EDCHM_GR4J** \insertCite{GR4J_Perrin_2003}{EDCHM}: 
//' 
//' Total same like original GR4J
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix EDCHM_GR4J(
    int n_time,
    int n_spat,
    NumericMatrix atmos_potentialEvatrans_mm,
    NumericMatrix atmos_precipitation_mm,
    NumericVector S_,
    NumericVector R_,
    NumericVector X_1, // x1
    NumericVector X_2, // x2
    NumericVector X_3, // x3
    NumericVector X_4 // x4
)
{

  int n_UH_land = ceil(max(X_4) * 2), n_UH_ground = ceil(max(X_4));
  NumericVector P_n, Q_r, E_s, P_s, Perc_, P_r, Pr_1, Pr_9,
  Q_1(n_spat), Q_9(n_spat), F_, Q_d, F_9, P_, E_, E_n, A_E, temp, v_temp(n_UH_land);
  NumericMatrix  Q_(n_time, n_spat),
  out_S(n_time, n_spat), out_Q9(n_time, n_spat), out_Q1(n_time, n_spat),
  out_Perc(n_time, n_spat), out_Pr(n_time, n_spat), out_AE(n_time, n_spat),
  out_R(n_time, n_spat), out_Qr(n_time, n_spat), out_Qd(n_time, n_spat),
  mat_Pr_1(n_UH_land, n_spat), mat_Pr_9(n_UH_ground, n_spat),
  UH_2(n_UH_land, n_spat), UH_1(n_UH_ground, n_spat);
  
  
  for (int j= 0; j < n_spat; j++) {
    resetVector(v_temp);
    temp = confluenIUH_GR4J2(X_4(j));
    std::copy(temp.begin(), temp.end(), v_temp.begin());
    UH_2(_, j) = v_temp;
    
    resetVector(v_temp);
    temp = confluenIUH_GR4J1(X_4(j));
    std::copy(temp.begin(), temp.end(), v_temp.begin());
    UH_1(_, j) = v_temp;
    
  }
  
  
  for (int i= 0; i < n_time; i++) {
    
    P_ = atmos_precipitation_mm(i, _);
    E_ = atmos_potentialEvatrans_mm(i, _);
    
    
    
    
    P_n = ifelse(P_ > E_, P_ - E_, 0.0);
    E_n = ifelse(P_ > E_, 0.0, E_ - P_);
    P_n = ifelse(P_n > 13 * X_1, 13 * X_1, P_n);
    E_n = ifelse(E_n > 13 * X_1, 13 * X_1, E_n);
    P_s = infilt_GR4J(P_n, S_, X_1);
    E_s = evatransActual_GR4J(E_n, X_1, S_);
    
    P_s = ifelse(P_ > E_, P_s, 0.0);
    E_s = ifelse(P_ > E_, 0.0, E_s);
    A_E = ifelse(P_ > E_, E_, E_s + P_);
    S_ += (P_s - E_s);
    
    Perc_ = percola_GR4J(S_, X_1);
    S_ +=  - Perc_;
    
    P_r = (P_n - P_s + Perc_);
    P_r = ifelse(P_r < 0, 0, P_r);
    Pr_1 = 0.1 * P_r; //
    Pr_9 = 0.9 * P_r;
    
    
    for (int ii = n_UH_land - 1; ii > 0; ii--) {
      mat_Pr_1(ii,_) = mat_Pr_1(ii-1, _);
      
    }
    mat_Pr_1(0, _) = Pr_1;
    
    for (int ii = n_UH_ground - 1; ii > 0; ii--) {
      mat_Pr_9(ii,_) = mat_Pr_9(ii-1, _);
      
    }
    mat_Pr_9(0, _) = Pr_9;
    
    for (int j= 0; j < n_spat; j++) {
      Q_9(j) = sum_product(mat_Pr_9(_, j), UH_1(_, j)); //as<double>(m1m1_mult(as<arma::mat>(mat_Pr_9(_, j)), as<arma::mat>(UH_1(_, j))));
      Q_1(j) = sum_product(mat_Pr_1(_, j), UH_2(_, j)); //as<double>(m1m1_mult(as<arma::mat>(mat_Pr_1(_, j)), as<arma::mat>(UH_2(_, j))));
      
    }
    
    F_ = lateral_GR4J(R_, X_3, X_2);
    Q_d = ifelse((Q_1 + F_) > 0.0, Q_1 + F_, 0) ;
    
    R_ += (Q_9 + F_);
    R_ = ifelse(R_ > 0.0, R_, 0.0) ;
    Q_r = baseflow_GR4J(R_, X_3);
    R_ +=  - Q_r;
    
    Q_(i,_) = Q_r + Q_d;
    
    
    // out_S(i,_) = S_;
    // out_Q9(i,_) = Q_9;
    // out_Q1(i,_) = Q_1;
    // out_Perc(i,_) = Perc_;
    // out_Pr(i,_) = Pr_1 + Pr_9;
    // out_AE(i,_) = A_E;
    // out_R(i,_) = R_;
    // out_Qr(i,_) = Q_r;
    // out_Qd(i,_) = Q_d;
    
    
    
  }
  // return List::create(
  //   _["S"] = out_S,
  //   _["Q9"] = mat_Pr_9,
  //   _["Q1"] = mat_Pr_1,
  //   _["Perc"] = out_Perc,
  //   _["Pr"] = out_Pr,
  //   _["AE"] = out_AE,
  //   _["R"] = out_R,
  //   _["Qr"] = out_Qr,
  //   _["Qd"] = out_Qd,
  //   _["Q"] = Q_,
  //   _["UH1"] = UH_1,
  //   _["UH2"] = UH_2
  // );
  return Q_;
}

