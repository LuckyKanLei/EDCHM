#include "00utilis.h"
// [[Rcpp::interfaces(r, cpp)]]

//' **potential evapotranspiration**
//' @name evatransPotential
//' @description The concept potential ET estimate mainly the confluence from atmosphere, 
//' when the water in the ET-area is always enough. 
//' Actually we should also consider the charterers of the ET-area, there is free water area or vegetation or bar soil,
//' But we may don't have so many information and effective method, in order to get the truly potential ET.
//' So there is the simplified method: reference ET, it define the ET-area with some fest charterers, e.g. [evatransPotential_FAO56()].
//' In this situation we should give some factor for different ET-area.
//' @references
//' \insertAllCited{}
//' @inheritParams all_vari
//' @details
//' - **_TurcWendling** \insertCite{ET_TurcWendling_1991}{EDCHM}: consider only the radiation and temperature as the main factors. 
//' \mjsdeqn{E_p = \frac{(100 R_s + 3.875 t_h k)\cdot(T + 22)}{150 (T + 123)}}
//' where
//'   - \mjseqn{E_p} is potential ET, `atmos_potentialEvatrans_mm`
//'   - \mjseqn{R_s} is solar radiation, `atmos_solarRadiat_MJ`
//'   - \mjseqn{t_h} is time step in hour, `time_step_h`
//'   - \mjseqn{T} is average air temperature, `atmos_temperature_Cel`
//'   - \mjseqn{k} is `param_evatrans_tur_k`
//' @param param_evatrans_tur_k <0.6, 1> parameter for [evatransPotential_TurcWendling()], higher value when closer to the sea
//' @return potential evapotranspiration (mm/m2)
//' @export
// [[Rcpp::export]]
NumericVector evatransPotential_TurcWendling(
    NumericVector atmos_temperature_Cel, 
    NumericVector atmos_solarRadiat_MJ, 
    NumericVector time_step_h,
    NumericVector param_evatrans_tur_k 
)
{
  return (atmos_solarRadiat_MJ * 100 + 3.875 * time_step_h * param_evatrans_tur_k) * (atmos_temperature_Cel + 22) / 150 / (atmos_temperature_Cel + 123);
}

//' @rdname evatransPotential
//' @details
//' - **_Linacre** \insertCite{ET_Linacre_1977}{EDCHM}: consider only the temperature as the main factors. 
//' \mjsdeqn{E_p = \frac{\frac{100(0.75 - \alpha)(T + 0.006 z)}{100 - \phi} + 15(T - T_d)}{80 - T}}
//' \mjsdeqn{T_d = T - 20 (1-H_R)}
//' where
//'   - \mjseqn{\alpha} is albedo, `land_albedo_1`
//'   - \mjseqn{z} is elevation, `land_elevation_m`
//'   - \mjseqn{T_d} is dewpoint temperature,
//'   - \mjseqn{H_R} is relative humidity, `atmos_relativeHumidity_1`
//' @export
// [[Rcpp::export]]
NumericVector evatransPotential_Linacre(
    NumericVector time_dayOfYear_,
    NumericVector atmos_temperature_Cel,
    NumericVector atmos_relativeHumidity_1,
    NumericVector land_latitude_Degree,
    NumericVector land_elevation_m,
    NumericVector land_albedo_1
)
{
  return ((0.75 - land_albedo_1) * 100 * (atmos_temperature_Cel + 0.006 * land_elevation_m) / (100 - land_latitude_Degree) + 3 * 100 * (1 - atmos_relativeHumidity_1)) / (80 - atmos_temperature_Cel);
}

//' @rdname evatransPotential
//' @details
//' - **_FAO56** \insertCite{ET_FAO56_1998}{EDCHM}: consider not only radiation and temperature but also other variable like wind speed
//' as the main factors. 
//' \mjsdeqn{E_p =\frac{0.408 \Delta\left(R_n - G\right)+\gamma \frac{900}{T+273} {u}_{2}\left({e}_{{s}}-{e}_{{a}}\right)}{\Delta+\gamma\left(1+0.34 {u}_{2}\right)}}
//' where
//'   - \mjseqn{\Delta} is slope vapour pressure curve (kPa °C-1)
//'   - \mjseqn{R_n} is net radiation, `atmos_netRadiat_MJ`
//'   - \mjseqn{G} is soil heat flux density
//'   - \mjseqn{u_2} is wind speed at 2 m height, `atmos_windSpeed2m_m_s`
//'   - \mjseqn{e_s} is saturation vapour pressure, `atmos_saturatVaporPress_hPa`
//'   - \mjseqn{e_a} is actual vapour pressure, `atmos_vaporPress_hPa`
//'   - \mjseqn{\gamma} is psychrometric constant
//' @export
// [[Rcpp::export]]
NumericVector evatransPotential_FAO56(
    NumericVector time_dayOfYear_,
    NumericVector atmos_temperature_Cel, 
    NumericVector atmos_vaporPress_hPa, 
    NumericVector atmos_saturatVaporPress_hPa, 
    NumericVector atmos_netRadiat_MJ,
    NumericVector atmos_windSpeed2m_m_s,
    NumericVector land_latitude_Degree,
    NumericVector land_elevation_m
)
{
  NumericVector Delta_, e_s, e_a, R_n, P_, gamma_, u_2, ET_o;
  
  R_n = atmos_netRadiat_MJ;
  u_2 = atmos_windSpeed2m_m_s;
  e_a = atmos_vaporPress_hPa;
  e_s = atmos_saturatVaporPress_hPa;
  //// Delta_,
  Delta_ = 4098 * (0.6108 * exp(17.27 * atmos_temperature_Cel / (atmos_temperature_Cel + 237.3))) / ((atmos_temperature_Cel + 237.3) * (atmos_temperature_Cel + 237.3)); // Eq.13
  
  //// e_s,
  // e_s = 0.3054 * exp(17.27 * atmos_temperature_max_Cel / (atmos_temperature_max_Cel + 237.3)) + 0.3054 * exp(17.27 * atmos_temperature_min_Cel / (atmos_temperature_min_Cel + 237.3)); // Eq.11, 12
  // 
  // //// e_a,
  // e_a = e_s * atmos_relative_humidity_Perc * 0.01; // Eq.19
  // //// R_n,
  // phi_ = M_PI / 180 * land_latitude_Degree; // Eq.22
  // d_r = 1 + 0.033 * cos(2 * M_PI / 365 * time_dayOfYear_); // Eq.23
  // delta_ = 0.409 * sin(2 * M_PI / 365 * time_dayOfYear_ - 1.39); // Eq.24
  // omega_s = acos(-tan(phi_) * tan(delta_)); // Eq.25
  // //// (24 (60)/M_PI) G_sc = 37.58603, G_sc <- 0.0820 MJ / (m^2 d)
  // R_a = 37.58603 * d_r * (omega_s * sin(phi_) * sin(delta_) + cos(phi_) * cos(delta_) * sin(omega_s)); // Eq.21
  // 
  // R_so = (0.75 + 2e-5 * land_elevation_m) * R_a; // Eq.37
  // R_ns = (1 - alpha_) * solar_radiation_MJ_m2d; // Eq.38
  // R_nl = sigma_ * (pow((atmos_temperature_max_Cel + 273.16), 4) + pow((atmos_temperature_min_Cel + 273.16), 4)) / 2 * (0.34 - 0.14 * sqrt(e_a)) * (1.35 * solar_radiation_MJ_m2d / R_so - 0.35); // Eq.39
  // R_n = R_ns - R_nl; // Eq.40
  
  //// G_,
  // G_ = 0.; // Eq.42
  
  //// gamma_,
  P_ = 101.3 * pow(((293 - 0.0065 * land_elevation_m) / 293), 5.26); // Eq.7
  gamma_ = 0.665e-3 * P_; // Eq.8
  
  //// u_2,
  // u_2 = wind_measure_height_m * 4.87 / (log(67.8 * wind_speed_m_s - 5.42)); // Eq.47
  // u_2 = wind_measure_height_m
  //// TE_o
  ET_o = (0.408 * Delta_ * (R_n - 0.) + gamma_ * 90 * u_2 * (e_s - e_a) / (atmos_temperature_Cel + 273)) / (Delta_ + gamma_ * (1 + 0.34 * u_2));
  return ET_o;
}

//' **actuall evapotranspiration**
//' @name evatransActual
//' @inheritParams all_vari
//' @description
//' \loadmathjax
//' Under the concept of the conceptional HM, the actually ET is always consider as a part of potential ET:
//' \mjsdeqn{E_a = f E_p}
//' where
//' - \mjseqn{E_a} is `land_evatrans_mm` or `soil_evatrans_mm`
//' - \mjseqn{E_p} is `atmos_potentialEvatrans_mm`
//' - \mjseqn{f} is estimated ratio.
//' Then the different `evatransActual` methods will estimate the ratio \mjseqn{f}.
//' @references
//' \insertAllCited{}
//' @return actually ET in (mm/m2/TS)
//' - evaporation in interception (landLy), `land_evatrans_mm`
//' - transpiration in root
//' - evaporation in soil (soilLy), `soil_evatrans_mm`
//' @details
//' - **_SupplyRatio**: the water content (the ratio to the maximal capacity) 
//' is considered as th main factors for the ratio \mjseqn{f}.
//' \mjsdeqn{f = k  \frac{W}{C}}
//' where
//'   - \mjseqn{W} is water volume in (mm/m2/TS), `water_mm`, `land_interceptWater_mm`, `soil_water_mm`
//'   - \mjseqn{C} is water capacity in (mm/m2), `capacity_mm`, `land_interceptCapacity_mm`, `soil_capacity_mm`
//'   - \mjseqn{k} is `param_evatrans_sur_k`
//' @param param_evatrans_sur_k <0.1, 1> parameter for [evatransActual_SupplyRatio()], ratio of potential ET, that is estimated as actually ET  
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_SupplyRatio(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_evatrans_sur_k
)
{
  NumericVector AET, k_;
  
  k_ = water_mm / capacity_mm * param_evatrans_sur_k;
  AET = atmos_potentialEvatrans_mm * k_;
  return ifelse(AET > water_mm, water_mm, AET);
}

//' @rdname evatransActual
//' @details
//' - **_SupplyPow**: the water content (the ratio to the maximal capacity) 
//' is considered as th main factors for the ratio \mjseqn{f}.
//' \mjsdeqn{f = k  \left(\frac{W}{C}\right)^\gamma}
//' where
//'   - \mjseqn{k} is `param_evatrans_sup_k`
//'   - \mjseqn{\gamma} is `param_evatrans_sup_gamma`
//' @param param_evatrans_sup_k <0.1, 1> parameter for [evatransActual_SupplyPow()], ratio of this method
//' @param param_evatrans_sup_gamma <1, 5> parameter for [evatransActual_SupplyPow()], exponent of this method
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_SupplyPow(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_evatrans_sup_k,
    NumericVector param_evatrans_sup_gamma
)
{
  NumericVector AET, k_;
  
  k_ = param_evatrans_sup_k * vecpow((water_mm / capacity_mm), param_evatrans_sup_gamma);
  AET = atmos_potentialEvatrans_mm * k_;
  return ifelse(AET > water_mm, water_mm, AET);
}



//' @rdname evatransActual
//' @details
//' - **_VIC** \insertCite{VIC_Wood_1992}{EDCHM}: is similar with [evatransActual_SupplyPow()], estimate the water content in the storage.
//' \mjsdeqn{f = 1-\left(1-\frac{W}{C}\right)^{\gamma}}
//' where
//'   - \mjseqn{\gamma} is `param_evatrans_vic_gamma`
//' @param param_evatrans_vic_gamma <0.2, 5> parameter for [evatransActual_VIC()]
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_VIC(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_evatrans_vic_gamma
)
{
  NumericVector AET, k_;
  
  k_ = 1 - vecpow((1- water_mm / capacity_mm), param_evatrans_vic_gamma);
  AET = atmos_potentialEvatrans_mm * k_;
  return ifelse(AET > water_mm, water_mm, AET);
}




//' @rdname evatransActual
//' @details
//' - **_GR4J** \insertCite{GR4J_Perrin_2003}{EDCHM}: is a little different than other method, it estimate not the ratio \mjseqn{f},
//' rather dieectly a equation with potential ET and water content.
//' And it need **no parameter**.
//' \mjsdeqn{E_a = \frac{W\left(2-\frac{W}{C}\right)\tanh \left(\frac{E_p}{C}\right)}{1 + \left(1-\frac{W}{C}\right)\tanh \left(\frac{E_p}{C}\right)}}
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_GR4J(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm
)
{
  NumericVector AET;
  
  
  AET = water_mm * (2 - water_mm / capacity_mm) * tanh(atmos_potentialEvatrans_mm / capacity_mm) / (1 + (1 - water_mm / capacity_mm) * tanh(atmos_potentialEvatrans_mm / capacity_mm));
  return ifelse(AET > water_mm, water_mm, AET);
  
}

//' @rdname evatransActual
//' @details
//' - **_UBC** \insertCite{VIC_Wood_1992}{EDCHM}: estimate the water content in the storage. 
//' (This is a little different than original, the parameter `P0AGEN` is replaced by \mjseqn{\frac{C}{\gamma}}.)
//' \mjsdeqn{f = 10^{\gamma \frac{W-C}{C}}}
//' where
//'   - \mjseqn{\gamma} is `param_evatrans_ubc_gamma`
//' @param param_evatrans_ubc_gamma <0.5, 2>parameter for [evatransActual_UBC()]
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_UBC(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_evatrans_ubc_gamma
)
{
  NumericVector diff_mm, AET, k_;
  diff_mm = capacity_mm - water_mm;
  
  
  
  k_ = vecpow10(- diff_mm / (param_evatrans_ubc_gamma * capacity_mm));
  AET = atmos_potentialEvatrans_mm * k_;
  return ifelse(AET > water_mm, water_mm, AET);
}

//' @rdname evatransActual
//' 
//' @details
//' - **Land_Liang** \insertCite{VIC2_Liang_1994}{EDCHM}: is also a similar method like [evatransActual_SupplyPow()], 
//' but it will estimate the supply ability agian, whwn the water is still not enough.
//' \mjsdeqn{E_l^* = \left(\frac{W}{C}\right)^\gamma E_p}
//' \mjsdeqn{E_l = \min \left(1, \frac{W}{E_l^*}\right) E_l^*}
//' where
//'   - \mjseqn{E_l^*} is the first estimated actuall ET
//'   - \mjseqn{E_l} is actuall ET from land, `land_evatrans_mm`
//'   - \mjseqn{\gamma} is `param_evatrans_lia_gamma`
//' @param param_evatrans_lia_gamma <0.4, 1> parameter for [evatransLand_Liang()]
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_LiangLand(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_evatrans_lia_gamma
)
{
  NumericVector AET, k_, f_;
  
  k_ = vecpow((water_mm / capacity_mm), param_evatrans_lia_gamma);
  AET = atmos_potentialEvatrans_mm * k_;
  f_ = water_mm / AET;
  f_ = ifelse(f_ > 1, 1, f_);
  AET = f_ * AET;
  return ifelse(AET > water_mm, water_mm, AET);
}


//' @rdname evatransActual
//' @details
//' - **Soil_Liang** \insertCite{VIC2_Liang_1994}{EDCHM}: estimate the water content in the storage. 
//' (This is a little different than original, the parameter `P0AGEN` is replaced by \mjseqn{\frac{C}{\gamma}}.)
//' \mjsdeqn{f = \int_{0}^{A_{s}} d A + \int_{A_{s}}^{1} \frac{i_{0}}{i_{m} [1-(1-A)^{1 / B} ]} d A }
//' where
//'   - \mjseqn{B} is `param_evatrans_lia_B`
//' @param param_evatrans_lia_B <0.01, 3> parameter for [evatransSoil_Liang()]
//' @export
// [[Rcpp::export]]
NumericVector evatransActual_LiangSoil(
    NumericVector atmos_potentialEvatrans_mm,
    NumericVector water_mm,
    NumericVector capacity_mm,
    NumericVector param_evatrans_lia_B
)
{
  NumericVector AET, k_, f_, i_0, i_m, B_1, B_p_1, A_s, A_s_1;
  
  i_m = capacity_mm * (param_evatrans_lia_B + 1);
  
  B_p_1 = (param_evatrans_lia_B + 1);
  B_1 = 1 / B_p_1;
  
  i_0 = i_m * (1 - vecpow(1 - water_mm / capacity_mm, B_1));
  
  // B_p_1 = (param_evatrans_lia_B + 1);
  // B_1 = 1 / param_evatrans_lia_B;
  // 
  // i_0 = i_m * (1 - vecpow(1 - water_mm * B_p_1 / i_m, B_1));
  A_s = 1 - vecpow((1 - i_0 / i_m), param_evatrans_lia_B);
  A_s_1 = (1 - A_s);
  
  k_ = A_s + i_0 / i_m * A_s_1 * (1 + 
    param_evatrans_lia_B / (1+param_evatrans_lia_B) * vecpow(A_s_1, 1 / param_evatrans_lia_B) +
    param_evatrans_lia_B / (2+param_evatrans_lia_B) * vecpow(A_s_1, 2 / param_evatrans_lia_B) +
    param_evatrans_lia_B / (3+param_evatrans_lia_B) * vecpow(A_s_1, 3 / param_evatrans_lia_B));
  
  AET = atmos_potentialEvatrans_mm * k_;
  
  return ifelse(AET > water_mm, water_mm, AET);
}

