#' variable define
#' @name all_vari
#' @description In this topic all of the variable in the EDCHM will be defined.
#' 
#' 
#' @param time_step_h <1, 24> (h) time step in hour
#' @param time_dayOfYear_ <1, 366> (d) the number of the day in the year between 1 (1 January) and 365 or 366 (31 December)
#' @param atmos_precipitation_mm (mm/m2/TS) precipitaion volum
#' @param atmos_rain_mm (mm/m2/TS) precipitation in rain form
#' @param atmos_snow_mm (mm/m2/TS) precipitation in snow form
#' @param atmos_temperature_Cel (Cel) the average air temperature in the time phase
#' @param atmos_temperatureMax_Cel (Cel) the maximal air temperature in the time phase
#' @param atmos_temperatureMin_Cel (Cel) the minimal air temperature in the time phase
#' @param atmos_solarRadiat_MJ (MJ/m2/TS) the solar radiation that actually reaches the earths surface
#' @param atmos_netRadiat_MJ	(MJ/m2/TS) the balance between the energy absorbed, reflected and emitted by the earths surface or the difference between the incoming net shortwave (Rns) and the net outgoing longwave (Rnl) radiation
#' @param atmos_relativeHumidity_1 (0, 1) relative humidity
#' @param atmos_vaporPress_hPa (hPa) actual vapour pressure, can be calculated by [atmos_VaporPress()]
#' @param atmos_saturatVaporPress_hPa (hPa) saturation vapour pressure at `atmos_temperature_Cel`, can be calculated by [atmos_SaturatVaporPress()]
#' @param atmos_windSpeed_m_s (m/s) measured wind speed at z m above ground surface
#' @param atmos_windMeasureHeight_m (m) height of measurement above ground surface
#' @param atmos_windSpeed2m_m_s (m/s) wind speed at 2 m above ground surface
#' @param atmos_potentialEvatrans_mm (mm/m2/TS) **potential / reference** evapotranspiration
#' @param land_albedo_1 <0, 1> albedo of the region
#' @param land_latitude_Degree (degree) average latitude
#' @param land_elevation_m (m) average elevation
#' @param land_water_mm (mm/m2) **pounded water** volume in `landLy` and there is no limit, different than `land_interceptWater_mm`
#' @param land_interceptWater_mm (mm/m2) initial water volume that can be **intercepted**
#' @param land_interceptCapacity_mm (mm/m2) average intercept Capacity (maximal storage capacity)
#' @param snow_ice_mm (mm/m2) water equivalent of **ice** in snowpack
#' @param soil_water_mm (mm/m2) water volume in `soilLy`
#' @param soil_capacity_mm (mm/m2) average soil Capacity (maximal storage capacity)
#' @param ground_water_mm (mm/m2/TS) water volume in `groundLy`
#' @param ground_capacity_mm (mm/m2) water storage capacity in `groundLy`
#' @param ground_lateral_mm (mm/m2/TS) lateral flow, exchange with outside region. It can be **NEGATIV**
#' @param ground_lateralPotential_mm (mm/m2/TS) potential lateral flow, that will as parameter evaluated
#' @param confluen_inputWater_mm,land_runoff_mm,ground_baseflow_mm (mm/m2) input water volum in every routeline
#' @param confluen_iuh_1,confluen_iuhLand_1,confluen_iuhGround_1 (vector of num, sume() = 1) the ratio in every timestep, can be calculated by [confluenIUH_GR4J1()], [confluenIUH_GR4J2()]
#' @param confluen_responseTime_TS,confluen_concentratTime_TS (TS) response or concentration time in every routeline
#' @param water_mm (mm/m2/TS) water volume in `soilLy` or interceptof `landLy`, same as `land_interceptWater_mm` or `soil_water_mm`
#' @param capacity_mm (mm/m2) water storage capacity in `soilLy` or interceptof `landLy`, same as `land_interceptCapacity_mm` or `soil_capacity_mm`
NULL




