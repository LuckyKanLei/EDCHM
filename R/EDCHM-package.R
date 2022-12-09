#' basic concept of EDCHM
#' @docType package
#' @name EDCHM-package
#' @description 
#' EDCHM is a conceptional hydrologist modelling framework, 
#' with that we can easy assembly a new hydrologist model with difficult process-method (modular).
#' 
#' # Layer and Boundary
#' 
#' Under EDCHM the **Layer** is defined layer with vertical space, in that have some common characters.
#' And the Layer always have top and bottom two **Boundary**, 
#' one Boundary in the same time the top Boundary of one Layer but also the bottom Layer of other Layer, that on the above.
#' 
#' So for EDCHM there are five basic Layers (now one additional snowLy) and six Boundaries:
#' 
#' - Layer:
#'    - atmosLy
#'    - landLy (snowLy)
#'    - soilLy
#'    - groundLy
#'    
#' - Boundary:
#'   - topBd
#'   - atmosBd
#'   - landBd
#'   - groundBd
#'   - btmBd
#'   
#' \if{html}{\figure{edchm_layer.svg}}
#' \if{latex}{\figure{edchm_layer.pdf}{options: width=165mm}}
#'   
#' # Variable define
#' 
#' All the names of the variables make up of three parties: **group-name**, **physical-name** and **variable-units**
#' and split with underline `_`.
#' - group-name: the Layer or Boundary name, which is defined in EDCHM. 
#' Sometimes it can also be the Process name or `time`, one important dimension.
#' But the group-name is no longer than 8 and all small letters.
#' - physical-name: the physical variable name, this part will with _camelCase_ to combine more words.
#' - variable-units: the physical units of the variable, in order to simplify the program, 
#' all the time dependent variable will use the `TS` as the time units, it will be defined by the model.
#' But the time units part will not appear in the name. 
#' And also all the variable meant they are homogeneous in the area, so all the variable will not give the square meter (m2) in the name.
#' 
#' e.g. **group_variableName_unit**
#' # Parameter define
#' The Parameter will defined in every function topic, but there will define the naming regulation:
#' Parameters make up with prefix `param`, Process name (sometimes same as Layer name), 
#' abbreviation in threee small letters of the Method, and the Parameter name in the original.
#' 
#' e.g. **param_processName_mtd_k**
#' 
#' \loadmathjax
#' 
#' The tables show the collection model variables and the formula symbols:
#' 
#' - some state variables:
#' 
#' | **Variable**                  | **Symbol**        | **Unit** | **Description**                            |
#' |-------------------------------|-------------------|----------|--------------------------------------------|
#' | `water_mm`                    | \mjseqn{W}        | mm/m2    | water volume in one `Layer`                |
#' | `land_water_mm`               | \mjseqn{W_{land}} | mm/m2    | .. in `landLy`                             |
#' | `land_interceptWater_mm`      | \mjseqn{W_{itcp}} | mm/m2    | .. in `landLy` (intercepted)               |
#' | `snow_ice_mm`                 | \mjseqn{W_{snow}} | mm/m2    | .. in `snowLy` (equal water)               |
#' | `soil_water_mm`               | \mjseqn{W_{soil}} | mm/m2    | .. in `soilLy`                             |
#' | `ground_water_mm`             | \mjseqn{W_{grnd}} | mm/m2    | .. in `groundLy`                           |
#' | `capacity_mm`                 | \mjseqn{C}        | mm/m2    | maximal capacity of storage in one `Layer` |
#' | `land_interceptCapacity_mm`   | \mjseqn{C_{itcp}} | mm/m2    | .. in `landLy` (intercepted)               |
#' | `soil_capacity_mm`            | \mjseqn{C_{soil}} | mm/m2    | .. in `soilLy`                             |
#' | `ground_capacity_mm`          | \mjseqn{C_{grnd}} | mm/m2    | . in `groundLy`                            |
#' 
#' - some flux variables:
#' 
#' | **Variable**                  | **Symbol**        | **Unit** | **Description**                            |
#' |-------------------------------|-------------------|----------|--------------------------------------------|
#' | `flux_mm`                     | \mjseqn{F}        | mm/m2/TS | flux or flow in unit area                  |
#' | `atmos_precipitation_mm`      | \mjseqn{P}        | mm/m2/TS | .. of precipitation                        |
#' | `atmos_rain_mm`               | \mjseqn{P_r}      | mm/m2/TS | .. of rain fall                            |
#' | `atmos_snow_mm`               | \mjseqn{P_s}      | mm/m2/TS | .. of snow fall                            |
#' | `atmos_evatrans_mm`           | \mjseqn{E_a}      | mm/m2/TS | .. of evapotranspiration                   |
#' | `land_intercept_mm`           | \mjseqn{F_{iflt}} | mm/m2/TS | .. of interception                         |
#' | `land_infilt_mm`              | \mjseqn{F_{iflt}} | mm/m2/TS | .. of infiltration                         |
#' | `land_runof_mm`               | \mjseqn{F_{roff}} | mm/m2/TS | .. of runoff                               |
#' | `snow_melt_mm`                | \mjseqn{F_{melt}} | mm/m2/TS | .. of snow melt                            |
#' | `soil_percola_mm`             | \mjseqn{F_{pecl}} | mm/m2/TS | .. of percolation                          |
#' | `soil_interflow_mm`           | \mjseqn{F_{intf}} | mm/m2/TS | .. of interflow                            |
#' | `ground_baseflow_mm`          | \mjseqn{F_{base}} | mm/m2/TS | .. of baseflow                             |
#' | `ground_capillarise_mm`       | \mjseqn{F_{capi}} | mm/m2/TS | .. of capillary rise                       |
#' | `ground_lateral_mm`           | \mjseqn{F_{ltrl}} | mm/m2/TS | .. of lateral flow                         |
#' | `potentialFlux_mm`            | \mjseqn{M}        | mm/m2/TS | potential (maximal) flux or flow           |
#' | `atmos_evatrans_mm`           | \mjseqn{E_p}      | mm/m2/TS | .. of evapotranspiration                   |
#' | `land_potentialInfilt_mm`     | \mjseqn{M_{iflt}} | mm/m2/TS | .. of infiltration                         |
#' | `soil_potentialPercola_mm`    | \mjseqn{M_{pecl}} | mm/m2/TS | .. of percolation                          |
#' | `soil_potentialInterflow_mm`  | \mjseqn{M_{intf}} | mm/m2/TS | .. of subsurface flow                      |
#' | `ground_potentialBaseflow_mm` | \mjseqn{M_{base}} | mm/m2/TS | .. of baseflow                             |
#' | `soil_potentialCapirise_mm`   | \mjseqn{M_{capi}} | mm/m2/TS | .. of capillary rise                       |
#' | `ground_potentialLateral_mm`  | \mjseqn{M_{ltrl}} | mm/m2/TS | .. of lateral flow                         |
#' 
#' - the stream flow will not in m3/TS or m3/s but also in flux dimension:
#' 
#' | **Variable**                  | **Symbol**        | **Unit** | **Description**                            |
#' |-------------------------------|-------------------|----------|--------------------------------------------|
#' | `streamflow_mm`               | \mjseqn{Q}        | mm/m2/TS | streamflow in flux dimension               |
#' | `flow_runoff_mm`              | \mjseqn{Q_{roff}} | mm/m2/TS | .. from runoff                             |
#' | `flow_interflow_mm`           | \mjseqn{Q_{itfl}} | mm/m2/TS | .. from interflow                          |
#' | `flow_baseflow_mm`            | \mjseqn{Q_{base}} | mm/m2/TS | .. from baseflow                           |
#' 
#' 
#'    
#' 
#' Additional there are also some symbols from the program-view:
#' 
#' - \mjseqn{D}: the collection of all data of one group or layer
#' 
#'    - \mjseqn{D_{atms}}: data in `atmosLy`
#'    - \mjseqn{D_{land}}: data in `landLy` 
#'    - \mjseqn{D_{snow}}: data in `snowLy`
#'    - \mjseqn{D_{soil}}: data in `soilLy`
#'    - \mjseqn{D_{grnd}}: data in `groundLy`
#'    - \mjseqn{D_{lssg}}: data in any (but the one) `landLy`, `snowLy`, `soilLy` or `groundLy`
#' 
#' - \mjseqn{f}: function or modular e.g. \mjseqn{f_{atmosSnow}} or \mjseqn{f_{inflt}}
#' 
#' 
#' 
#' # Process
#' There are now 10 Process available:
#' 
#' \if{html}{\figure{hydro_modula_structure.svg}}
#' \if{latex}{\figure{hydro_modula_structure.pdf}{options: width=165mm}}
#' 
#' - atmos: caculate the basic variable (data) in the `atmosLy`, sometimes the variable komm directly from the original data
#'    - atmosSnow: rain and snoe to split
#' - evatrans: evapotranspiration
#'    - evatransPotential: potential ET calculate
#'    - evatransActual: actual ET calculate
#' - snow: snow
#'    - snowMelt:
#' - intercep: interception
#' - infilt: infiltration `landLy` to `soilLy`
#' - capirise: capital rise, `groundLy` to `soilLy`
#' - percola: percolation, `soilLy` to `groundLy`
#' - baseflow: interception, `groundLy`
#' - lateral: lateral flow, exchange with outside region, `groundLy`
#' - confluen: confluence from catchment (or subcatchment) to river
#'    - confluenIUH: IUH to calculate
#' 
#' 
#' 
#' @references
#' \insertAllCited{}
#' @useDynLib EDCHM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
#' @import mathjaxr
NULL