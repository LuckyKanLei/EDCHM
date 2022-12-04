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
#' - Layer:
#'    - atmosLy
#'    - landLy (snowLy)
#'    - soilLy
#'    - groundLy
#' - Boundary:
#'   - topBd
#'   - atmosBd
#'   - landBd
#'   - groundBd
#'   - btmBd
#'   
#' ![](edchm_layer.pdf){options: width=165mm}
#' 
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
#' For the formulas symbols there are also some conventions:
#' 
#' - \mjseqn{W}: water volume in one Layer (mm/m2)
#' 
#'    - \mjseqn{W_{land}}: pounded water in `landLy` 
#'    - \mjseqn{W_{itcp}}: intercepted water in `landLy`
#'    - \mjseqn{W_{snow}}: equal water volume in `snowLy`
#'    - \mjseqn{W_{soil}}: storage in `soilLy`
#'    - \mjseqn{W_{grnd}}: storage in `groundLy`
#'    
#' - \mjseqn{C}: maximal capacity of water volume in one Layer (mm/m2)
#' 
#'    - \mjseqn{C_{itcp}}: intercepted water capacity in `landLy`
#'    - \mjseqn{C_{soil}}: storage capacity in `soilLy`
#'    - \mjseqn{C_{grnd}}: storage capacity in `groundLy`
#' 
#' - \mjseqn{P}: flux of precipitation (mm/m2/TS)
#' 
#'    - \mjseqn{P_r}: rain fall
#'    - \mjseqn{P_s}: snow fall
#'    
#' - \mjseqn{E}: flux of evapotranspiration (mm/m2/TS)
#' 
#'    - \mjseqn{E_p}: potential evapotranspiration
#'    - \mjseqn{E_a}: actual evapotranspiration from one specific Layer (e.g. `soilLy` or `landLy`)
#'    
#' - \mjseqn{I}: flux of infiltration (mm/m2/TS)
#' - \mjseqn{M}: flux or snow melt (mm/m2/TS) 
#' - \mjseqn{R}: flux of runoff (mm/m2/TS), that from precipitation travels over the soil surface to the nearest stream channel
#' - \mjseqn{S}: flux of subsurface flow (mm/m2/TS), that from `soilLy` travels over the soil subsurface to the nearest stream channel
#' - \mjseqn{B}: flux of baseflow (mm/m2/TS), that from `groundLy` travels groundwater to the nearest stream channel
#' 
#' - \mjseqn{A}: flux of capillary rise (mm/m2/TS) (A from German Aufstieg)
#' - \mjseqn{V}: flux of percolation (mm/m2/TS) (V from German Versinkung)
#' - \mjseqn{L}: flux of lateral flow (mm/m2/TS), the exchange with the other catchment, it can be positive (into) or negative (go out)
#' 
#' - \mjseqn{Q}: streamflow in the river section (mm/m2/TS) (PS: **the results from EDCHM has't transform the unit from flux (mm/m2/TS) to the streamflow (m3/TS) or (m3/s)**)
#' 
#'    - \mjseqn{Q_r}: proportion of \mjseqn{Q} from runoff (\mjseqn{R})
#'    - \mjseqn{Q_s}: proportion of \mjseqn{Q} from subsurface flow (\mjseqn{S})
#'    - \mjseqn{Q_b}: proportion of \mjseqn{Q} from baseflow (\mjseqn{B})
#' 
#' e.g. **group_variableName_unit**
#' # Parameter define
#' The Parameter will defined in every function topic, but there will define the naming regulation:
#' Parameters make up with prefix `param`, Process name (sometimes same as Layer name), 
#' abbreviation in threee small letters of the Method, and the Parameter name in the original.
#' 
#' e.g. **param_processName_mtd_k**
#' 
#' # Process
#' There are now 10 Process available:
#' 
#' ![](hydro_modula_structure.pdf){options: width=165mm}
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