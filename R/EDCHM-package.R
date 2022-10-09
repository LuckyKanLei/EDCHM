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
#' # Process
#' 
#' There are now 10 Process available:
#' 
#' - atmos: caculate the basic variable (data) in the `atmosLy`, sometimes the variable komm directly from the original data
#'    - atmosSnow: rain and snoe to split
#' - evatrans: evapotranspiration
#'    - evatransPotential: potential ET calculate
#'    - evatransActual: actual ET calculate
#'    - evatransLand: actual ET in `landLy` (evaporation from interception)
#'    - evatransSoil: actual ET in `soilLy` (transpiration from root zone + evaporation from top soil)
#' - snow: snow
#'    - snowMelt:
#' - intercep: interception
#' - infilt: infiltration `landLy` to `soilLy`
#' - capirise: capital rise, `groundLy` to `soilLy`
#' - percola: percolation, `soilLy` to `groundLy`
#' - baseflow: interception, `groundLy`
#' - lateral: lateral flow, exchange with outside region, `groundLy`
#' - confluen: confluence to river
#'    - confluenIUH: IUH to calculate
#' 
#' 
#' @useDynLib EDCHM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
#' @import mathjaxr
NULL