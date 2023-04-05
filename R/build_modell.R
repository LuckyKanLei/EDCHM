#' build model with **EDCHM-standard-structure**
#' @description Build a hydrological model (based on c++) with EDCHM modules.
#' The **EDCHM-standard-structure** include 15 process:
#'
#' \if{html}{\figure{structure_standard.svg}}
#' \if{latex}{\figure{structure_standard.pdf}{options: width=140mm}}
#' 
#' Some of them can be set in `"NULL"`, means without this process.
#' - `stmosSnow` with input data
#' - `stmosSnow` and `snowMelt`
#' - `evatransPotential` with input data
#' - `intercep` and `evatransLand` (`evatransActual`)
#' - `intflow` and `confluenSoil` (`confluenIUH`)
#' - `capirise`
#' - `lateral`
#' @importFrom stringr str_remove str_remove_all str_replace str_replace_all str_split str_split_fixed str_which
#' @importFrom utils read.csv
#' @importFrom purrr reduce map
#' @param process_method named char vector, some like example.
#' It must contian all the 15 processes. But you can set process in `"NULL"`, when it not necessary.
#' @param path_model char of path, path to space the c++ source files
#' @param name_model char, name of the model
#' @examples
#' my_process_method <- c( 
#' atmosSnow = "atmosSnow_ThresholdT",
#' evatransPotential = "NULL",
#' evatransLand = "evatransActual_SupplyRatio",
#' evatransSoil = "evatransActual_SupplyRatio",
#' intercep = "intercep_Full",
#' snowMelt = "snowMelt_Kustas",
#' infilt = "infilt_AcceptRatio",
#' percola = "percola_SupplyRatio",
#' inteflow = "inteflow_GR4Jfix",
#' capirise = "capirise_HBV",
#' baseflow = "baseflow_SupplyRatio",
#' lateral = "lateral_GR4Jfix",
#' confluenLand = "confluenIUH_Kelly",
#' confluenSoil = "confluenIUH_Kelly",
#' confluenGround = "confluenIUH_Nash"
#' )
#' @return the range of parameters 
#' @export
build_modell <- function(process_method, name_model, path_model = NULL) {
  
  vari_initial <- c("land_interceptWater_mm", "soil_water_mm", "ground_water_mm", "snow_ice_mm")
  vari_boundary <- c("atmos_precipitation_mm",
                     "atmos_temperature_Cel",
                     "atmos_relativeHumidity_1", "atmos_saturatVaporPress_hPa", "atmos_vaporPress_hPa",
                     "atmos_solarRadiat_MJ", "atmos_netRadiat_MJ",
                     "atmos_windSpeed2m_m_s")
  
  process_vari <- c(atmosSnow = "atmos_snow_mm",
                    evatransPotential = "atmos_potentialEvatrans_mm",
                    evatransLand = "land_evatrans_mm",
                    evatransSoil = "soil_evatrans_mm",
                    intercep = "land_intercept_mm",
                    snowMelt = "snow_melt_mm",
                    infilt = "soil_infilt_mm",
                    percola = "soil_percolation_mm",
                    inteflow = "soil_interflow_mm",
                    capirise = "ground_capilarise_mm",
                    baseflow = "ground_baseflow_mm",
                    lateral = "ground_lateral_mm",
                    confluenLand = "confluenLand_iuh_1",
                    confluenSoil = "confluenSoil_iuh_1",
                    confluenGround = "confluenGround_iuh_1")
  process_after <- c(atmosSnow = "atmos_precipitation_mm = atmos_precipitation_mm - atmos_snow_mm;",
                     evatransPotential = "",
                     evatransLand = "land_interceptWater_mm += - land_evatrans_mm;",
                     evatransSoil = "soil_water_mm += - soil_evatrans_mm;",
                     intercep = "land_interceptWater_mm += land_intercept_mm;\nland_water_mm = atmos_precipitation_mm - land_intercept_mm;",
                     snowMelt = "land_water_mm += snow_melt_mm;\nsnow_ice_mm += -snow_melt_mm;\nsnow_ice_mm += atmos_snow_mm;",
                     infilt = "soil_water_mm += soil_infilt_mm;\nland_runoff_mm = land_water_mm - soil_infilt_mm;",
                     percola = "ground_water_mm += soil_percolation_mm;\nsoil_water_mm += - soil_percolation_mm;",
                     inteflow = "soil_water_mm += - soil_interflow_mm;",
                     capirise = "ground_water_mm += - ground_capilarise_mm;\nsoil_water_mm += ground_capilarise_mm;",
                     baseflow = "ground_water_mm += - ground_baseflow_mm;",
                     lateral = "ground_water_mm += ground_lateral_mm;",
                     confluenLand = "",
                     confluenSoil = "",
                     confluenGround = "")
  
  # MODULE ALL ----------
  lines_h <- readLines(file.path(path.package("EDCHM"), "include/EDCHM_RcppExports.h"))
  lines_process <- lines_h[str_which(lines_h, "    inline")] |>
    str_remove_all("NumericVector ") |>
    str_remove_all("double ") |>
    str_remove("    inline ") |>
    str_remove(" \\{")
  module_name <- str_split_fixed(lines_process, "\\(", 2)[,1]
  
  
  # CHECK NULL ----------
  if (sum(process_method[c("evatransLand", "intercep")] == "NULL") == 1) stop("The process `evatransLand` and `intercep` kann in the same tiem with `NULL`, but can't only one of them is `NULL`.")
  if (sum(process_method[c("inteflow", "confluenSoil")] == "NULL") == 1) stop("The process `inteflow` and `confluenSoil` kann in the same tiem with `NULL`, but can't only one of them is `NULL`.")
  if (process_method["atmosSnow"] != "NULL" & process_method["snowMelt"] == "NULL" ) stop("The process `snowMelt` can't be `NULL`, when `atmosSnow` is not `NULL`.")
  
  # PROCESS ----------
  idx_process <- which(!(process_method == "NULL"))
  names(idx_process) <- names(process_vari[idx_process])
  
  idx_module <- match(process_method[idx_process], module_name)
  names(idx_module) <- names(process_vari[idx_process])
  
  
  
  
  
  ## process lines ---------
  lines_process_select <- paste0(process_vari[idx_process], " = ", lines_process[idx_module], ";")
  names(lines_process_select) <- names(process_vari[idx_process])
  
  
  ## fix evatransActual ------------
  if(process_method["evatransLand"] != "NULL") {
    
    if(idx_module["evatransLand"] == idx_module["evatransSoil"]) {
      lines_process_select["evatransLand"] <- lines_process_select["evatransLand"] |> str_replace("param_evatrans_", "param_evatransLand_")
    }
    
    process_after["evatransSoil"] <- paste0(process_after["evatransSoil"], "\nland_water_mm = atmos_precipitation_mm;")
    lines_process_select["evatransLand"] <- lines_process_select["evatransLand"] |> str_replace("water_mm", "land_interceptWater_mm") |> str_replace("capacity_mm", "land_interceptCapacity_mm")
  }
  lines_process_select["evatransSoil"] <- lines_process_select["evatransSoil"] |> str_replace("water_mm", "soil_water_mm") |> str_replace("capacity_mm", "soil_capacity_mm")
  
  
  ## fix confluen ------------
  lines_process_select["confluenLand"] <- lines_process_select["confluenLand"] |> str_replace_all("confluen_", "confluenLand_")
  lines_process_select["confluenGround"] <- lines_process_select["confluenGround"] |> str_replace_all("confluen_", "confluenGround_")
  if (process_method["confluenSoil"] != "NULL")   lines_process_select["confluenSoil"] <- lines_process_select["confluenSoil"] |> str_replace_all("confluen_", "confluenSoil_")
  if (process_method["intercep"] == "NULL" & process_method["atmosSnow"] == "NULL") {
    lines_process_select["infilt"] <- lines_process_select["infilt"] |> str_replace_all("land_water_mm", "atmos_precipitation_mm")
    process_after["infilt"] <- process_after["infilt"] |> str_replace_all("land_water_mm", "atmos_precipitation_mm")
  } 
  
  # VARIABLE --------------
  argu_select <- str_split_fixed(lines_process_select, "\\(", 2)[,2] |>
    str_remove("\\);") |> str_split(", ") |>
    reduce(c) |> unique()
  if(process_method["intercep"] == "NULL" & process_method["atmosSnow"] == "NULL") argu_select <- c(argu_select, "atmos_precipitation_mm")
  if(process_method["evatransPotential"] == "NULL" ) {
    vari_boundary <- c("atmos_potentialEvatrans_mm", vari_boundary)
  }
  if(process_method["atmosSnow"] == "NULL" & process_method["snowMelt"] != "NULL") {
    vari_boundary <- c("atmos_snow_mm", vari_boundary)
  }
  
  
  argu_select <- unique(argu_select)
  argu_matrix <- argu_select[argu_select %in% vari_boundary] |> sort()
  param_ori <- argu_select[str_which(argu_select, "^param_")] |> sort()
  argu_param <- c("confluenLand_responseTime_TS", "confluenGround_responseTime_TS", param_ori) |> sort()
  if (process_method["confluenSoil"] != "NULL") argu_param <- c(argu_param, "confluenSoil_responseTime_TS") |> sort()
  vari_declare_matrxi <- c("land_runoff_mm", "ground_baseflow_mm", "confluen_streamflow_mm")
  if (process_method["confluenSoil"] != "NULL") vari_declare_matrxi <- c(vari_declare_matrxi, "soil_interflow_mm")
  if (process_method["intercep"] == "NULL" & process_method["atmosSnow"] == "NULL") {
    vari_declare_vector <- process_vari[idx_process] |> unique() |> sort()
  } else {
    vari_declare_vector <- process_vari[idx_process] |> c("land_water_mm") |> unique() |> sort()
  }
  
  
  argu_vector <- argu_select[!(argu_select %in% c(argu_matrix, argu_param, vari_declare_matrxi, vari_declare_vector))] |> sort()
  vari_matrix <- c(argu_matrix, vari_declare_matrxi)
  
  # BUILD UP -----------
  ## header --------
  lines_head <- paste0("#include <Rcpp.h>
// [[Rcpp::depends(EDCHM)]]
#include <EDCHM.h>
using namespace Rcpp;
using namespace EDCHM;
// [[Rcpp::export]]")
  
  lines_head_full <- paste0("#include \"00utilis.h\"
//' @rdname full_modells
//' @export
// [[Rcpp::export]]")
  
  ## argument ------------
  lines_argu <- paste0("\n\nNumericMatrix EDCHM_", name_model,
                       "(\nint n_time, \nint n_spat,\n",
                       c(paste0("NumericMatrix ", argu_matrix),
                         paste0("NumericVector ", c(argu_vector, argu_param))) |> paste0(collapse = ", \n"), "\n)\n{\n")
  
  ## declare ------------
  lines_declare_vector <- paste0("NumericVector ", paste0(vari_declare_vector[!(vari_declare_vector %in% c("soil_interflow_mm", "land_runoff_mm", "ground_baseflow_mm", "confluen_streamflow_mm"))], collapse = ", "), ";")
  # lines_declare_matrix <- paste0("NumericMatrix ", paste0(c("land_runoff_mm(n_time, n_spat)", "ground_baseflow_mm(n_time, n_spat)", "confluen_streamflow_mm(n_time, n_spat)"), collapse = ", "), ";")
  lines_declare_matrix <- paste0("NumericMatrix ", paste0(paste0(vari_declare_matrxi, "(n_time, n_spat)"), collapse = ", "), ";")
  
  ## time loop ------------
  lines_for_i <- "for (int i= 0; i < n_time; i++) {\n"
  idx_process_i <- (1:length(idx_process))[-str_which(names(process_vari[idx_process]), "^confluen")]
  lines_process_i <- paste0(lines_process_select[idx_process_i], "\n", process_after[idx_process[idx_process_i]], "\n")
  
  for (i in vari_matrix) {
    lines_process_i <- str_replace_all(lines_process_i, i, paste0(i, "\\(i, _\\)"))
  }
  
  ## spat loop ------------
  lines_for_j <- "}\nfor (int j= 0; j < n_spat; j++) {"
  idx_process_j <- str_which(names(process_vari[idx_process]), "^confluen")
  
  
  lines_process_j <- lines_process_select[idx_process_j] |>
    str_replace_all("_TS", "_TS(j)") |> str_replace_all("nas_n(?=[\\)|,])", "nas_n(j)") |> str_replace_all("kel_k(?=[\\)|,])", "kel_k(j)")
  if (process_method["confluenSoil"] != "NULL") {
    lines_end <- "\nconfluen_streamflow_mm(_, j) = confluen_IUH3S(land_runoff_mm(_, j), soil_interflow_mm(_, j), ground_baseflow_mm(_, j), confluenLand_iuh_1, confluenSoil_iuh_1, confluenGround_iuh_1);\n}\nreturn confluen_streamflow_mm;\n}"
    
  } else {
    lines_end <- "\nconfluen_streamflow_mm(_, j) = confluen_IUH2S(land_runoff_mm(_, j), ground_baseflow_mm(_, j), confluenLand_iuh_1, confluenGround_iuh_1);\n}\nreturn confluen_streamflow_mm;\n}"
    
  }
  
  
  # EXPORT -----------
  # write(c(lines_head,
  #         lines_argu,
  #         lines_declare_vector,
  #         lines_declare_matrix,
  #         lines_for_i,
  #         lines_process_i,
  #         lines_for_j,
  #         lines_process_j,
  #         lines_end),
  #       file.path(path_model, paste0("EDCHM_", name_model, ".cpp")))
  line_Model <- c(lines_head,
                  lines_argu,
                  lines_declare_vector,
                  lines_declare_matrix,
                  lines_for_i,
                  lines_process_i,
                  lines_for_j,
                  lines_process_j,
                  lines_end)
  
  # Param RANGE --------------------
  # lines_all_parameter <- readLines("E:\\Kan_Lei\\PACKAGE\\EDCHM\\R/all_parameter.R")
  # idx_Param <- str_which(lines_all_parameter, "@param")
  # lines_parameter_range <- (lines_all_parameter[idx_Param] |> str_split_fixed("[<|>]", 3))[,2]
  # lines_parameter_Name <- (lines_all_parameter[idx_Param] |> str_split_fixed("[<]", 3))[,1] |> str_remove("#' @param") |> str_remove_all(" ")
  # Paramrange <- read.csv(text = lines_parameter_range, header = F)
  # colnames(Paramrange) <- c("min", "max")
  # rownames(Paramrange) <- lines_parameter_Name
  
  
  param_ori_ori <- param_ori |>  str_replace("param_confluenGround", "param_confluen") |> 
    str_replace("param_confluenSoil", "param_confluen") |>  
    str_replace("param_confluenLand", "param_confluen") |>  
    str_replace("param_evatransLand_", "param_evatrans_")
  # idx_param_select <- map(param_ori_ori, \(x) str_which(lines_all_parameter, x)) |> unlist()
  # lines_parameter_range <- (lines_all_parameter[idx_param_select] |> str_split_fixed("[<|>]", 3))[,2]
  # df_range <- read.csv(text = lines_parameter_range, header = F)
  # colnames(df_range) <- c("min", "max")
  # rownames(df_range) <- param_ori
  df_range <- Paramrange[param_ori_ori,]
  rownames(df_range) <- param_ori
  
  
  if(is.null(path_model)) {
    return(list(code_Model = line_Model, range_Parameter = df_range))
  } else {
    write(line_Model, file.path(path_model, paste0("EDCHM_", name_model, ".cpp")))
    return(df_range)
  }
}
