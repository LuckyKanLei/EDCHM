#' Calibrate algorithms
#' @name cali
#' @description
#' # DDS \insertCite{DDS_Tolson_2007}{EDCHM}
#' dynamically dimensioned search (DDS) calibration algorithm
#' @param fitness fitness function
#' @param x_Min,x_Max,x_Init minimal, maximal initial parameter
#' @param max_iter maximal number of iteration
#' @param ... other parameters
#' @param r parameter for algorithm
#' @importFrom stats rbinom rnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @import hydroGOF
#' @export
cali_DDS <- function(fitness, x_Min, x_Max, x_Init = NA, max_iter = 100, r = 0.2, ...){
  ## S.1 Set initial parameters ####
  if (any(is.na(x_Init))) x_Init <- (x_Min + x_Max) / 2
  
  x_Best <- x_Init
  
  ## S.2 Evaluate initial ####
  y_Init <- fitness(x_Init, ...)
  y_Best <- y_Init
  
  ## S.3 Select peturb J of the D decision
  n_x <- length(x_Min)
  P_i <- 1 - log(1:max_iter) / log(max_iter)
  randm_Para <- apply(as.matrix(P_i), 1, function(x) as.logical(rbinom(n_x, 1, x)))
  idx_null <- which(colSums(randm_Para) == 0)
  idx_new4null <- sample(1:n_x, length(idx_null), replace = T)
  randm_Para[idx_new4null, idx_null] <- TRUE
  lst_Cali_x <- apply(randm_Para, 2, which)
  
  ## S.4 Peturb each entry by N(0,1)*r(x_max - x_min) reflecting
  sigma_ <- x_Max - x_Min
  
  bar_progress <- txtProgressBar(style = 3)
  for(i in 2:max_iter){
    setTxtProgressBar(bar_progress, i/max_iter, title = paste(i/max_iter,"% of Calibration"))
    
    x_New <- x_Best
    idx <- lst_Cali_x[[i]]
    N_01 <- rnorm(n_x)
    x_New0 <- (x_Best + r * N_01 * sigma_)
    x_New1 <- cbind(2 * x_Min - x_New0, x_Min) |> apply(1, max)
    x_New2 <- cbind(2 * x_Max - x_New0, x_Max) |> apply(1, min)
    x_New0[x_New0 < x_Min] <- x_New1[x_New0 < x_Min]
    x_New0[x_New0 > x_Max] <- x_New2[x_New0 > x_Max]
    x_New[idx] <- x_New0[idx]
    ## S.5 Evaluate objective function
    y_New <- fitness(x_New, ...)
    if(y_New < y_Best){
      x_Best <- x_New
      y_Best <- y_New
    }
  }
  close(bar_progress)
  list(x_Best = x_Best, y_Best = y_Best)
}



#' @rdname cali
#' @description
#' # UVS
#' Univariable search algorithmus
#' @param max_n_loop maximal number of loops
#' @param n_1_search maximal number in one loop
#' @export
cali_UVS <- function(fitness, x_Min, x_Max, x_Init = NA, max_n_loop = 10, n_1_search = 8, ...){
  ## S.1 Set initial parameters ####
  if (any(is.na(x_Init))) x_Init <- (x_Min + x_Max) / 2
  x_Best <- x_Init
  
  ## S.2 Evaluate initial ####
  y_Init <- fitness(x_Init, ...)
  y_Best <- y_Init
  
  
  ## loop for every variable -----------
  n_param <- length(x_Min)
  bar_progress <- txtProgressBar(style = 3)
  for (j in 1:max_n_loop) {
    for (i in 1:n_param) {
      best_i <- uvs_1(fitness, x_Min[i], x_Max[i], i, x_Best, n_1_search, ...)
      
      if (best_i$y_Best < y_Best) {
        x_Best[i] <- best_i$x_Best
        y_Best <- best_i$y_Best
      }
      
    }
    
  }
  close(bar_progress)
  
  list(x_Best = x_Best, y_Best = y_Best)
  
}







uvs_1 <- function(fitness, x_min, x_max, i_x, x_Best, n_1_search = 8, ...){
  x_Best <- x_Best
  x_Best[i_x] <- x_min
  y_left <- fitness(x_Best, ...)
  
  
  x_Best[i_x] <- x_max
  y_right <- fitness(x_Best, ...)
  
  x_left <- x_min
  x_right <- x_max
  
  
  for(i in 1:n_1_search){
    
    
    
    if(y_right < y_left){
      x_left <- x_right + 0.618*(x_left - x_right)
      
      x_Best[i_x] <- x_left
      y_left <- fitness(x_Best, ...)
      
    } else {
      x_right <- x_left + 0.618*(x_right - x_left)
      
      x_Best[i_x] <- x_right
      y_right <- fitness(x_Best, ...)
      
    }
  }
  if (y_left < y_right) {
    return(list(x_Best = x_left, y_Best = y_left))
  } else {
    return(list(x_Best = x_right, y_Best = y_right))
  }
}
