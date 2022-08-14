#' create **IUH** (Instant Unit Graphy)
#' @name confluenIUH
#' @param confluen_resposeTime (TS) input water volum in every routeline
#' @return IUH (list of num vector) 
#' @export
confluenIUH_GR4J1 <- function(
  confluen_resposeTime
)
{
  t_max <- ceiling(confluen_resposeTime)
  SH_1 <- (1:t_max / confluen_resposeTime)^2.5
  SH_1[t_max] <- 1
  c(SH_1[1], diff(SH_1))
}



#' @rdname confluenIUH
#' @export
confluenIUH_GR4J2 <- function(
  confluen_resposeTime
)
{
  t_max_1 <- ceiling(confluen_resposeTime)
  t_max_2 <- ceiling(2 * confluen_resposeTime)
  SH_2_1 <- .5 * (1:(t_max_1 - 1) / confluen_resposeTime)^2.5
  SH_2_2 <- 1 - .5 * (2 - t_max_1:(t_max_2 - 1) / confluen_resposeTime)^2.5
  SH_2 <- c(SH_2_1, SH_2_2, 1)
  c(SH_2_1[1], diff(SH_2))
}
