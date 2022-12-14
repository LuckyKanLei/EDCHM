% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{atmosSnow}
\alias{atmosSnow}
\alias{atmosSnow_ThresholdT}
\alias{atmosSnow_UBC}
\title{caculate \strong{snowfall}}
\usage{
atmosSnow_ThresholdT(
  atmos_precipitation_mm,
  atmos_temperature_Cel,
  param_atmos_thr_Ts
)

atmosSnow_UBC(
  atmos_precipitation_mm,
  atmos_temperature_Cel,
  param_atmos_ubc_A0FORM
)
}
\arguments{
\item{atmos_precipitation_mm}{(mm/m2/TS) precipitaion volum}

\item{atmos_temperature_Cel}{(Cel) the average air temperature in the time phase}

\item{param_atmos_thr_Ts}{<-1, 3> (Cel) threshold air temperature that snow, parameter for \code{\link[=atmosSnow_ThresholdT]{atmosSnow_ThresholdT()}}}

\item{param_atmos_ubc_A0FORM}{<0.01, 3> (Cel) threshold air temperature that snow, it can not equal or small than 0, parameter for \code{\link[=atmosSnow_UBC]{atmosSnow_UBC()}}}
}
\value{
atmos_snow_mm (mm/m2/TS) snowfall volume
}
\description{
\loadmathjax

Under the concept of the conceptual HM, the amount of snowfall is always calculated by the temperature \mjseqn{T} and the precipitation \mjseqn{P} availability.
The proportion of snowfall is always determined by the air temperature.

So we can give the function from:

\mjsdeqn{P_s = f_{atmosSnow}(D_{atms})}

to:

\mjsdeqn{P_s = f_{atmosSnow}(P, T) = k^*P}
\mjsdeqn{0 \leq k^* \leq 1}
where
\itemize{
\item \mjseqn{P} is \code{atmos_precpitation_mm}
\item \mjseqn{T} is \code{atmos_teperature_Cel}
\item \mjseqn{k^*} is estimated portion
}

Then the different \code{atmosSnow} methods will estimate the portion \mjseqn{k^*}.

The output density distribution from 2 methods:

\if{html}{\figure{mdl_atmosSnow.svg}}
\if{latex}{\figure{mdl_atmosSnow.pdf}{options: width=140mm}}
}
\section{\strong{_ThresholdT}:}{
\if{html}{\figure{mdl_atmosSnow_thr.svg}}
\if{latex}{\figure{mdl_atmosSnow_thr.pdf}{options: width=140mm}}

Only a temperature is as the threshold defined, so estimate the portion \mjseqn{k^*} as:
\mjsdeqn{k^{*}=1, \quad T \leq T_s}
where
\itemize{
\item \mjseqn{T_s} is \code{param_atmos_thr_Ts}
}
}

\section{\strong{_UBC} \insertCite{UBC_Quick_1977}{EDCHM}:}{
\if{html}{\figure{mdl_atmosSnow_ubc.svg}}
\if{latex}{\figure{mdl_atmosSnow_ubc.pdf}{options: width=140mm}}

estimate the portion \mjseqn{k^*}{} as:
\mjsdeqn{k^* = 1- \frac{T}{T_0}}
\mjsdeqn{k^* \geq 0}
where
\itemize{
\item \mjseqn{T_0} is \code{param_atmos_ubc_A0FORM}
}
}

\references{
\insertAllCited{}
}
