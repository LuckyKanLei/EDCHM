% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{baseflow}
\alias{baseflow}
\alias{baseflow_GR4J}
\alias{baseflow_GR4Jfix}
\alias{baseflow_SupplyRatio}
\alias{baseflow_SupplyPow}
\alias{baseflow_MaxPow}
\alias{baseflow_ThreshPow}
\alias{baseflow_Arno}
\title{\strong{baseflow}}
\usage{
baseflow_GR4J(ground_water_mm, ground_capacity_mm)

baseflow_GR4Jfix(ground_water_mm, ground_capacity_mm, param_baseflow_grf_gamma)

baseflow_SupplyRatio(ground_water_mm, param_baseflow_sur_k)

baseflow_SupplyPow(
  ground_water_mm,
  param_baseflow_sup_k,
  param_baseflow_sup_gamma
)

baseflow_MaxPow(
  ground_water_mm,
  ground_capacity_mm,
  ground_potentialBaseflow_mm,
  param_baseflow_map_gamma
)

baseflow_ThreshPow(
  ground_water_mm,
  ground_capacity_mm,
  ground_potentialBaseflow_mm,
  param_baseflow_thp_thresh,
  param_baseflow_thp_gamma
)

baseflow_Arno(
  ground_water_mm,
  ground_capacity_mm,
  ground_potentialBaseflow_mm,
  param_baseflow_arn_thresh,
  param_baseflow_arn_k
)
}
\arguments{
\item{ground_water_mm}{(mm/m2/TS) water volume in \code{groundLy}}

\item{ground_capacity_mm}{(mm/m2) water storage capacity in \code{groundLy}}

\item{param_baseflow_grf_gamma}{<2, 7> exponential parameter for \code{\link[=baseflow_GR4Jfix]{baseflow_GR4Jfix()}}}

\item{param_baseflow_sur_k}{<0.01, 1> coefficient parameter for \code{\link[=baseflow_SupplyRatio]{baseflow_SupplyRatio()}}}

\item{param_baseflow_sup_k}{<0.01, 1> coefficient parameter for \code{\link[=baseflow_SupplyPow]{baseflow_SupplyPow()}}}

\item{param_baseflow_sup_gamma}{<0, 1> exponential parameter for \code{\link[=baseflow_SupplyPow]{baseflow_SupplyPow()}}}

\item{ground_potentialBaseflow_mm}{<0.01, 7> (mm/m2/TS) \strong{potential} baseflow}

\item{param_baseflow_map_gamma}{<0.1, 5> exponential parameter for \code{\link[=baseflow_MaxPow]{baseflow_MaxPow()}}}

\item{param_baseflow_thp_thresh}{<0.1, 0.9> coefficient parameter for \code{\link[=baseflow_ThreshPow]{baseflow_ThreshPow()}}}

\item{param_baseflow_thp_gamma}{<0.1, 5> exponential parameter for \code{\link[=baseflow_ThreshPow]{baseflow_ThreshPow()}}}

\item{param_baseflow_arn_thresh}{<0.1, 0.9> coefficient parameter for \code{\link[=baseflow_ThreshPow]{baseflow_ThreshPow()}}}

\item{param_baseflow_arn_k}{<0.1, 1> exponential parameter for \code{\link[=baseflow_ThreshPow]{baseflow_ThreshPow()}}}
}
\value{
ground_baseflow_mm (mm/m2/TS)
}
\description{
\loadmathjax

In hydrological modeling, baseflow refers to the flow of water in rivers and streams that is sustained by the release of water from the groundwater.
Or baseflow refers to the flow of water from an aquifer or deeper soil horizon to surface water, typically due to a head gradient between fully saturated soil and stream  \insertCite{Raven_Manual_35}{EDCHM}.
It may be considered the sum of the contribution of deep groundwater exchange with a river and delayed storage  \insertCite{Raven_Manual_35}{EDCHM}.

It is always calculated (only) by the water in the ground layer \mjseqn{W_{grnd}}, which can also be treated as part of \mjseqn{W_{grnd}}.
However, the impact of other RUs (response units) on the route to the river will be ignored.

So we can give the function from:

\mjsdeqn{F_{base} = f_{baseflow}(D_{grnd})}

to:

\mjsdeqn{F_{base} = f_{baseflow}(W_{grnd}, C_{grnd}, M_{base}, ...)}
\mjsdeqn{F_{base} = k^* W_{grnd} \quad {\rm or} \quad F_{base} = k^* M_{base}}
\mjsdeqn{0 \leq k^* \leq 1}

where
\itemize{
\item \mjseqn{W_{grnd}} is \code{ground_water_mm}
\item \mjseqn{M_{base}} is \code{ground_potentialBaseflow_mm}
\item \mjseqn{C_{grnd}} is \code{ground_capacity_mm}, but not all the methods need the \mjseqn{C_{grnd}}
\item \mjseqn{k^*} is estimated ratio
}

The output density distribution from 7 methods:

\if{html}{\figure{mdl_baseflow.svg}}
\if{latex}{\figure{mdl_baseflow.pdf}{options: width=140mm}}
}
\section{\strong{_GR4J} \insertCite{GR4J_Perrin_2003}{EDCHM}:}{
\if{html}{\figure{mdl_baseflow_gr4.svg}}
\if{latex}{\figure{mdl_baseflow_gr4.pdf}{options: width=140mm}}

\mjsdeqn{F_{base} = k^* W_{grnd}}
\mjsdeqn{k^* = 1 - \left[ 1 + \left(\frac{W_{grnd}}{C_{grnd}} \right)^4 \right]^{-1/4}}
where
\itemize{
\item \mjseqn{k^*} is estimated ratio
}
}

\section{\strong{_GR4Jfix} \insertCite{GR4J_Perrin_2003}{EDCHM}:}{
\if{html}{\figure{mdl_baseflow_grf.svg}}
\if{latex}{\figure{mdl_baseflow_grf.pdf}{options: width=140mm}}

This method based on \verb{_GR4J} use a new parameter to replace the numer 4:
\mjsdeqn{F_{base} = k^* W_{grnd}}
\mjsdeqn{k^* = 1 - \left[ 1 + \left(\frac{W_{grnd}}{C_{grnd}} \right)^\gamma \right]^{-1/\gamma}}
where
\itemize{
\item \mjseqn{\gamma} is \code{param_baseflow_grf_gamma}
}
}

\section{\strong{_SupplyRatio}:}{
\if{html}{\figure{mdl_baseflow_sur.svg}}
\if{latex}{\figure{mdl_baseflow_sur.pdf}{options: width=140mm}}

\mjsdeqn{F_{base} = k W_{grnd}}
where
\itemize{
\item \mjseqn{k} is \code{param_baseflow_sur_k}
}
}

\section{\strong{_SupplyPow}:}{
\if{html}{\figure{mdl_baseflow_sup.svg}}
\if{latex}{\figure{mdl_baseflow_sup.pdf}{options: width=140mm}}

\mjsdeqn{F_{base} = k(W_{grnd})^\gamma}
where
\itemize{
\item \mjseqn{k} is \code{param_baseflow_sup_k}
\item \mjseqn{\gamma} is \code{param_baseflow_sup_gamma}
}
}

\section{\strong{_MaxPow}:}{
\if{html}{\figure{mdl_baseflow_map.svg}}
\if{latex}{\figure{mdl_baseflow_map.pdf}{options: width=140mm}}

\mjsdeqn{F_{base} = M_{base} \left(\frac{W_{grnd}}{C_{grnd}} \right)^\gamma}
where
\itemize{
\item \mjseqn{M_{base}} is \code{ground_potentialBaseflow_mm}
\item \mjseqn{\gamma} is \code{param_baseflow_map_gamma}
}
}

\section{\strong{_ThreshPow}}{
\if{html}{\figure{mdl_baseflow_thp.svg}}
\if{latex}{\figure{mdl_baseflow_thp.pdf}{options: width=140mm}}

This method based on the \verb{_MaxPow} and add the one threshold \mjseqn{\phi_b}:
\mjsdeqn{F_{base} = 0, \quad \frac{W_{grnd}}{C_{grnd}} < \phi_b}
\mjsdeqn{F_{base} = M_{base} \left(\frac{\frac{W_{grnd}}{C_{grnd}} - \phi_b}{1-\phi_b} \right)^\gamma, \quad \frac{W_{grnd}}{C_{grnd}} \geq \phi_b}
where
\itemize{
\item \mjseqn{\phi_b} is \code{param_baseflow_thp_thresh}
\item \mjseqn{\gamma} is \code{param_baseflow_thp_gamma}
}
}

\section{\strong{_Arno} \insertCite{baseflow_Arno_1991,VIC2_Liang_1994}{EDCHM}:}{
\if{html}{\figure{mdl_baseflow_arn.svg}}
\if{latex}{\figure{mdl_baseflow_arn.pdf}{options: width=140mm}}

This method has also in two cases divided by a threshold water content \mjseqn{\phi_b}:
\mjsdeqn{F_{base} = k M_{base} \frac{W_{grnd}}{C_{grnd}}, \quad \frac{W_{grnd}}{C_{grnd}} < \phi_b}
\mjsdeqn{F_{base} = k M_{base} \frac{W_{grnd}}{C_{grnd}} + (1-k) M_{base} \left(\frac{W_{grnd} - W_s}{C_{grnd} - W_s} \right)^2, \quad \frac{W_{grnd}}{C_{grnd}} \geq \phi_b}
\mjsdeqn{W_s = k C_{grnd}}
where
\itemize{
\item \mjseqn{\phi_b} is \code{param_baseflow_arn_thresh}
\item \mjseqn{k} is \code{param_baseflow_arn_k}
}
}

\references{
\insertAllCited{}
}
