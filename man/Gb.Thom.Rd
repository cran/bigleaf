% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boundary_layer_conductance.r
\name{Gb.Thom}
\alias{Gb.Thom}
\title{Boundary Layer Conductance according to Thom 1972}
\usage{
Gb.Thom(ustar, Sc = NULL, Sc_name = NULL, constants = bigleaf.constants())
}
\arguments{
\item{ustar}{Friction velocity (m s-1)}

\item{Sc}{Optional: Schmidt number of additional quantities to be calculated}

\item{Sc_name}{Optional: Name of the additional quantities, has to be of same length than 
\code{Sc_name}}

\item{constants}{k - von-Karman constant \cr
Sc_CO2 - Schmidt number for CO2 \cr 
Pr - Prandtl number (if \code{Sc} is provided)}
}
\value{
a data.frame with the following columns:
 \item{Gb_h}{Boundary layer conductance for heat transfer (m s-1)}
 \item{Rb_h}{Boundary layer resistance for heat transfer (s m-1)}
 \item{kB_h}{kB-1 parameter for heat transfer}
 \item{Gb_Sc_name}{Boundary layer conductance for \code{Sc_name} (m s-1). Only added if \code{Sc_name} and 
                   \code{Sc_name} are provided}
}
\description{
An empirical formulation for the canopy boundary layer conductance 
             for heat transfer based on a simple ustar dependency.
}
\details{
The empirical equation for Rb suggested by Thom 1972 is:
 
     \deqn{Rb = 6.2ustar^-0.67}
 
   Gb (=1/Rb) for water vapor and heat are assumed to be equal in this package.
   Gb for other quantities x is calculated as (Hicks et al. 1987):
 
     \deqn{Gb_x = Gb / (Sc_x / Pr)^0.67}
 
 where Sc_x is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).
}
\examples{
Gb.Thom(seq(0.1,1.4,0.1))

## calculate Gb for SO2 as well
Gb.Thom(seq(0.1,1.4,0.1),Sc=1.25,Sc_name="SO2")

}
\references{
Thom, A., 1972: Momentum, mass and heat exchange of vegetation.
            Quarterly Journal of the Royal Meteorological Society 98, 124-134.
            
            Hicks, B.B., Baldocchi, D.D., Meyers, T.P., Hosker, J.R., Matt, D.R., 1987:
            A preliminary multiple resistance routine for deriving dry deposition velocities
            from measured quantities. Water, Air, and Soil Pollution 36, 311-330.
}
\seealso{
\code{\link{Gb.Choudhury}}, \code{\link{Gb.Su}}, \code{\link{aerodynamic.conductance}}
}
