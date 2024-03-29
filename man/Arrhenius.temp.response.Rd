% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigleaf_physiology.r
\name{Arrhenius.temp.response}
\alias{Arrhenius.temp.response}
\title{(Modified) Arrhenius Temperature Response Function}
\usage{
Arrhenius.temp.response(
  param,
  Temp,
  Ha,
  Hd,
  dS,
  constants = bigleaf.constants()
)
}
\arguments{
\item{param}{Parameter measured at measurement temperature (umol m-2 s-1)}

\item{Temp}{Measurement temperature (degC)}

\item{Ha}{Activation energy for param (kJ mol-1)}

\item{Hd}{Deactivation energy for param (kJ mol-1)}

\item{dS}{Entropy term for param (kJ mol-1 K-1)}

\item{constants}{Kelvin - conversion degree Celsius to Kelvin \cr
Rgas - universal gas constant (J mol-1 K-1) \cr
kJ2J - conversion kilojoule (kJ) to joule (J)}
}
\value{
param25 - value of the input parameter at the reference temperature of 25degC (umol m-2 s-1)
}
\description{
(Modified) Arrhenius function describing
             the temperature response of biochemical parameters.
}
\details{
The function returns the biochemical rate at a reference
         temperature of 25degC given a predefined temperature response function.
         This temperature response is given by a modified form of the Arrhenius
         function:

            \deqn{param25 = param / 
                        ( exp(Ha * (Temp - Tref) / (Tref*Rgas*Temp)) *
                        (1 + exp((Tref*dS - Hd) / (Tref * Rgas))) /
                        (1 + exp((Temp*dS - Hd) / (Temp * Rgas)))
                        )
                 }
                 
         where param is the value/rate of the parameter at measurement temperature,
         Temp is temperature in K, Tref is reference temperature (298.15K), and Rgas
         is the universal gas constant (8.314 J K-1 mol-1). Ha is the activation
         energy (kJ mol-1), Hd is the deactivation energy (kJ mol-1), and dS the
         entropy term (kJ mol-1 K-1) of the respective parameter.
         
         If either Hd or dS or both are not provided, the equation above reduces
         to the first term (i.e. the common Arrhenius equation without the deactivation
         term.)
}
\references{
Johnson F.H., Eyring H., Williams R.W. 1942: 
            The nature of enzyme inhibitions in bacterial luminescence: sulfanilamide,
            urethane, temperature and pressure. Journal of cellular and comparative
            physiology 20, 247-268.

            Kattge J., Knorr W., 2007: Temperature acclimation in a biochemical
            model of photosynthesis: a reanalysis of data from 36 species.
            Plant, Cell and Environment 30, 1176-1190.
}
