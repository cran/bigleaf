#############################
#### Surface conditions  ####
#############################

#' Big-Leaf Surface Conditions
#' 
#' @description Calculates meteorological conditions at the big-leaf surface
#'              by inverting bulk transfer equations for water, energy, and carbon
#'              fluxes.
#' 
#' @param data             Data.frame or matrix containing all required input variables
#' @param Tair             Air temperature (deg C)
#' @param pressure         Atmospheric pressure (kPa)
#' @param H                Sensible heat flux (W m-2)
#' @param LE               Latent heat flux (W m-2)
#' @param VPD              Vapor pressure deficit (kPa)
#' @param Ga               Aerodynamic conductance for heat/water vapor (m s-1)
#' @param calc.surface.CO2 Calculate surface CO2 concentration? Defaults to \code{FALSE}.
#' @param Ca               Atmospheric CO2 concentration (mol mol-1). Required if \code{calc.surface.CO2 = TRUE}.
#' @param NEE              Net ecosystem exchange (umol m-2 s-1). Required if \code{calc.surface.CO2 = TRUE}.
#' @param Ga_CO2           Aerodynamic conductance for CO2 (m s-1). Required if \code{calc.surface.CO2 = TRUE}.
#' @param Esat.formula     Optional: formula to be used for the calculation of esat and the slope of esat.
#'                         One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}.
#'                         See \code{\link{Esat.slope}}.           
#' @param constants        cp - specific heat of air for constant pressure (J K-1 kg-1) \cr 
#'                         eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                         Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' @details Canopy surface temperature and humidity are calculated by inverting bulk transfer equations of
#'          sensible and latent heat, respectively. 'Canopy surface' in this case refers to 
#'          the surface of the big-leaf (i.e. at height d + z0h; the apparent sink of sensible heat and water vapor).
#'          Aerodynamic canopy surface temperature is given by:
#'          
#'            \deqn{Tsurf = Tair + H / (\rho * cp * Ga)}
#'          
#'          where \eqn{\rho} is air density (kg m-3). 
#'          Vapor pressure at the canopy surface is:
#'          
#'            \deqn{esurf = e + (LE * \gamma)/(Ga * \rho * cp)}
#'          
#'          where \eqn{\gamma} is the psychrometric constant (kPa K-1).
#'          Vapor pressure deficit (VPD) at the canopy surface is calculated as:
#'          
#'            \deqn{VPD_surf = Esat_surf - esurf}
#'            
#'          CO2 concentration at the canopy surface is given by:
#'          
#'            \deqn{Ca_surf = Ca + NEE / Ga_CO2}
#'          
#'          Note that Ga is assumed to be equal for water vapor and sensible heat.
#'          Ga is further assumed to be the inverse of the sum of the turbulent part
#'          and the canopy boundary layer conductance (1/Ga = 1/Ga_m + 1/Gb; 
#'          see \code{\link{aerodynamic.conductance}}). Ga_CO2, the aerodynamic conductance
#'          for CO2 is also calculated by \code{\link{aerodynamic.conductance}}.
#'          If Ga is replaced by Ga_m (i.e. only the turbulent conductance part), 
#'          the results of the functions represent conditions outside the canopy
#'          boundary layer, i.e. in the canopy airspace.
#' 
#' @note The following sign convention for NEE is employed (relevant if 
#'       \code{calc.surface.CO2 = TRUE}): 
#'       negative values of NEE denote net CO2 uptake by the ecosystem.
#' 
#' @return a data.frame with the following columns:
#'         \item{Tsurf}{Surface temperature (deg C)} \cr
#'         \item{esat_surf}{Saturation vapor pressure at the surface (kPa)} \cr
#'         \item{esurf}{vapor pressure at the surface (kPa)} \cr
#'         \item{VPD_surf}{vapor pressure deficit at the surface (kPa)} \cr
#'         \item{qsurf}{specific humidity at the surface (kg kg-1)} \cr
#'         \item{rH_surf}{relative humidity at the surface (-)} \cr
#'         \item{Ca_surf}{CO2 concentration at the surface (umol mol-1)}             
#'         
#' @examples
#' # calculate surface temperature, water vapor, VPD etc. at the surface
#' # for a given temperature and turbulent fluxes, and under different 
#' # aerodynamic conductance.
#' surface.conditions(Tair=25,pressure=100,LE=100,H=200,VPD=1.2,Ga=c(0.02,0.05,0.1)) 
#'          
#' # now calculate also surface CO2 concentration
#' surface.conditions(Tair=25,pressure=100,LE=100,H=200,VPD=1.2,Ga=c(0.02,0.05,0.1),
#'                    Ca=400,Ga_CO2=c(0.02,0.05,0.1),NEE=-20,calc.surface.CO2=TRUE)
#'                    
#' @references Knauer, J. et al., 2018: Towards physiologically meaningful water-use efficiency estimates
#'             from eddy covariance data. Global Change Biology 24, 694-710.
#'             
#'             Blanken, P.D. & Black, T.A., 2004: The canopy conductance of a boreal aspen forest,
#'             Prince Albert National Park, Canada. Hydrological Processes 18, 1561-1578.
#'             
#'             Shuttleworth, W. J., Wallace, J.S., 1985: Evaporation from sparse crops-
#'             an energy combination theory. Quart. J. R. Met. Soc. 111, 839-855.
#'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
#' @export 
surface.conditions <- function(data,Tair="Tair",pressure="pressure",LE="LE",H="H",
                               VPD="VPD",Ga="Ga_h",calc.surface.CO2=FALSE,Ca="Ca",Ga_CO2="Ga_CO2",
                               NEE="NEE",Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                               constants=bigleaf.constants()){
  
  check.input(data,list(Tair,pressure,LE,H,VPD,Ga))
  
  rho   <- air.density(Tair,pressure,constants)
  gamma <- psychrometric.constant(Tair,pressure,constants)
  
  # 1) Temperature
  Tsurf <- Tair + H / (rho * constants$cp * Ga)
  
  # 2) Humidity
  esat      <- Esat.slope(Tair,Esat.formula,constants)[,"Esat"]
  e         <- esat - VPD
  esat_surf <- Esat.slope(Tsurf,Esat.formula,constants)[,"Esat"]
  esurf     <- e + (LE * gamma)/(Ga * rho * constants$cp)
  VPD_surf  <- pmax(esat_surf - esurf,0)
  qsurf     <- VPD.to.q(VPD_surf,Tsurf,pressure,Esat.formula,constants)
  rH_surf   <- VPD.to.rH(VPD_surf,Tsurf,Esat.formula)
  
  # 3) CO2 concentration
  if (calc.surface.CO2){
    check.input(data,Ca,NEE,Ga_CO2)
    Ca_surf <- surface.CO2(Ca,NEE,Ga_CO2,Tair,pressure)
  } else {
    Ca_surf <- rep(NA_integer_,length(Tair))
  }
  
  return(data.frame(Tsurf,esat_surf,esurf,VPD_surf,qsurf,rH_surf,Ca_surf))
}



#' CO2 Concentration at the Canopy Surface
#'
#' @description the CO2 concentration at the canopy surface derived from net ecosystem
#'              CO2 exchange and measured atmospheric CO2 concentration.
#'              
#' @param Ca       Atmospheric CO2 concentration (umol mol-1)
#' @param NEE      Net ecosystem exchange (umol CO2 m-2 s-1)
#' @param Ga_CO2   Aerodynamic conductance for CO2 (m s-1)
#' @param Tair     Air temperature (degC)
#' @param pressure Atmospheric pressure (kPa)
#' 
#' @details CO2 concentration at the canopy surface is calculated as:
#' 
#'          \deqn{Ca_surf = Ca + NEE / Ga_CO2}
#'        
#'        Note that this equation can be used for any gas measured (with NEE
#'        replaced by the net exchange of the respective gas and Ga_CO2 by the Ga of 
#'        that gas).
#' 
#' @note the following sign convention is employed: negative values of NEE denote
#'       net CO2 uptake by the ecosystem.
#' 
#' @return \item{Ca_surf -}{CO2 concentration at the canopy surface (umol mol-1)}
#' 
#' @examples 
#' surface.CO2(Ca=400,NEE=-30,Ga_CO2=0.05,Tair=25,pressure=100)
#' 
#' @export
surface.CO2 <- function(Ca,NEE,Ga_CO2,Tair,pressure){
  
  Ga_CO2 <- ms.to.mol(Ga_CO2,Tair,pressure)
  
  Ca_surf <- Ca + NEE/Ga_CO2
  
  return(Ca_surf)
  
}



#' Radiometric Surface Temperature
#' 
#' @description Radiometric surface temperature from longwave radiation
#'              measurements.
#'              
#' @param data        Data.frame or matrix containing all required input variables            
#' @param LW_up       Longwave upward radiation (W m-2)
#' @param LW_down     Longwave downward radiation (W m-2)
#' @param emissivity  Emissivity of the surface (-)
#' @param constants   sigma - Stefan-Boltzmann constant (W m-2 K-4) \cr
#'                    Kelvin - conversion degree Celsius to Kelvin 
#' 
#' @details Radiometric surface temperature (Trad) is calculated as:
#' 
#'            \deqn{Trad = ((LW_up - (1 - \epsilon)*LW_down) / (\sigma \epsilon))^(1/4)}   
#' 
#' @return a data.frame with the following columns:
#'         \item{Trad_K}{Radiometric surface temperature (K)} \cr
#'         \item{Trad_degC}{Radiometric surface temperature (degC)} 
#' 
#' @examples 
#' # determine radiometric surface temperature for the site DE-Tha in June 2014 
#' # assuming an emissivity of 0.98.
#' # (Note that variable 'LW_down' was only included for the DE-Tha example dataset
#' # and not for the others due restrictions on file size) 
#' Trad <- radiometric.surface.temp(DE_Tha_Jun_2014,emissivity=0.98)
#' summary(Trad)
#' 
#' @references Wang, W., Liang, S., Meyers, T. 2008: Validating MODIS land surface
#'             temperature products using long-term nighttime ground measurements.
#'             Remote Sensing of Environment 112, 623-635.
#' 
#' @export
radiometric.surface.temp <- function(data,LW_up="LW_up",LW_down="LW_down",
                                     emissivity,constants=bigleaf.constants()){
  
  check.input(data,list(LW_up,LW_down))
  
  Trad_K    <- ((LW_up - (1 - emissivity)*LW_down) / (constants$sigma * emissivity))^(1/4)
  Trad_degC <- Trad_K - constants$Kelvin
  
  return(data.frame(Trad_K,Trad_degC))
}


### can be added at some point, but the emissivity values seem to be fairly low
# emissivity.from.albedo <- function(albedo){
#   
#   emissivity <- -0.16*albedo + 0.99
#   
#   return(emissivity)
# }

