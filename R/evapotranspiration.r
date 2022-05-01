###########################
#### Evapotranspiration ###
###########################

#' Potential Evapotranspiration
#' 
#' @description Potential evapotranspiration according to Priestley & Taylor 1972 or
#'              the Penman-Monteith equation with a prescribed surface conductance.
#' 
#' @param data      Data.frame or matrix containing all required variables; optional
#' @param Tair      Air temperature (degC)
#' @param pressure  Atmospheric pressure (kPa)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
#' @param VPD       Vapor pressure deficit (kPa); only used if \code{approach = "Penman-Monteith"}.
#' @param Ga        Aerodynamic conductance to heat/water vapor (m s-1); only used if \code{approach = "Penman-Monteith"}.
#' @param approach  Approach used. Either \code{"Priestley-Taylor"} (default), or \code{"Penman-Monteith"}.
#' @param alpha     Priestley-Taylor coefficient; only used if \code{approach = "Priestley-Taylor"}.
#' @param Gs_pot    Potential/maximum surface conductance (mol m-2 s-1); defaults to 0.6 mol m-2 s-1;
#'                  only used if \code{approach = "Penman-Monteith"}.
#' @param missing.G.as.NA  if \code{TRUE}, missing G are treated as \code{NA}s, otherwise set to 0. 
#' @param missing.S.as.NA  if \code{TRUE}, missing S are treated as \code{NA}s, otherwise set to 0. 
#' @param Esat.formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}.
#'                      See \code{\link{Esat.slope}}. 
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air \cr
#'                  Pa2kPa - conversion pascal (Pa) to kilopascal (kPa) \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) (only used if \code{approach = "Penman-Monteith"}) \cr
#'                  Rgas - universal gas constant (J mol-1 K-1) (only used if \code{approach = "Penman-Monteith"}) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin (only used if \code{approach = "Penman-Monteith"}) \cr
#' 
#' @details Potential evapotranspiration is calculated according to Priestley & Taylor, 1972
#'          if \code{approach = "Priestley-Taylor"} (the default):
#' 
#'            \deqn{LE_pot,PT = (\alpha * \Delta * (Rn - G - S)) / (\Delta + \gamma)}
#'
#'          \eqn{\alpha} is the Priestley-Taylor coefficient, \eqn{\Delta} is the slope 
#'          of the saturation vapor pressure curve (kPa K-1), and \eqn{\gamma} is the 
#'          psychrometric constant (kPa K-1).
#'          if \code{approach = "Penman-Monteith"}, potential evapotranspiration is calculated according
#'          to the Penman-Monteith equation:
#' 
#'          \deqn{LE_pot,PM = (\Delta * (Rn - G - S) + \rho * cp * VPD * Ga) / (\Delta + \gamma * (1 + Ga/Gs_pot)}
#'          
#'          where \eqn{\Delta} is the slope of the saturation vapor pressure curve (kPa K-1),
#'          \eqn{\rho} is the air density (kg m-3), and \eqn{\gamma} is the psychrometric constant (kPa K-1).
#'          The value of \code{Gs_pot} is typically a maximum value of Gs observed at the site, e.g. the 90th
#'          percentile of Gs within the growing season.
#'          
#' @return a data.frame with the following columns:
#'         \item{ET_pot}{Potential evapotranspiration (kg m-2 s-1)}
#'         \item{LE_pot}{Potential latent heat flux (W m-2)}
#'         
#' @note If the first argument \code{data} is provided (either a matrix or a data.frame),
#'       the following variables can be provided as character (in which case they are interpreted as
#'       the column name of \code{data}) or as numeric vectors, in which case they are taken
#'       directly for the calculations. If \code{data} is not provided, all input variables have to be
#'       numeric vectors.        
#'   
#' @references Priestley, C.H.B., Taylor, R.J., 1972: On the assessment of surface heat flux
#'             and evaporation using large-scale parameters. Monthly Weather Review 100, 81-92.  
#'             
#'             Allen, R.G., Pereira L.S., Raes D., Smith M., 1998: Crop evapotranspiration -
#'             Guidelines for computing crop water requirements - FAO Irrigation and drainage
#'             paper 56.
#'              
#'             Novick, K.A., et al. 2016: The increasing importance of atmospheric demand
#'             for ecosystem water and carbon fluxes. Nature Climate Change 6, 1023 - 1027.
#'             
#' @seealso \code{\link{surface.conductance}}
#'                                
#' @examples 
#' # Calculate potential ET of a surface that receives a net radiation of 500 Wm-2
#' # using Priestley-Taylor:
#' potential.ET(Tair=30,pressure=100,Rn=500,alpha=1.26,approach="Priestley-Taylor")    
#' 
#' # Calculate potential ET for a surface with known Gs (0.5 mol m-2 s-1) and Ga (0.1 m s-1)
#' # using Penman-Monteith:
#' LE_pot_PM <- potential.ET(Gs_pot=0.5,Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,
#'                           approach="Penman-Monteith")[,"LE_pot"]
#' LE_pot_PM
#' 
#' # now cross-check with the inverted equation
#' surface.conductance(Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,LE=LE_pot_PM)
#' @export
potential.ET <- function(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,
                         VPD="VPD",Ga="Ga_h",approach=c("Priestley-Taylor","Penman-Monteith"),
                         alpha=1.26,Gs_pot=0.6,missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                         Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                         constants=bigleaf.constants()){
  
  approach <- match.arg(approach)
  
  check.input(data,list(Tair,pressure,Rn,G,S))
  
  if(!is.null(G)){
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    cat("Ground heat flux G is not provided and set to 0.",fill=TRUE)
    G <- 0
  }
  
  if(!is.null(S)){
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    cat("Energy storage fluxes S are not provided and set to 0.",fill=TRUE)
    S <- 0
  }
  
  gamma  <- psychrometric.constant(Tair,pressure,constants)
  Delta  <- Esat.slope(Tair,Esat.formula,constants)[,"Delta"]
  
  
  if (approach == "Priestley-Taylor"){
    
    LE_pot <- (alpha * Delta * (Rn - G - S)) / (Delta + gamma)
    ET_pot <- LE.to.ET(LE_pot,Tair)
    
  } else if (approach == "Penman-Monteith"){
    
    check.input(data,list(Gs_pot,VPD,Ga))
    
    Gs_pot <- mol.to.ms(Gs_pot,Tair=Tair,pressure=pressure,constants=constants)
    rho    <- air.density(Tair,pressure,constants)
    
    LE_pot <- (Delta * (Rn - G - S) + rho * constants$cp * VPD * Ga) / 
      (Delta + gamma * (1 + Ga / Gs_pot))
    ET_pot <- LE.to.ET(LE_pot,Tair)
  }
  
  return(data.frame(ET_pot,LE_pot))
  
}




#' Reference Evapotranspiration
#' 
#' @description Reference evapotranspiration calculated from the Penman-Monteith
#'              equation with a prescribed surface conductance.
#'              This function is deprecated. Use potential.ET(...,approach="Penman-Monteith") instead.
#' 
#' @param data      Data.frame or matrix containing all required variables; optional
#' @param Gs_ref    Reference surface conductance (m s-1); defaults to 0.0143 m s-1.
#' @param Tair      Air temperature (degC)
#' @param pressure  Atmospheric pressure (kPa)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Ga        Aerodynamic conductance to heat/water vapor (m s-1)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
#' @param missing.G.as.NA  if \code{TRUE}, missing G are treated as \code{NA}s, otherwise set to 0. 
#' @param missing.S.as.NA  if \code{TRUE}, missing S are treated as \code{NA}s, otherwise set to 0. 
#' @param Esat.formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}.
#'                      See \code{\link{Esat.slope}}. 
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) (only if \code{approach = "Penman-Monteith"}) \cr
#'                  Rgas - universal gas constant (J mol-1 K-1) (only if \code{approach = "Penman-Monteith"}) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin (only if \code{approach = "Penman-Monteith"}) \cr
#' 
#' @export                            
reference.ET <- function(data,Gs_ref=0.0143,Tair="Tair",pressure="pressure",VPD="VPD",Rn="Rn",Ga="Ga_h",
                         G=NULL,S=NULL,missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                         Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                         constants=bigleaf.constants()){
  
  stop("this function is deprecated (since bigleaf version 0.6.0). For the calculation of potential ET from the Penman-Monteith equation (as formerly calculated with this function), use function potential.ET() with the argument approach='Penman-Monteith'. Note that the default value for argument 'Gs_pot' is expressed now in mol m-2 s-1 for simplicity (0.6 mol m-2 s-1).")
  
}





#' Equilibrium and Imposed Evapotranspiration
#' 
#' @description Evapotranspiration (ET) split up into imposed ET and equilibrium ET.
#' 
#' @param data      Data.frame or matrix containing all required input variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param VPD       Air vapor pressure deficit (kPa)
#' @param Gs        surface conductance to water vapor (m s-1)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
#' @param missing.G.as.NA  if \code{TRUE}, missing G are treated as \code{NA}s, otherwise set to 0. 
#' @param missing.S.as.NA  if \code{TRUE}, missing S are treated as \code{NA}s, otherwise set to 0.
#' @param Esat.formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}.
#'                      See \code{\link{Esat.slope}}. 
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                  Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#'                  
#' @details Total evapotranspiration can be written in the form (Jarvis & McNaughton 1986):
#' 
#'            \deqn{ET = \Omega ET_eq + (1 - \Omega)ET_imp}
#'          
#'          where \eqn{\Omega} is the decoupling coefficient as calculated from
#'          \code{\link{decoupling}}. \code{ET_eq} is the equilibrium evapotranspiration rate,
#'          the ET rate that would occur under uncoupled conditions, where the heat budget
#'          is dominated by radiation (when Ga -> 0):
#'          
#'            \deqn{ET_eq = (\Delta * (Rn - G - S) * \lambda) / (\Delta + \gamma)}
#'          
#'          where \eqn{\Delta} is the slope of the saturation vapor pressure curve (kPa K-1),
#'          \eqn{\lambda} is the latent heat of vaporization (J kg-1), and \eqn{\gamma}
#'          is the psychrometric constant (kPa K-1).
#'          \code{ET_imp} is the imposed evapotranspiration rate, the ET rate
#'          that would occur under fully coupled conditions (when Ga -> inf):
#'          
#'            \deqn{ET_imp = (\rho * cp * VPD * Gs * \lambda) / \gamma}
#'          
#'          where \eqn{\rho} is the air density (kg m-3).
#' 
#' @note Surface conductance (Gs) can be calculated with \code{\link{surface.conductance}}.
#'       Aerodynamic conductance (Ga) can be calculated using \code{\link{aerodynamic.conductance}}.
#'       
#' @return A data.frame with the following columns:
#'         \item{ET_eq}{Equilibrium ET (kg m-2 s-1)}
#'         \item{ET_imp}{Imposed ET (kg m-2 s-1)}
#'         \item{LE_eq}{Equilibrium LE (W m-2)}
#'         \item{LE_imp}{Imposed LE (W m-2)}      
#' 
#' @references Jarvis, P.G., McNaughton, K.G., 1986: Stomatal control of transpiration:
#'             scaling up from leaf to region. Advances in Ecological Research 15, 1-49.
#'             
#'             Monteith, J.L., Unsworth, M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London. 
#'             
#' @seealso \code{\link{decoupling}}            
#'             
#' @examples 
#' df <- data.frame(Tair=20,pressure=100,VPD=seq(0.5,4,0.5),
#'                  Gs_ms=seq(0.01,0.002,length.out=8),Rn=seq(50,400,50))            
#' equilibrium.imposed.ET(df)            
#'             
#' @export
equilibrium.imposed.ET <- function(data,Tair="Tair",pressure="pressure",VPD="VPD",Gs="Gs_ms",
                                   Rn="Rn",G=NULL,S=NULL,missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                                   Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                                   constants=bigleaf.constants()){
  
  check.input(data,list(Tair,pressure,VPD,Rn,Gs,G,S))
  
  if(!is.null(G)){
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    cat("Ground heat flux G is not provided and set to 0.",fill=TRUE)
    G <- 0
  }
  
  if(!is.null(S)){
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    cat("Energy storage fluxes S are not provided and set to 0.",fill=TRUE)
    S <- 0
  }
  
  rho    <- air.density(Tair,pressure,constants)
  gamma  <- psychrometric.constant(Tair,pressure,constants)
  Delta  <- Esat.slope(Tair,Esat.formula,constants)[,"Delta"]
  
  LE_eq  <- (Delta * (Rn - G - S)) / (gamma + Delta)
  LE_imp <- (rho * constants$cp * Gs * VPD) / gamma
  
  ET_imp <- LE.to.ET(LE_imp,Tair)
  ET_eq  <- LE.to.ET(LE_eq,Tair)
  
  return(data.frame(ET_eq,ET_imp,LE_eq,LE_imp))
}
