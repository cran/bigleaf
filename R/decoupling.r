#######################################
#### Canopy-Atmosphere Decoupling  ####
#######################################

#' Canopy-Atmosphere Decoupling Coefficient
#' 
#' @description The canopy-atmosphere decoupling coefficient 'Omega'. 
#' 
#' @param data        Data.frame or matrix containing all required input variables
#' @param Tair        Air temperature (deg C)
#' @param pressure    Atmospheric pressure (kPa)
#' @param Ga          Aerodynamic conductance to heat/water vapor (m s-1)
#' @param Gs          Surface conductance (m s-1)
#' @param approach    Approach used to calculate omega. Either \code{"Jarvis&McNaughton_1986"} (default)
#'                    or \code{"Martin_1989"}.
#' @param LAI         Leaf area index (m2 m-2), only used if \code{approach = "Martin_1989"}.
#' @param Esat.formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}.
#'                      See \code{\link{Esat.slope}}. 
#' @param constants   Kelvin - conversion degree Celsius to Kelvin \cr
#'                    cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                    eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                    sigma - Stefan-Boltzmann constant (W m-2 K-4) \cr
#'                    Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' @details The decoupling coefficient Omega ranges from 0 to 1 and quantifies the
#'          linkage of the conditions (foremost humidity and temperature) at the canopy surface
#'          to the ambient air. Values close to 0 indicate well coupled conditions
#'          characterized by high physiological (i.e. stomatal) control on transpiration
#'          and similar conditions at the canopy surface compared to the atmosphere above
#'          the canopy. Values close to 1 indicate the opposite, i.e. decoupled conditions and 
#'          a low stomatal control on transpiration (Jarvis & McNaughton 1986). \cr
#'          The \code{"Jarvis&McNaughton_1986"} approach (default option) is the original
#'          formulation for the decoupling coefficient, given by (for an amphistomatous 
#'          canopy):
#'          
#'          \deqn{\Omega = \frac{\epsilon + 1}{\epsilon + 1 + \frac{Ga}{Gc}}}{%
#'          \Omega = (\epsilon + 1) / ( \epsilon + 1 + Ga/Gc)}
#'          
#'          where \eqn{\epsilon = \frac{s}{\gamma}}{\epsilon = s/\gamma} is a dimensionless coefficient
#'          with s being the slope of the saturation vapor pressure curve (Pa K-1), and \eqn{\gamma} the 
#'          psychrometric constant (Pa K-1).
#'          
#'          The approach \code{"Martin_1989"} by Martin 1989 additionally takes radiative coupling
#'          into account:
#'          
#'          \deqn{\Omega = \frac{\epsilon + 1 + \frac{Gr}{Ga}}{\epsilon + (1 + \frac{Ga}{Gs}) (1 + \frac{Gr}{Ga})}}{%
#'          \Omega = (\epsilon + 1 + Gr/Ga) / (\epsilon + (1 + Ga/Gs) (1 + Gr/Ga))}
#' 
#' @return \item{\eqn{\Omega} -}{the decoupling coefficient Omega (-)}
#' 
#' @references Jarvis P.G., McNaughton K.G., 1986: Stomatal control of transpiration:
#'             scaling up from leaf to region. Advances in Ecological Research 15, 1-49. 
#'             
#'             Martin P., 1989: The significance of radiative coupling between
#'             vegetation and the atmosphere. Agricultural and Forest Meteorology 49, 45-53.
#' 
#' @seealso \code{\link{aerodynamic.conductance}}, \code{\link{surface.conductance}},
#'          \code{\link{equilibrium.imposed.ET}}
#' 
#' @examples 
#' # Omega calculated following Jarvis & McNaughton 1986
#' set.seed(3)
#' df <- data.frame(Tair=rnorm(20,25,1),pressure=100,Ga_h=rnorm(20,0.06,0.01),
#'                  Gs_ms=rnorm(20,0.005,0.001))
#' decoupling(df,approach="Jarvis&McNaughton_1986")
#' 
#' # Omega calculated following Martin 1989 (requires LAI)
#' decoupling(df,approach="Martin_1989",LAI=4)
#' 
#' @export
decoupling <- function(data,Tair="Tair",pressure="pressure",Ga="Ga_h",Gs="Gs_ms",
                       approach=c("Jarvis&McNaughton_1986","Martin_1989"),
                       LAI,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                       constants=bigleaf.constants()){
  
  approach    <- match.arg(approach)
  
  check.input(data,list(Tair,pressure,Ga,Gs))
  
  Delta   <- Esat.slope(Tair,Esat.formula,constants)[,"Delta"]
  gamma   <- psychrometric.constant(Tair,pressure,constants)
  epsilon <- Delta/gamma
  
  if (approach == "Jarvis&McNaughton_1986"){
    
    Omega <- (epsilon + 1) / (epsilon + 1 + Ga/Gs)
    
  } else if (approach == "Martin_1989") {
    
    if (is.null(LAI)){
      
      stop("LAI is not provided!")
      
    } else {
      
      Gr    <- longwave.conductance(Tair,LAI,constants)
      Omega <- (epsilon + 1 + Gr/Ga) / (epsilon + 1 + Ga/Gs + Gr/Gs + Gr/Ga)
      
    }
  }
  
  return(Omega)
  
}



#' Longwave Radiative Transfer Conductance of the Canopy
#' 
#' @param Tair      Air temperature (deg C)
#' @param LAI       Leaf area index (m2 m-2)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  sigma - Stefan-Boltzmann constant (W m-2 K-4) \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1)
#'                  
#' @details the following formula is used (Martin, 1989):
#' 
#'          \deqn{Gr = 4 \sigma Tair^3 LAI / cp}                             
#'                                       
#' @return \item{Gr -}{longwave radiative transfer conductance of the canopy (m s-1)}
#'                  
#' @references Martin P., 1989: The significance of radiative coupling between
#'             vegetation and the atmosphere. Agricultural and Forest Meteorology 49, 45-53.
#'          
#' @examples 
#' longwave.conductance(25,seq(1,8,1))            
#'                  
#' @export             
longwave.conductance <- function(Tair,LAI,constants=bigleaf.constants()){
  
  Tair <- Tair + constants$Kelvin
  
  Gr <- 4 * constants$sigma * Tair^3 * LAI / constants$cp
  
  return(Gr)
}
