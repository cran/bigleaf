#####################################################
### Stability parameters and stability correction ### ---------------------------------------------------
#####################################################

#' Monin-Obukhov Length
#' 
#' @description calculates the Monin-Obukhov length.
#' 
#' @param data      Data.frame or matrix containing all required variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param ustar     Friction velocity (m s-1)
#' @param H         Sensible heat flux (W m-2)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  k - von Karman constant (-) \cr
#'                  g - gravitational acceleration (m s-2)
#' 
#' @details The Monin-Obukhov length (L) is given by:
#' 
#'              \deqn{L = - (\rho * cp * ustar^3 * Tair) / (k * g * H)}
#'              
#'              where \eqn{rho} is air density (kg m-3).
#' 
#' @return \item{L -}{Monin-Obukhov length (m)}
#' 
#' @note Note that L gets very small for very low ustar values with implications
#'       for subsequent functions using L as input. It is recommended to filter
#'       data and exclude low ustar values (ustar < ~0.2) beforehand. 
#' 
#' @references Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany. 
#' 
#' @seealso \code{\link{stability.parameter}}
#' 
#' @examples 
#' Monin.Obukhov.length(Tair=25,pressure=100,ustar=seq(0.2,1,0.1),H=seq(40,200,20))
#' 
#' @export
Monin.Obukhov.length <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",
                                 H="H",constants=bigleaf.constants()){
  
  check.input(data,list(Tair,pressure,ustar,H))
  
  rho  <- air.density(Tair,pressure,constants)
  Tair <- Tair + constants$Kelvin
  MOL  <- (-rho*constants$cp*ustar^3*Tair) / (constants$k*constants$g*H)
  
  return(MOL)
}



#' Stability Parameter "zeta"
#' 
#' @description calculates "zeta", a parameter characterizing stratification in 
#'              the lower atmosphere.
#' 
#' @param data      Data.frame or matrix containing all required variables
#' @param Tair      Air temperature (degC)
#' @param pressure  Atmospheric pressure (kPa)
#' @param ustar     Friction velocity (m s-1)
#' @param H         Sensible heat flux (W m-2)
#' @param zr        Instrument (reference) height (m)
#' @param d         Zero-plane displacement height (m)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  k - von Karman constant (-) \cr
#'                  g - gravitational acceleration (m s-2)
#' 
#' @details The stability parameter \eqn{\zeta} is given by:
#' 
#'            \deqn{\zeta = (zr - d) / L}
#'          
#'          where L is the Monin-Obukhov length (m), calculated from the function
#'          \code{\link{Monin.Obukhov.length}}. The displacement height d can 
#'          be estimated from the function \code{\link{roughness.parameters}}.
#'          
#' @return \item{\eqn{\zeta} - }{stability parameter (-)}
#' 
#' @examples 
#' df <- data.frame(Tair=25,pressure=100,ustar=seq(0.2,1,0.1),H=seq(40,200,20))
#' stability.parameter(df,zr=40,d=15)
#' 
#' @export           
stability.parameter <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",
                                H="H",zr,d,constants=bigleaf.constants()){
  
  check.input(data,list(Tair,pressure,ustar,H))
  
  MOL  <- Monin.Obukhov.length(data,Tair,pressure,ustar,H,constants)
  zeta <- (zr - d) / MOL
  
  return(zeta)
  
}



#' Integrated Stability Correction Functions for Heat and Momentum
#' 
#' @description dimensionless stability functions needed to correct deviations
#'              from the exponential wind profile under non-neutral conditions.
#'              
#' @param zeta         Stability parameter zeta (-)
#' @param formulation  Formulation for the stability function. Either \code{"Dyer_1970"}, 
#'                     or \code{"Businger_1971"}
#'
#' @details The functions give the integrated form of the universal functions. They
#'          depend on the value of the stability parameter \eqn{\zeta},
#'          which can be calculated from the function \code{\link{stability.parameter}}.
#'          The integration of the universal functions is:
#'          
#'            \deqn{\psi = -x * zeta} 
#'          
#'          for stable atmospheric conditions (\eqn{\zeta} >= 0), and
#'          
#'            \deqn{\psi = 2 * log( (1 + y) / 2) }
#'          
#'          for unstable atmospheric conditions (\eqn{\zeta} < 0).
#'          
#'          The different formulations differ in their value of x and y.
#'   
#' @return a data.frame with the following columns:
#'          \item{psi_h}{the value of the stability function for heat and water vapor (-)}
#'          \item{psi_m}{the value of the stability function for momentum (-)}
#' 
#' @references Dyer, A.J., 1974: A review of flux-profile relationships. 
#'             Boundary-Layer Meteorology 7, 363-372.
#'             
#'             Dyer, A. J., Hicks, B.B., 1970: Flux-Gradient relationships in the
#'             constant flux layer. Quart. J. R. Meteorol. Soc. 96, 715-721.
#'             
#'             Businger, J.A., Wyngaard, J. C., Izumi, I., Bradley, E. F., 1971:
#'             Flux-Profile relationships in the atmospheric surface layer. 
#'             J. Atmospheric Sci. 28, 181-189.
#'             
#'             Paulson, C.A., 1970: The mathematical representation of wind speed
#'             and temperature profiles in the unstable atmospheric surface layer.
#'             Journal of Applied Meteorology 9, 857-861.
#' 
#'             Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
#'
#' @examples 
#' zeta <- seq(-2,0.5,0.05)
#' stability.correction(zeta)
#' stability.correction(zeta,formulation="Businger_1971")                          
#'             
#' @export   
stability.correction <- function(zeta,formulation=c("Dyer_1970","Businger_1971")){
  
  formulation  <- match.arg(formulation)
  
  check.input(NULL,zeta)
  
  psi_h = psi_m <- rep(NA_real_,length(zeta))
  
  # universal functions
  if (formulation == "Businger_1971"){
    x_h <- -7.8
    x_m <- -6
    y_h <- 0.95 * ( 1 - 11.6 * zeta)^0.5
    y_m <- (1 - 19.3*zeta)^0.25
  } else if (formulation == "Dyer_1970"){
    x_h = x_m <- -5
    y_h       <- (1 - 16 * zeta)^0.5
    y_m       <- (1 - 16 * zeta)^0.25
  }
  
  # integration of universal functions (after Paulson_1970 and Foken 2008)
  # stable
  stable <- zeta >= 0 & !is.na(zeta)
  psi_h[stable] <- x_h * zeta[stable]
  psi_m[stable] <- x_m * zeta[stable]
  # unstable
  unstable <- zeta < 0 & !is.na(zeta)
  psi_h[unstable] <- 2 * log( (1 + y_h[unstable] ) / 2)
  psi_m[unstable] <- 2 * log( (1 + y_m[unstable] ) / 2) +
                     log( ( 1 + y_m[unstable]^2 ) / 2)
                     -2 * atan(y_m[unstable]) + pi/2
  
  return(data.frame(psi_h,psi_m))
  
} 
