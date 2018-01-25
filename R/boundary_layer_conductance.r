################################################
### Boundary layer conductance formulations #### 
################################################

#' Boundary Layer Conductance according to Thom 1972
#' 
#' @description An empirical formulation for the canopy boundary layer conductance 
#'              to heat/water vapor based on a simple ustar dependency.
#' 
#' @param ustar     Friction velocity (m s-1)
#' @param constants k - von-Karman constant (-) \cr
#'                  Rbwc - Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-)
#'  
#' @details The empirical equation for Rb to water suggested by Thom 1972 is:
#'  
#'    \deqn{Rb = 6.2ustar^-0.67}
#'  
#'  Rb for water vapor and heat is assumed to be equal. Rb for CO2 (Rb_CO2) is given as:
#'  
#'    \deqn{Rb_CO2 = 1.37 * Rb}
#'  
#'  The factor 1.37 arises due the lower molecular diffusivity of CO2 compared to water.
#'  It is lower than the ratio of the molecular diffusivities (Dw/DCO2 = 1.6), as movement
#'  across the boundary layer is assumed to be partly by diffusion and partly by turbulent
#'  mixing (Nobel 2005).
#'  
#' @return a data.frame with the following columns:
#'  \item{Rb}{Boundary layer resistance for heat and water (s m-1)}
#'  \item{Rb_CO2}{Boundary layer resistance for CO2 (s m-1)}
#'  \item{Gb}{Boundary layer conductance (m s-1)}
#'  \item{kB}{kB-1 parameter (-)}
#' 
#' @references Thom, A., 1972: Momentum, mass and heat exchange of vegetation.
#'             Quarterly Journal of the Royal Meteorological Society 98, 124-134.
#'             
#'             Nobel, P. S., 2005: Physicochemical and Environmental Plant Physiology. Third 
#'             Edition. Elsevier Academic Press, Burlington, USA.
#' 
#' @seealso \code{\link{Gb.Choudhury}}, \code{\link{Gb.Su}}, \code{\link{aerodynamic.conductance}}
#' 
#' @examples 
#' Gb.Thom(seq(0.1,1.4,0.1))
#' 
#' @export
Gb.Thom <- function(ustar,constants=bigleaf.constants()){
  Rb     <- 6.2*ustar^-0.667
  Gb     <- 1/Rb
  kB     <- Rb*constants$k*ustar
  Rb_CO2 <- constants$Rbwc * Rb
  
  return(data.frame(Rb,Rb_CO2,Gb,kB))
}




#' Boundary Layer Conductance according to Choudhury & Monteith 1988
#' 
#' @description A formulation for the canopy boundary layer conductance 
#'              to heat/water vapor according to Choudhury & Monteith 1988.
#'              
#' @param data             Data.frame or matrix containing all required variables
#' @param Tair             Air temperature (degC)
#' @param pressure         Atmospheric pressure (kPa)
#' @param wind             Wind speed at sensor height (m s-1)
#' @param ustar            Friction velocity (m s-1)
#' @param H                Sensible heat flux (W m-2)
#' @param leafwidth        Leaf width (m)
#' @param LAI              One-sided leaf area index
#' @param zh               Canopy height (m)
#' @param zr               Instrument (reference) height (m)
#' @param d                Zero-plane displacement height (-), can be calculated using \code{roughness.parameters}
#' @param stab_formulation Stability correction function used (If \code{stab_correction = TRUE}).
#'                         Either \code{"Dyer_1970"} or \code{"Businger_1971"}.
#' @param constants        k - von-Karman constant (-) \cr
#'                         Rbwc - Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-)
#' 
#' @return A data frame with the following columns:
#'     \item{Rb}{Boundary layer resistance for heat and water (s m-1)}
#'     \item{Rb_CO2}{Boundary layer resistance for CO2 (s m-1)}
#'     \item{Gb}{Boundary layer conductance (m s-1)}
#'     \item{kB}{kB-1 parameter (-)}
#' 
#' @details Boundary layer conductance according to Choudhury & Monteith 1988 is
#'          given by:
#'          
#'            \deqn{Gb = LAI((2a/\alpha)*sqrt(u(h)/w)*(1-exp(-\alpha/2)))}
#'          
#'          where u(zh) is the wind speed at the canopy surface, approximated from
#'          measured wind speed at sensor height zr and a wind extinction coefficient \eqn{\alpha}:
#'          
#'            \deqn{u(zh) = u(zr) / (exp(\alpha(zr/zh -1)))}.
#'          
#'          \eqn{\alpha} is modeled as an empirical relation to LAI (McNaughton & van den Hurk 1995):
#'          
#'            \deqn{\alpha = 4.39 - 3.97*exp(-0.258*LAI)}
#'          
#'          Rb for water vapor and heat is assumed to be equal. Rb for CO2 (Rb_CO2) is given as:
#'  
#'            \deqn{Rb_CO2 = 1.37 * Rb}
#'  
#'          The factor 1.37 arises due the lower molecular diffusivity of CO2 compared to water.
#'          It is lower than the ratio of the molecular diffusivities (Dw/DCO2 = 1.6), as movement
#'          across the boundary layer is assumed to be partly by diffusion and partly by turbulent
#'          mixing (Nobel 2005).
#'          
#' @references Choudhury, B. J., Monteith J.L., 1988: A four-layer model for the heat
#'             budget of homogeneous land surfaces. Q. J. R. Meteorol. Soc. 114, 373-398.
#'             
#'             McNaughton, K. G., Van den Hurk, B.J.J.M., 1995: A 'Lagrangian' revision of
#'             the resistors in the two-layer model for calculating the energy budget of a
#'             plant canopy. Boundary-Layer Meteorology 74, 261-288.
#'             
#'             Nobel, P. S., 2005: Physicochemical and Environmental Plant Physiology. Third 
#'             Edition. Elsevier Academic Press, Burlington, USA.
#'             
#' @seealso \code{\link{Gb.Thom}}, \code{\link{Gb.Su}}, \code{\link{aerodynamic.conductance}}
#'    
#' @examples 
#' ## bulk canopy boundary layer resistance for a closed canopy (LAI=5) 
#' ## with large leaves (leafwdith=0.1)            
#' df <- data.frame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250))    
#' Gb.Choudhury(data=df,leafwidth=0.1,LAI=5,zh=25,d=17.5,zr=40)
#' 
#' ## same conditions, but smaller leaves (leafwidth=0.01)
#' Gb.Choudhury(data=df,leafwidth=0.01,LAI=5,zh=25,d=17.5,zr=40) 
#' 
#' @export                                                                                                                                                                                                                                                                                    
Gb.Choudhury <- function(data,Tair="Tair",pressure="pressure",wind="wind",ustar="ustar",H="H",
                         leafwidth,LAI,zh,zr,d,stab_formulation=c("Dyer_1970","Businger_1971"),
                         constants=bigleaf.constants()){
  
  stab_formulation <- match.arg(stab_formulation)
  
  check.input(data,list(Tair,pressure,wind,ustar,H))
  
  alpha   <- 4.39 - 3.97*exp(-0.258*LAI)
  wind_zh <- wind.profile(data=data,heights=zh,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                          zr=zr,zh=zh,d=d,stab_correction=TRUE,
                          stab_formulation=stab_formulation)[,1]
  
  ## avoid zero windspeed
  wind_zh <- pmax(0.01,wind_zh)
  
  Gb      <- LAI*((0.02/alpha)*sqrt(wind_zh/leafwidth)*(1-exp(-alpha/2)))
  Rb      <- 1/Gb
  kB      <- Rb*constants$k*ustar
  Rb_CO2  <- constants$Rbwc * Rb
  
  return(data.frame(Rb,Rb_CO2,Gb,kB))
}






#' Boundary Layer Conductance according to Su et al. 2001
#' 
#' @description A physically based formulation for the canopy boundary layer conductance
#'              to heat/water vapor according to Su et al. 2001. 
#'
#' @param data     Data.frame or matrix containing all required variables
#' @param Tair     Air temperature (degC)
#' @param pressure Atmospheric pressure (kPa)
#' @param ustar    Friction velocity (m s-1)
#' @param wind     Wind speed (m s-1)
#' @param H        Sensible heat flux (W m-2)
#' @param zh       Canopy height (m)
#' @param zr       Reference height (m)
#' @param d        Zero-plane displacement height (-), can be calculated using \code{roughness.parameters}
#' @param Dl       Leaf characteristic dimension (m)
#' @param fc       Fractional vegetation cover [0-1] (if not provided, calculated from LAI)
#' @param LAI      One-sided leaf area index (-)
#' @param N        Number of leaf sides participating in heat exchange (defaults to 2)
#' @param Cd       Foliage drag coefficient (-)
#' @param hs       Roughness height of the soil (m)
#' @param stab_formulation Stability correction function used (If \code{stab_correction = TRUE}).
#'                         Either \code{"Dyer_1970"} or \code{"Businger_1971"}.
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Tair0 - reference air temperature (K) \cr
#'                  Rbwc - Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-)
#' 
#' @return A data.frame with the following columns:
#'     \item{Rb}{Boundary layer resistance for heat and water (s m-1)}
#'     \item{Rb_CO2}{Boundary layer resistance for CO2 (s m-1)}
#'     \item{Gb}{Boundary layer conductance (m s-1)}
#'     \item{kB}{kB-1 parameter (-)}
#'     
#' @details The formulation is based on the kB-1 model developed by Massman 1999. 
#'          Su et al. 2001 derived the following approximation:
#'           
#'            \deqn{kB-1 = (k Cd fc^2) / (4Ct ustar/u(zh)) + kBs-1(1 - fc)^2}
#'          
#'          If fc (fractional vegetation cover) is missing, it is estimated from LAI:
#' 
#'            \deqn{fc = 1 - exp(-LAI/2)}
#'          
#'          The wind speed at the top of the canopy is calculated using function
#'          \code{\link{wind.profile}}.
#'          
#'          Ct is the heat transfer coefficient of the leaf (Massman 1999):
#'          
#'            \deqn{Ct = Pr^-2/3 Reh^-1/2 N}
#'          
#'          where Pr is the Prandtl number (set to 0.71), and Reh is the Reynolds number for leaves:
#'          
#'            \deqn{Dl wind(zh) / v}
#'           
#'          kBs-1, the kB-1 value for bare soil surface, is calculated according 
#'          to Su et al. 2001:
#'          
#'            \deqn{kBs^-1 = 2.46(Re)^0.25 - ln(7.4)}
#'          
#'          Rb for water vapor and heat is assumed to be equal. Rb for CO2 (Rb_CO2) is given as:
#'  
#'            \deqn{Rb_CO2 = 1.37 * Rb}
#'  
#'          The factor 1.37 arises due the lower molecular diffusivity of CO2 compared to water.
#'          It is lower than the ratio of the molecular diffusivities (Dw/DCO2 = 1.6), as movement
#'          across the boundary layer is assumed to be partly by diffusion and partly by turbulent
#'          mixing (Nobel 2005).
#' 
#' @references Su, Z., Schmugge, T., Kustas, W. & Massman, W., 2001: An evaluation of
#'             two models for estimation of the roughness height for heat transfer between
#'             the land surface and the atmosphere. Journal of Applied Meteorology 40, 1933-1951.
#' 
#'             Massman, W., 1999: A model study of kB H- 1 for vegetated surfaces using
#'            'localized near-field' Lagrangian theory. Journal of Hydrology 223, 27-43.
#'            
#'             Nobel, P. S., 2005: Physicochemical and Environmental Plant Physiology. Third 
#'             Edition. Elsevier Academic Press, Burlington, USA.
#' 
#' @seealso \code{\link{Gb.Thom}}, \code{\link{Gb.Choudhury}}, \code{\link{aerodynamic.conductance}}
#' 
#' @examples 
#' # Canopy boundary layer resistance (and kB-1 parameter) for a set of meteorological conditions,
#' # a leaf characteristic dimension of 1cm, and an LAI of 5
#' df <- data.frame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250)) 
#' Gb.Su(data=df,zh=25,zr=40,d=17.5,Dl=0.01,LAI=5)
#' 
#' # the same meteorological conditions, but larger leaves
#' Gb.Su(data=df,zh=25,zr=40,d=17.5,Dl=0.1,LAI=5)
#' 
#' # same conditions, large leaves, and sparse canopy cover (LAI = 1.5)
#' Gb.Su(data=df,zh=25,zr=40,d=17.5,Dl=0.1,LAI=1.5)
#' 
#' @export
Gb.Su <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",wind="wind",
                  H="H",zh,zr,d,Dl,fc=NULL,LAI=NULL,N=2,Cd=0.2,hs=0.01,
                  stab_formulation=c("Dyer_1970","Businger_1971"),
                  constants=bigleaf.constants()){
  
  stab_formulation <- match.arg(stab_formulation)
  
  check.input(data,list(Tair,pressure,ustar,wind,H))
  
  if (is.null(fc)) {
    if (is.null(LAI)){
      stop("one of 'fc' or 'LAI' must be provided",call.=FALSE)
    } else {
      fc  <- (1-exp(-LAI/2)) 
    }
  } 
  
  wind_zh <- wind.profile(data=data,heights=zh,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                          zr=zr,zh=zh,d=d,stab_correction=TRUE,
                          stab_formulation=stab_formulation)[,1]
  
  v   <- kinematic.viscosity(Tair,pressure,constants)
  Re  <- Reynolds.Number(Tair,pressure,ustar,hs,constants)
  kBs <- 2.46 * (Re)^0.25 - log(7.4)
  Reh <- Dl * wind_zh / v
  Ct  <- 1*0.71^-0.6667*Reh^-0.5*N   # 0.71 = Prandtl number
  
  kB     <- (constants$k*Cd)/(4*Ct*ustar/wind_zh)*fc^2 + kBs*(1 - fc)^2
  Rb     <- kB/(constants$k*ustar)
  Rb_CO2 <- constants$Rbwc * Rb
  Gb     <- 1/Rb
  
  
  return(data.frame(Rb,Rb_CO2,Gb,kB))
}
