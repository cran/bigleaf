################################################
### Boundary layer conductance formulations #### 
################################################

#' Boundary Layer Conductance according to Thom 1972
#' 
#' @description An empirical formulation for the canopy boundary layer conductance 
#'              for heat transfer based on a simple ustar dependency.
#' 
#' @param ustar     Friction velocity (m s-1)
#' @param Sc        Optional: Schmidt number of additional quantities to be calculated
#' @param Sc_name   Optional: Name of the additional quantities, has to be of same length than 
#'                  \code{Sc_name}
#' @param constants k - von-Karman constant \cr
#'                  Sc_CO2 - Schmidt number for CO2 \cr 
#'                  Pr - Prandtl number (if \code{Sc} is provided)
#'
#'  
#' @details The empirical equation for Rb suggested by Thom 1972 is:
#'  
#'      \deqn{Rb = 6.2ustar^-0.67}
#'  
#'    Gb (=1/Rb) for water vapor and heat are assumed to be equal in this package.
#'    Gb for other quantities x is calculated as (Hicks et al. 1987):
#'  
#'      \deqn{Gb_x = Gb / (Sc_x / Pr)^0.67}
#'  
#'  where Sc_x is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).
#'  
#' @return a data.frame with the following columns:
#'  \item{Gb_h}{Boundary layer conductance for heat transfer (m s-1)}
#'  \item{Rb_h}{Boundary layer resistance for heat transfer (s m-1)}
#'  \item{kB_h}{kB-1 parameter for heat transfer}
#'  \item{Gb_Sc_name}{Boundary layer conductance for \code{Sc_name} (m s-1). Only added if \code{Sc_name} and 
#'                    \code{Sc_name} are provided}
#'  
#' @references Thom, A., 1972: Momentum, mass and heat exchange of vegetation.
#'             Quarterly Journal of the Royal Meteorological Society 98, 124-134.
#'             
#'             Hicks, B.B., Baldocchi, D.D., Meyers, T.P., Hosker, J.R., Matt, D.R., 1987:
#'             A preliminary multiple resistance routine for deriving dry deposition velocities
#'             from measured quantities. Water, Air, and Soil Pollution 36, 311-330.
#' 
#' @seealso \code{\link{Gb.Choudhury}}, \code{\link{Gb.Su}}, \code{\link{aerodynamic.conductance}}
#' 
#' @examples 
#' Gb.Thom(seq(0.1,1.4,0.1))
#' 
#' ## calculate Gb for SO2 as well
#' Gb.Thom(seq(0.1,1.4,0.1),Sc=1.25,Sc_name="SO2")
#' 
#' @export
Gb.Thom <- function(ustar,Sc=NULL,Sc_name=NULL,constants=bigleaf.constants()){
  
  check.input(NULL,ustar)
  
  Rb_h <- 6.2*ustar^-0.667
  Gb_h <- 1/Rb_h
  kB_h <- Rb_h*constants$k*ustar
  
  if (!is.null(Sc) | !is.null(Sc_name)){
    if (length(Sc) != length(Sc_name)){
      stop("arguments 'Sc' and 'Sc_name' must have the same length")
    }
    if (!is.numeric(Sc)){
      stop("argument 'Sc' must be numeric")
    }
  }
  
  Sc   <- c(constants$Sc_CO2,Sc)
  Gb_x <- data.frame(lapply(Sc,function(x) Gb_h / (x/constants$Pr)^0.67))
  colnames(Gb_x) <- paste0("Gb_",c("CO2",Sc_name))
  
  return(data.frame(Gb_h,Rb_h,kB_h,Gb_x))
}




#' Boundary Layer Conductance according to Choudhury & Monteith 1988
#' 
#' @description A formulation for the canopy boundary layer conductance 
#'              for heat transfer according to Choudhury & Monteith 1988.
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
#' @param z0m              Roughness length for momentum (m). If not provided, calculated from \code{roughness.parameters} 
#'                         within \code{wind.profile}
#' @param stab_formulation Stability correction function used (If \code{stab_correction = TRUE}).
#'                         Either \code{"Dyer_1970"} or \code{"Businger_1971"}.
#' @param Sc               Optional: Schmidt number of additional quantities to be calculated
#' @param Sc_name          Optional: Name of the additonal quantities, has to be of same length than 
#'                         \code{Sc_name}
#' @param constants        k - von-Karman constant \cr
#'                         Sc_CO2 - Schmidt number for CO2 \cr 
#'                         Pr - Prandtl number (if \code{Sc} is provided)
#'                         
#' @return A data frame with the following columns:
#'  \item{Gb_h}{Boundary layer conductance for heat transfer (m s-1)}
#'  \item{Rb_h}{Boundary layer resistance for heat transfer (s m-1)}
#'  \item{kB_h}{kB-1 parameter for heat transfer}
#'  \item{Gb_Sc_name}{Boundary layer conductance for \code{Sc_name} (m s-1). Only added if \code{Sc_name} and 
#'                    \code{Sc_name} are provided}
#' 
#' @details Boundary layer conductance according to Choudhury & Monteith 1988 is
#'          given by:
#'          
#'            \deqn{Gb_h = LAI((2a/\alpha)*sqrt(u(h)/w)*(1-exp(-\alpha/2)))}
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
#'          Gb (=1/Rb) for water vapor and heat are assumed to be equal in this package.
#'          Gb for other quantities x is calculated as (Hicks et al. 1987):
#'  
#'            \deqn{Gb_x = Gb / (Sc_x / Pr)^0.67}
#'  
#'          where Sc_x is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).
#'          
#' @note If the roughness length for momentum (\code{z0m}) is not provided as input, it is estimated 
#'       from the function \code{roughness.parameters} within \code{wind.profile}. This function
#'       estimates a single \code{z0m} value for the entire time period! If a varying \code{z0m} value 
#'       (e.g. across seasons or years) is required, \code{z0m} should be provided as input argument.
#'          
#' @references Choudhury, B. J., Monteith J.L., 1988: A four-layer model for the heat
#'             budget of homogeneous land surfaces. Q. J. R. Meteorol. Soc. 114, 373-398.
#'             
#'             McNaughton, K. G., Van den Hurk, B.J.J.M., 1995: A 'Lagrangian' revision of
#'             the resistors in the two-layer model for calculating the energy budget of a
#'             plant canopy. Boundary-Layer Meteorology 74, 261-288.
#'             
#'             Hicks, B.B., Baldocchi, D.D., Meyers, T.P., Hosker, J.R., Matt, D.R., 1987:
#'             A preliminary multiple resistance routine for deriving dry deposition velocities
#'             from measured quantities. Water, Air, and Soil Pollution 36, 311-330.
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
                         leafwidth,LAI,zh,zr,d,z0m=NULL,stab_formulation=c("Dyer_1970","Businger_1971"),
                         Sc=NULL,Sc_name=NULL,constants=bigleaf.constants()){
  
  stab_formulation <- match.arg(stab_formulation)
  
  check.input(data,list(Tair,pressure,wind,ustar,H))
  
  alpha   <- 4.39 - 3.97*exp(-0.258*LAI)

  if (is.null(z0m)){
    estimate_z0m <- TRUE
    z0m <- NULL
  } else {
    estimate_z0m <- FALSE
  }
  
  wind_zh <- wind.profile(data=data,z=zh,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                          zr=zr,estimate_z0m=estimate_z0m,zh=zh,d=d,z0m=z0m,frac_z0m=NULL,
                          stab_correction=TRUE,stab_formulation=stab_formulation)
  
  ## avoid zero windspeed
  wind_zh <- pmax(0.01,wind_zh)
  
  if (!is.null(Sc) | !is.null(Sc_name)){
    if (length(Sc) != length(Sc_name)){
      stop("arguments 'Sc' and 'Sc_name' must have the same length")
    }
    if (!is.numeric(Sc)){
      stop("argument 'Sc' must be numeric")
    }
  }
  
  Gb_h <- LAI*((0.02/alpha)*sqrt(wind_zh/leafwidth)*(1-exp(-alpha/2)))
  Rb_h <- 1/Gb_h
  kB_h <- Rb_h*constants$k*ustar
  
  Sc   <- c(constants$Sc_CO2,Sc)
  Gb_x <- data.frame(lapply(Sc,function(x) Gb_h / (x/constants$Pr)^0.67))
  colnames(Gb_x) <- paste0("Gb_",c("CO2",Sc_name))
  
  
  return(data.frame(Gb_h,Rb_h,kB_h,Gb_x))
}






#' Boundary Layer Conductance according to Su et al. 2001
#' 
#' @description A physically based formulation for the canopy boundary layer conductance
#'              to heat transfer according to Su et al. 2001. 
#'
#' @param data      Data.frame or matrix containing all required variables
#' @param Tair      Air temperature (degC)
#' @param pressure  Atmospheric pressure (kPa)
#' @param ustar     Friction velocity (m s-1)
#' @param wind      Wind speed (m s-1)
#' @param H         Sensible heat flux (W m-2)
#' @param zh        Canopy height (m)
#' @param zr        Reference height (m)
#' @param d         Zero-plane displacement height (-), can be calculated using \code{roughness.parameters}
#' @param z0m       Roughness length for momentum (m). If not provided, calculated from \code{roughness.parameters} 
#'                  within \code{wind.profile}
#' @param Dl        Leaf characteristic dimension (m)
#' @param fc        Fractional vegetation cover [0-1] (if not provided, calculated from LAI)
#' @param LAI       One-sided leaf area index (-)
#' @param N         Number of leaf sides participating in heat exchange (defaults to 2)
#' @param Cd        Foliage drag coefficient (-)
#' @param hs        Roughness height of the soil (m)
#' @param stab_formulation Stability correction function used (If \code{stab_correction = TRUE}).
#'                         Either \code{"Dyer_1970"} or \code{"Businger_1971"}.
#' @param Sc        Optional: Schmidt number of additional quantities to be calculated
#' @param Sc_name   Optional: Name of the additional quantities, has to be of same length than 
#'                  \code{Sc_name}
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Tair0 - reference air temperature (K) \cr
#'                  Sc_CO2 - Schmidt number for CO2 \cr 
#'                  Pr - Prandtl number (if \code{Sc} is provided)
#' 
#' @return A data.frame with the following columns:
#'  \item{Gb_h}{Boundary layer conductance for heat transfer (m s-1)}
#'  \item{Rb_h}{Boundary layer resistance for heat transfer (s m-1)}
#'  \item{kB_h}{kB-1 parameter for heat transfer}
#'  \item{Gb_Sc_name}{Boundary layer conductance for \code{Sc_name} (m s-1). Only added if \code{Sc_name} and 
#'                    \code{Sc_name} are provided}
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
#'            \deqn{Reh = Dl wind(zh) / v}
#'           
#'          kBs-1, the kB-1 value for bare soil surface, is calculated according 
#'          to Su et al. 2001:
#'          
#'            \deqn{kBs^-1 = 2.46(Re)^0.25 - ln(7.4)}
#'          
#'          Gb (=1/Rb) for water vapor and heat are assumed to be equal in this package.
#'          Gb for other quantities x is calculated as (Hicks et al. 1987):
#'  
#'            \deqn{Gb_x = Gb / (Sc_x / Pr)^0.67}
#'  
#'          where Sc_x is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).
#' 
#' @note If the roughness length for momentum (\code{z0m}) is not provided as input, it is estimated 
#'       from the function \code{roughness.parameters} within \code{wind.profile}. This function
#'       estimates a single \code{z0m} value for the entire time period! If a varying \code{z0m} value 
#'       (e.g. across seasons or years) is required, \code{z0m} should be provided as input argument.
#' 
#' 
#' @references Su, Z., Schmugge, T., Kustas, W. & Massman, W., 2001: An evaluation of
#'             two models for estimation of the roughness height for heat transfer between
#'             the land surface and the atmosphere. Journal of Applied Meteorology 40, 1933-1951.
#' 
#'             Massman, W., 1999: A model study of kB H- 1 for vegetated surfaces using
#'            'localized near-field' Lagrangian theory. Journal of Hydrology 223, 27-43.
#'            
#'             Hicks, B.B., Baldocchi, D.D., Meyers, T.P., Hosker, J.R., Matt, D.R., 1987:
#'             A preliminary multiple resistance routine for deriving dry deposition velocities
#'             from measured quantities. Water, Air, and Soil Pollution 36, 311-330.
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
                  H="H",zh,zr,d,z0m=NULL,Dl,fc=NULL,LAI=NULL,N=2,Cd=0.2,hs=0.01,
                  stab_formulation=c("Dyer_1970","Businger_1971"),
                  Sc=NULL,Sc_name=NULL,constants=bigleaf.constants()){
  
  stab_formulation <- match.arg(stab_formulation)
  
  check.input(data,list(Tair,pressure,ustar,wind,H))
  
  if (is.null(fc)){
    if (is.null(LAI)){
      stop("one of 'fc' or 'LAI' must be provided",call.=FALSE)
    } else {
      fc <- (1-exp(-LAI/2)) 
    }
  } 
  
  if (is.null(z0m)){
    estimate_z0m <- TRUE
    z0m <- NULL
  } else {
    estimate_z0m <- FALSE
  }
  
  wind_zh <- wind.profile(data=data,z=zh,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                          zr=zr,estimate_z0m=estimate_z0m,zh=zh,d=d,z0m=z0m,frac_z0m=NULL,
                          stab_correction=TRUE,stab_formulation=stab_formulation)
  
  v   <- kinematic.viscosity(Tair,pressure,constants)
  Re  <- Reynolds.Number(Tair,pressure,ustar,hs,constants)
  kBs <- 2.46 * (Re)^0.25 - log(7.4)
  Reh <- Dl * wind_zh / v
  Ct  <- 1*constants$Pr^-0.6667*Reh^-0.5*N
  
  kB_h <- (constants$k*Cd)/(4*Ct*ustar/wind_zh)*fc^2 + kBs*(1 - fc)^2
  Rb_h <- kB_h/(constants$k*ustar)
  Gb_h <- 1/Rb_h
  
  if (!is.null(Sc) | !is.null(Sc_name)){
    if (length(Sc) != length(Sc_name)){
      stop("arguments 'Sc' and 'Sc_name' must have the same length")
    }
    if (!is.numeric(Sc)){
      stop("argument 'Sc' must be numeric")
    }
  }
  
  Sc   <- c(constants$Sc_CO2,Sc)
  Gb_x <- data.frame(lapply(Sc,function(x) Gb_h / (x/constants$Pr)^0.67))
  colnames(Gb_x) <- paste0("Gb_",c("CO2",Sc_name))
  
  return(data.frame(Gb_h,Rb_h,kB_h,Gb_x))
}




#' Roughness length for heat
#' 
#' @description Roughness length for heat (thermal roughness length, z0h) from the 
#'              kB-1 parameter and roughness length for momentum (z0m).
#'
#' @param z0m   Roughness length for momentum (m)
#' @param kB_h  kB-1 parameter for heat transfer
#' 
#' @return Roughness length for heat, z0h (m)
#' 
#' @details The roughness length for heat (z0h) can be calculated from the 
#'          following relationship (e.g. Verma 1989):
#'       
#'          \deqn{kB_h = ln(z0m/z0h)} 
#'       
#'           it follows:
#'       
#'           \deqn{z0h = z0m / exp(kB_h)}
#' 
#' @note    If unknown, \code{z0m} can be calculated from \code{\link{roughness.parameters}}.
#'          \code{kB_h} can be calculated from \code{\link{Gb.Thom}}, \code{\link{Gb.Choudhury}}, \code{\link{Gb.Su}}
#'          or \code{\link{aerodynamic.conductance}}. 
#'           
#'           
#' @references Verma, S., 1989: Aerodynamic resistances to transfers of heat, mass and momentum.
#'             In: Estimation of areal evapotranspiration, IAHS Pub, 177, 13-20.
#'             
#'             Rigden, A., Li, D., Salvucci, G., 2018: Dependence of thermal roughness length on friction 
#'             velocity across land cover types: A synthesis analysis using AmeriFlux data. Agricultural 
#'             and Forest Meteorology 249, 512-519.
#'           
#' @examples 
#' roughness.length.heat(2,2.5)
#' 
#' @export  
roughness.length.heat <- function(z0m,kB_h){
  
  z0h <- z0m / exp(kB_h)
  
  return(z0h)
}

