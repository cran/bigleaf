###########################
#### Surface roughness ####
###########################

#' Roughness Reynolds Number
#' 
#' @description calculates the Roughness Reynolds Number.
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param ustar     Friction velocity (m s-1)
#' @param z0m       Roughness length (m)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Tair0 - reference air temperature (K)
#'                  
#' @details The Roughness Reynolds Number is calculated as in Massman 1999a:
#'          
#'            \deqn{Re = z0m * ustar / v}
#'          
#'          where \code{v} is the kinematic viscosity (m2 s-1).
#'          
#' @return \item{Re -}{Roughness Reynolds Number (-)}
#' 
#' @references Massman, W.J., 1999a: A model study of kB H- 1 for vegetated surfaces using
#'            'localized near-field' Lagrangian theory. Journal of Hydrology 223, 27-43.
#' 
#' @examples 
#' Reynolds.Number(25,100,0.5,z0m=0.5)                             
#' 
#' @export
Reynolds.Number <- function(Tair,pressure,ustar,z0m,constants=bigleaf.constants()){
  
  v  <- kinematic.viscosity(Tair,pressure,constants)
  Re <- z0m*ustar/v
  
  return(Re)
}



#' Roughness Parameters
#' 
#' @description A simple approximation of the two roughness parameters displacement height (d)
#'              and roughness length for momentum (z0m).
#'              
#' @param method    Method to use, one of \code{"canopy_height","canopy_height&LAI","wind_profile"} \cr
#'                  NOTE: if \code{method = "canopy_height"}, only the following three arguments
#'                  are used. If \code{method = "canopy_height&LAI"}, only \code{zh, LAI, cd}, 
#'                  and \code{hs} are required.     
#' @param zh        Vegetation height (m)          
#' @param frac_d    Fraction of displacement height on canopy height (-)
#' @param frac_z0m  Fraction of roughness length on canopy height (-)
#' @param LAI       Leaf area index (-) 
#' @param zr        Instrument (reference) height (m)
#' @param cd        Mean drag coefficient for individual leaves. Defaults to 0.2. 
#'                  Only needed if \code{method = "canopy_height&LAI"}.
#' @param hs        roughness length of the soil surface (m). Only needed if \code{method = "canopy_height&LAI"}
#'                  The following arguments are only needed if \code{method = "wind_profile"}!
#' @param data      Data.frame or matrix containing all required variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param wind      Wind speed at height zr (m s-1)
#' @param ustar     Friction velocity (m s-1)
#' @param H         Sensible heat flux (W m-2)
#' @param d         Zero-plane displacement height (m); optional
#' @param z0m       Roughness length for momentum (m); optional
#' @param stab_roughness   Should stability correction be considered? Default is \code{TRUE}.
#' @param stab_formulation Stability correction function used (If \code{stab_correction = TRUE}).
#'                         Either \code{"Dyer_1970"} or \code{"Businger_1971"}.
#' @param constants k - von-Karman constant (-) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  g - gravitational acceleration (m s-2) \cr
#'                  se_median - conversion standard error (SE) of the mean to SE of the median
#'                  
#' 
#' @details The two main roughness parameters, the displacement height (d)
#'          and the roughness length for momentum (z0m) can be estimated from simple
#'          empirical relationships with canopy height (zh). If \code{method = "canopy_height"},
#'          the following formulas are used:  
#'          
#'            \deqn{d = frac_d * zh}
#'          
#'            \deqn{z0m = frac_z0m * zh}
#'          
#'          where frac_d defaults to 0.7 and frac_z0m to 0.1.
#'          
#'          Alternatively, d and z0m can be estimated from both canopy height and LAI
#'          (If \code{method = "canopy_height&LAI"}).
#'          Based on data from Shaw & Pereira 1982, Choudhury & Monteith 1988 proposed 
#'          the following semi-empirical relations:
#'          
#'            \deqn{X = cd * LAI}
#'          
#'            \deqn{d = 1.1 * zh * ln(1 + X^(1/4))} 
#'          
#'            \deqn{z0m = hs + 0.3 * zh * X^(1/2)   for 0 <= X <= 0.2}
#'          
#'            \deqn{z0m = hs * zh * (1 - d/zh)   for 0.2 < X} 
#'          
#'          If \code{method = "wind_profile"}, z0m is estimated by solving
#'          the wind speed profile for z0m:
#'          
#'            \deqn{z0m = median((zr - d) * exp(-k*wind / ustar - psi_m)}
#'                  
#'          By default, d in this equation is fixed to 0.7*zh, but can be set to any
#'          other value. psi_m is 0 if \code{stab_roughness = FALSE}.       
#' 
#' @return a data.frame with the following columns:
#'         \item{d}{Zero-plane displacement height (m)}
#'         \item{z0m}{Roughness length for momentum (m)}
#'         \item{z0m_se}{Only if \code{method = wind_profile}: Standard Error of the median for z0m (m)}
#'
#' @references Choudhury, B. J., Monteith J.L., 1988: A four-layer model for the heat
#'             budget of homogeneous land surfaces. Q. J. R. Meteorol. Soc. 114, 373-398.
#'             
#'             Shaw, R. H., Pereira, A., 1982: Aerodynamic roughness of a plant canopy: 
#'             a numerical experiment. Agricultural Meteorology, 26, 51-65.
#'    
#' @seealso \code{\link{wind.profile}}
#'     
#' @examples 
#' # estimate d and z0m from canopy height for a dense (LAI=5) and open (LAI=2) canopy
#' roughness.parameters(method="canopy_height&LAI",zh=25,LAI=5)
#' roughness.parameters(method="canopy_height&LAI",zh=25,LAI=2)   
#'    
#' # fix d to 0.7*zh and estimate z0m from the wind profile
#' df <- data.frame(Tair=c(25,25,25),pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=200)
#' roughness.parameters(method="wind_profile",zh=25,zr=40,frac_d=0.7,data=df)
#' 
#' # assume d = 0.8*zh
#' roughness.parameters(method="wind_profile",zh=25,zr=40,frac_d=0.8,data=df) 
#' 
#' @importFrom stats median sd complete.cases 
#' @export                                  
roughness.parameters <- function(method=c("canopy_height","canopy_height&LAI","wind_profile"),zh,
                                 frac_d=0.7,frac_z0m=0.1,LAI,zr,cd=0.2,hs=0.01,data,Tair="Tair",pressure="pressure",
                                 wind="wind",ustar="ustar",H="H",d=NULL,z0m=NULL,
                                 stab_roughness=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                                 constants=bigleaf.constants()){
  
  method           <- match.arg(method)
  stab_formulation <- match.arg(stab_formulation)
  
  if (method == "canopy_height"){
    
    d      <- frac_d*zh
    z0m    <- frac_z0m*zh
    z0m_se <- NA
    
  } else if (method == "canopy_height&LAI"){
    
    X <- cd * LAI
    d <- 1.1 * zh * log(1 + X^(1/4))
    
    if (X >= 0 & X <= 0.2){
      z0m <- hs + 0.3 * X^(1/2)
    } else {
      z0m <- 0.3 * zh * (1 - d/zh)
    }
    z0m_se <- NA
    
  } else if (method == "wind_profile"){
    
    check.input(data,Tair,pressure,wind,ustar,H)
    
    if (is.null(d)){
      
      d <- frac_d * zh
      
    }
    
    if (stab_roughness){
      
      zeta  <- stability.parameter(data=data,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                                   zr=zr,d=d,constants=constants)
      psi_m <- stability.correction(zeta,formulation=stab_formulation)[,"psi_m"]
      z0m_all <- (zr - d) * exp(-constants$k*wind / ustar - psi_m)
      
    } else {
      
      z0m_all <- (zr - d) * exp(-constants$k*wind / ustar)
      
    }
    
    z0m_all[z0m_all > zh] <- NA
    
    z0m    <- median(z0m_all,na.rm=TRUE)
    z0m_se <- constants$se_median * (sd(z0m_all,na.rm=TRUE) / sqrt(length(z0m_all[complete.cases(z0m_all)])))
    
  }
  
  return(data.frame(d,z0m,z0m_se))
}



#' Wind Speed at Given Heights in the Surface Layer
#' 
#' @description Wind speed at a given height above the canopy estimated from single-level
#'              measurements of wind speed.
#'          
#' @param data      Data.frame or matrix containing all required variables
#' @param heights   Vector with heights for which wind speed is to be 
#'                  calculated.
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)                                                                                  
#' @param ustar     Friction velocity (m s-1)
#' @param H         Sensible heat flux (W m-2)
#' @param wind      Wind speed at height zr (m s-1); only used if \code{stab_correction = TRUE}
#' @param zr        Instrument (reference) height (m)
#' @param zh        Canopy height (m)
#' @param d         Zero-plane displacement height (-)
#' @param frac_d    Fraction of displacement height on canopy height (-);
#'                  only used if \code{d} is not available
#' @param z0m       Roughness length (m), optional; only used if \code{stab_correction = FALSE} (default=0.1) 
#' @param frac_z0m  Fraction of roughness length on canopy height (-), optional; only used 
#'                  if \code{stab_correction = FALSE} (default=0.1), only used if \code{z0m} is not available
#' @param estimate_z0m Should \code{z0m} be estimated from the logarithmic wind profile? If \code{TRUE} (the default),
#'                     arguments \code{z0m} and \code{frac_z0m} are ignored.
#'                     See \code{\link{roughness.parameters}} for details. 
#' @param stab_correction Should stability correction be applied? Defaults to \code{TRUE}
#' @param stab_formulation Stability correction function used (If \code{stab_correction = TRUE}).
#'                         Either \code{"Dyer_1970"} or \code{"Businger_1971"}.
#' @param constants k - von-Karman constant (-) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin \cr
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  g - gravitational acceleration (m s-2)
#'
#' @details The underlying assumption is the existence of a logarithmic wind profile
#'          above the height d + z0m (the height at which wind speed mathematically reaches zero
#'          according to the Monin-Obhukov similarity theory).
#'          In this case, the wind speed at a given height z is given by:
#'          
#'            \deqn{u(z) = (ustar/k) * (ln((z - d) / z0m) - \psi{m}}
#'          
#'          The roughness parameters zero-plane displacement height (d) and roughness length (z0m)
#'          can be approximated from \code{\link{roughness.parameters}}. \eqn{\psi{m}} is omitted
#'          if \code{stab_correction = FALSE} (not recommended). If \code{estimate_z0m = TRUE},
#'          z0m is first estimated from the wind profile equation and then used in the equation
#'          above for the calculation of \code{u(z)} (see e.g. Newman & Klein 2014).        
#'                                                             
#' @note Note that this equation is only valid for z >= d + z0m, and it is not 
#'       meaningful to calculate values closely above d + z0m. All values in \code{heights}
#'       smaller than d + z0m will return 0.                                 
#'                                  
#' @return A data.frame with rows representing time and columns representing heights 
#'         as specified in \code{heights}.
#'         
#' @references Monteith, J.L., Unsworth, M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London. 
#'             
#'             Newman, J.F., Klein, P.M., 2014: The impacts of atmospheric stability on
#'             the accuracy of wind speed extrapolation methods. Resources 3, 81-105.
#'         
#' @seealso \code{\link{roughness.parameters}}
#' 
#' @examples 
#' df <- data.frame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250)) 
#' wind.profile(df,heights=seq(18,40,2),zr=40,zh=25,d=16)
#' 
#' @export                                                                                                                          
wind.profile <- function(data,heights,Tair="Tair",pressure="pressure",ustar="ustar",
                         H="H",wind="wind",zr,zh,d=NULL,frac_d=0.7,z0m=NULL,frac_z0m=NULL,
                         estimate_z0m=TRUE,stab_correction=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                         constants=bigleaf.constants()){
  
  stab_formulation <- match.arg(stab_formulation)
  
  check.input(data,ustar)
  
  # data frame containing the results
  wind_heights <- data.frame(matrix(NA,ncol=length(heights),nrow=length(ustar)))
  colnames(wind_heights) <- paste0(heights,"m")
  
  
  ## determine roughness parameters
  if (is.null(d)){
    if (is.null(frac_d)){
      stop("Either 'd' or 'frac_d' must be specified")
    }
    d <- frac_d * zh
  }
  
  if (is.null(z0m) & !estimate_z0m){
    if (is.null(frac_z0m)){
      stop("Either 'z0m' or 'frac_z0m' must be specified if 'estimate_z0m' = FALSE")
    }
    z0m <- frac_z0m * zh
  }
  
  
  if (estimate_z0m){
    
    if (!is.null(z0m) | !is.null(frac_z0m)){
      cat("Note that arguments 'z0m' and 'frac_z0m' are ignored if 'estimate_z0m' = TRUE. z0m is
           calculated from the logarithmic wind_profile equation.",fill=TRUE)
    }
    
    check.input(data,Tair,pressure,wind,ustar,H)
    
    z0m <- roughness.parameters(method="wind_profile",zh=zh,zr=zr,d=d,data=data,
                                Tair=Tair,pressure=pressure,wind=wind,ustar=ustar,H=H,
                                stab_roughness=TRUE,stab_formulation=stab_formulation,
                                constants=constants)[,"z0m"]
    
    #cat("calculated z0m =",round(z0m,2),"m",fill=TRUE)
  }
  

  if ( any(heights < (d + z0m) & !is.na(d + z0m)) ){
    warning("function is only valid for heights above d + z0m! Wind speed for heights below d + z0m will return 0!") 
  }
  
  ## calculate wind speeds at given heights z
  for (z in heights){
    i <- which(heights == z)
    
    if (stab_correction){
      
      zeta  <- stability.parameter(data=data,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                                   zr=z,d=d,constants=constants)
      psi_m <- stability.correction(zeta,formulation=stab_formulation)[,"psi_m"]
      wind_heights[,i] <- pmax(0,(ustar / constants$k) * (log(pmax(0,(z - d)) / z0m) - psi_m))
      
    } else {
      
      wind_heights[,i] <- pmax(0,(ustar / constants$k) * (log(pmax(0,(z - d)) / z0m)))
      
    }
  }
  
  return(wind_heights)
}

