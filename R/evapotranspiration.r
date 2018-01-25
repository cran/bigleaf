###########################
#### Evapotranspiration ###
###########################

#' Potential Evapotranspiration
#' 
#' @description Potential evapotranspiration according to Priestley & Taylor 1972.
#' 
#' @param data      Data.frame or matrix containing all required variables; optional
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
#' @param alpha     Priestley-Taylor coefficient (-)
#' @param missing.G.as.NA  if \code{TRUE}, missing G are treated as \code{NA}s, otherwise set to 0. 
#' @param missing.S.as.NA  if \code{TRUE}, missing S are treated as \code{NA}s, otherwise set to 0. 
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
#' 
#' @details Potential evapotranspiration is calculated according to Priestley & Taylor, 1972:
#' 
#'            \deqn{LE_pot = (\alpha * \Delta * (Rn - G - S)) / (\Delta + \gamma)}
#'
#'          \eqn{\alpha} is the Priestley-Taylor coefficient, \eqn{\Delta} is the slope 
#'          of the saturation vapor pressure curve (kPa K-1),
#'          and \eqn{\gamma} is the psychrometric constant (kPa K-1).
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
#' @references Priestley C.H.B., Taylor R.J., 1972: On the assessment of surface heat flux
#'             and evaporation using large-scale parameters. Monthly Weather Review 100, 81-92.  
#'          
#' @seealso \code{\link{reference.ET}}
#'                                
#' @examples 
#' # Calculate potential ET from a surface that receives Rn of 400 Wm-2
#' potential.ET(Tair=30,pressure=100,Rn=400,alpha=1.26)    
#' 
#' @export
potential.ET <- function(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,alpha=1.26,
                         missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,constants=bigleaf.constants()){
  
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
  Delta  <- Esat.slope(Tair)[,"Delta"]
  
  LE_pot <- (alpha * Delta * (Rn - G - S)) / (Delta + gamma)
  ET_pot <- LE.to.ET(LE_pot,Tair)
  
  return(data.frame(ET_pot,LE_pot))
}




#' Reference Evapotranspiration
#' 
#' @description Reference evapotranspiration calculated from the Penman-Monteith
#'              equation with a prescribed surface conductance.
#' 
#' @param data      Data.frame or matrix containing all required variables
#' @param Gs        Surface conductance (m s-1); defaults to 0.0143 m s-1 (~ 0.58 mol m-2 s-1)
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Ga        Aerodynamic conductance (m s-1)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
#' @param missing.G.as.NA  if \code{TRUE}, missing G are treated as \code{NA}s, otherwise set to 0. 
#' @param missing.S.as.NA  if \code{TRUE}, missing S are treated as \code{NA}s, otherwise set to 0. 
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                  Rgas - universal gas constant (J mol-1 K-1) \cr
#'                  Kelvin - conversion degree Celsius to Kelvin \cr
#'                                
#' @details Reference evapotranspiration is calculated according to the Penman-Monteith
#'          equation:
#' 
#'          \deqn{LE_0 = (\Delta * (Rn - G - S) * \rho * cp * VPD * Ga) / (\Delta + \gamma * (1 + Ga/Gs)}
#'          
#'          where \eqn{\Delta} is the slope of the saturation vapor pressure curve (kPa K-1),
#'          \eqn{\rho} is the air density (kg m-3), and \eqn{\gamma} is the psychrometric constant (kPa K-1).
#'          The reference evapotranspiration is calculated with respect to a 'reference surface',
#'          which is typically a well-watered grass/crop of 0.12m height, an albedo of 0.23 and 
#'          a surface conductance of ~ 0.6 mol m-2 s-1 (Allen et al. 1998), but can be calculated for any other
#'          surface.
#'
#' @return a data.frame with the following columns:
#'         \item{ET_0}{Reference evapotranspiration (kg m-2 s-1)}
#'         \item{LE_0}{Reference latent heat flux (W m-2)}              
#'                  
#' @references  Allen R.G., Pereira L.S., Raes D., Smith M., 1998: Crop evapotranspiration -
#'              Guidelines for computing crop water requirements - FAO Irrigation and drainage
#'              paper 56.
#' 
#' @seealso \code{\link{potential.ET}}
#' 
#' @examples 
#' # Calculate ET_ref for a surface with known Gs (0.5 mol m-2 s-1) and Ga (0.1 m s-1)
#' 
#' # Gs is required in m s-1
#' Gs_ms <- mol.to.ms(0.5,Tair=20,pressure=100)
#' ET_ref <- reference.ET(Gs=Gs_ms,Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400)
#' 
#' # now cross-check with the inverted version
#' surface.conductance(Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,LE=ET_ref[,"LE_ref"])
#' 
#' @export                 
reference.ET <- function(data,Gs=0.0143,Tair="Tair",pressure="pressure",VPD="VPD",Rn="Rn",Ga="Ga",
                         G=NULL,S=NULL,missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                         constants=bigleaf.constants()){
  
  check.input(data,list(Tair,pressure,VPD,Rn,Ga,G,S))
  
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
  Delta  <- Esat.slope(Tair)[,"Delta"]
  rho    <- air.density(Tair,pressure)
  
  LE_ref <- (Delta * (Rn - G - S) + rho * constants$cp * VPD * Ga) / 
            (Delta + gamma * (1 + Ga / Gs))
  
  ET_ref <- LE.to.ET(LE_ref,Tair)
  
  return(data.frame(ET_ref,LE_ref))
  
}





#' Equilibrium and Imposed Evapotranspiration
#' 
#' @description Evapotranspiration (ET) split up into imposed ET and equilibrium ET.
#' 
#' @param data      Data.frame or matrix containing all required input variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param VPD       Air vapor pressure deficit (kPa)
#' @param Gs        surface conductance (m s-1)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
#' @param missing.G.as.NA  if \code{TRUE}, missing G are treated as \code{NA}s, otherwise set to 0. 
#' @param missing.S.as.NA  if \code{TRUE}, missing S are treated as \code{NA}s, otherwise set to 0.
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
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
#' @references Jarvis P.G., McNaughton K.G., 1986: Stomatal control of transpiration:
#'             scaling up from leaf to region. Advances in Ecological Research 15, 1-49.
#'             
#'             Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London. 
#'             
#' @seealso \code{\link{decoupling}}            
#'             
#' @examples 
#' df <- data.frame(Tair=20,pressure=100,VPD=seq(0.5,4,0.5),
#'                  Gs=seq(0.01,0.002,length.out=8),Rn=seq(50,400,50))            
#' equilibrium.imposed.ET(df)            
#'             
#' @export
equilibrium.imposed.ET <- function(data,Tair="Tair",pressure="pressure",VPD="VPD",Gs="Gs",
                                   Rn="Rn",G=NULL,S=NULL,missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                                   constants=bigleaf.constants()){
  
  check.input(data,list(Tair,pressure,VPD,Rn,Gs,G,S))
  
  if(!is.null(G)){
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    cat("ground heat flux G is not provided and set to 0.",fill=TRUE)
    G <- 0
  }
  
  if(!is.null(S)){
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    cat("Energy storage fluxes S are not provided and set to 0.",fill=TRUE)
    S <- 0
  }
  
  rho    <- air.density(Tair,pressure)
  gamma  <- psychrometric.constant(Tair,pressure,constants)
  Delta  <- Esat.slope(Tair)[,"Delta"]
  
  LE_eq  <- (Delta * (Rn - G - S)) / (gamma + Delta)
  LE_imp <- (rho * constants$cp * Gs * VPD) / gamma
  
  ET_imp <- LE.to.ET(LE_imp,Tair)
  ET_eq  <- LE.to.ET(LE_eq,Tair)
  
  return(data.frame(ET_eq,ET_imp,LE_eq,LE_imp))
}
