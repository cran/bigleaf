##############################
#### Surface conductance  ####
##############################

#' Surface Conductance to Water Vapor
#' 
#' @description Calculates surface conductance to water vapor from the inverted Penman-Monteith
#'              equation (by default) or from a simple flux-gradient approach.
#' 
#' @param data      Data.frame or matrix containing all required input variables
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param Rn        Net radiation (W m-2)
#' @param G         Ground heat flux (W m-2); optional
#' @param S         Sum of all storage fluxes (W m-2); optional
#' @param LE        Latent heat flux (W m-2)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Ga        Aerodynamic conductance to heat/water vapor (m s-1)
#' @param missing.G.as.NA  if \code{TRUE}, missing G are treated as \code{NA}s, otherwise they are set to 0.
#'                         Only used if \code{formulation = "Penman-Monteith"}.
#' @param missing.S.as.NA  if \code{TRUE}, missing S are treated as \code{NA}s, otherwise they are set to 0. 
#'                          Only used if \code{formulation = "Penman-Monteith"}.
#' @param formulation Formulation used. Either \code{"Penman-Monteith"} (the default) 
#'                    using the inverted Penman-Monteith equation, or \code{"Flux-Gradient"},
#'                    for a simple flux-gradient approach requiring ET, pressure, and VPD only.
#' @param Esat.formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}. 
#'                      Only used if \code{formulation = "Penman-Monteith"}. See \code{\link{Esat.slope}}.
#' @param constants   cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                    eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                    Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                    Rgas - universal gas constant (J mol-1 K-1) \cr
#'                    Kelvin - conversion degree Celsius to Kelvin \cr
#'                    Mw - molar mass of water vapor (kg mol-1) \cr
#'                    Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' 
#' @details If \code{formulation = "Penman-Monteith"} (the default), surface conductance (Gs) in m s-1 
#'          is calculated from the inverted Penman-Monteith equation:
#' 
#'     \deqn{Gs = ( LE * Ga * \gamma ) / ( \Delta * A + \rho * cp * Ga * VPD - LE * ( \Delta + \gamma ) )}
#'  
#'  Where \eqn{\gamma} is the psychrometric constant (kPa K-1), \eqn{\Delta} is the slope of the 
#'  saturation vapor pressure curve (kPa K-1), and \eqn{\rho} is air density (kg m-3).
#'  Available energy (A) is defined as A = Rn - G - S. If G and/or S are not provided, A = Rn.
#'  
#'  By default, any missing data in G and S are set to 0. If \code{missing.S.as.NA = TRUE}
#'  or \code{missing.S.as.NA = TRUE}, Gs will give \code{NA} for these timesteps.
#'  
#'  If \code{formulation="Flux-Gradient"}, Gs (in mol m-2 s-1) is calculated from VPD and ET only:
#'  
#'     \deqn{Gs = ET/pressure * VPD}
#'  
#'  where ET is in mol m-2 s-1. Note that this formulation assumes fully coupled conditions (i.e. Ga = inf).
#'  This formulation is equivalent to the inverted form of Eq.6 in McNaughton & Black 1973:
#'  
#'     \deqn{Gs = LE * \gamma / (\rho * cp * VPD)}
#'  
#'  which gives Gs in m s-1. Note that Gs > Gc (canopy conductance) under conditions 
#'  when a significant fraction of ET comes from interception or soil evaporation. 
#'  
#'  If \code{pressure} is not available, it can be approximated by elevation using the 
#'  function \code{\link{pressure.from.elevation}}
#'  
#' @return a dataframe with the following columns: 
#'  \item{Gs_ms}{Surface conductance in m s-1}
#'  \item{Gs_mol}{Surface conductance in mol m-2 s-1}
#' 
#' @examples 
#' ## filter data to ensure that Gs is a meaningful proxy to canopy conductance (Gc)
#' DE_Tha_Jun_2014_2 <- filter.data(DE_Tha_Jun_2014,quality.control=FALSE,
#'                                  vars.qc=c("Tair","precip","VPD","H","LE"),
#'                                  filter.growseas=FALSE,filter.precip=TRUE,
#'                                  filter.vars=c("Tair","PPFD","ustar","LE"),
#'                                  filter.vals.min=c(5,200,0.2,0),
#'                                  filter.vals.max=c(NA,NA,NA,NA),NA.as.invalid=TRUE,
#'                                  quality.ext="_qc",good.quality=c(0,1),
#'                                  missing.qc.as.bad=TRUE,GPP="GPP",doy="doy",
#'                                  year="year",tGPP=0.5,ws=15,min.int=5,precip="precip",
#'                                  tprecip=0.1,precip.hours=24,records.per.hour=2)
#' 
#' # calculate Gs based on a simple gradient approach
#' Gs_gradient <- surface.conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                                    VPD="VPD",formulation="Flux-Gradient")
#' summary(Gs_gradient)
#' 
#' # calculate Gs from the the inverted PM equation (now Rn, and Ga are needed),
#' # using a simple estimate of Ga based on Thom 1972
#' Ga <- aerodynamic.conductance(DE_Tha_Jun_2014_2,Rb_model="Thom_1972")[,"Ga_h"]
#' 
#' # if G and/or S are available, don't forget to indicate (they are ignored by default).
#' # Note that Ga is not added to the data.frame 'DE_Tha_Jun_2014'
#' Gs_PM <- surface.conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                              Rn="Rn",G="G",S=NULL,VPD="VPD",Ga=Ga,
#'                              formulation="Penman-Monteith")
#' summary(Gs_PM)
#' 
#'                               
#' # now add Ga to the data.frame 'DE_Tha_Jun_2014' and repeat
#' DE_Tha_Jun_2014_2$Ga <- Ga
#' Gs_PM2 <- surface.conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                               Rn="Rn",G="G",S=NULL,VPD="VPD",Ga="Ga",
#'                               formulation="Penman-Monteith")
#' # note the difference to the previous version (Ga="Ga")
#' summary(Gs_PM2)
#' 
#' @references Monteith, J., 1965: Evaporation and environment. In Fogg, G. E. (Ed.),
#'             The state and movement of water in living organisms (pp.205-234).
#'             19th Symp. Soc. Exp. Biol., Cambridge University Press, Cambridge
#' 
#'             McNaughton, K.G., Black, T.A., 1973: A study of evapotranspiration
#'             from a Douglas Fir forest using the energy balance approach. Water
#'             Resources Research 9, 1579-1590.
#'             
#' @export
surface.conductance <- function(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,
                                VPD="VPD",LE="LE",Ga="Ga_h",missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                                formulation=c("Penman-Monteith","Flux-Gradient"),
                                Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                                constants=bigleaf.constants()){ 
  
  formulation <- match.arg(formulation)
  
  if (formulation == "Flux-Gradient"){
  
    check.input(data,list(Tair,pressure,VPD,LE))
    
    Gs_mol <- (LE.to.ET(LE,Tair)/constants$Mw) * pressure / VPD
    Gs_ms  <- mol.to.ms(Gs_mol,Tair,pressure)
    
  } else if (formulation == "Penman-Monteith"){
    
    check.input(data,list(Tair,pressure,VPD,LE,Rn,Ga,G,S))
    
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
    
    Delta <- Esat.slope(Tair,Esat.formula,constants)[,"Delta"]
    gamma <- psychrometric.constant(Tair,pressure,constants)
    rho   <- air.density(Tair,pressure,constants)
    
    Gs_ms  <- ( LE * Ga * gamma ) / ( Delta * (Rn-G-S) + rho * constants$cp * Ga * VPD - LE * ( Delta + gamma ) )
    Gs_mol <- ms.to.mol(Gs_ms,Tair,pressure)
    
  }
  
  return(data.frame(Gs_ms,Gs_mol))

}