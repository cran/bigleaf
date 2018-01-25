#######################
## Unit Conversions ### 
#######################

#' Conversion between Latent Heat Flux and Evapotranspiration
#' 
#' @description converts evaporative water flux from mass (ET=evapotranspiration)
#'              to energy (LE=latent heat flux) units, or vice versa.
#'              
#' @aliases LE.to.ET ET.to.LE
#' 
#' @param LE   Latent heat flux (W m-2)
#' @param ET   Evapotranspiration (kg m-2 s-1)
#' @param Tair Air temperature (deg C)
#' 
#' @details 
#' The conversions are given by:
#' 
#' \deqn{ET = LE/\lambda}
#' 
#' \deqn{LE = \lambda ET}
#' 
#' where \eqn{\lambda} is the latent heat of vaporization (J kg-1) as calculated by
#' \code{\link{latent.heat.vaporization}}.
#' 
#' @examples 
#' # LE of 200 Wm-2 and air temperature of 25degC
#' LE.to.ET(200,25)
#' 
#' @export
LE.to.ET <- function(LE,Tair){
  
  lambda <- latent.heat.vaporization(Tair)
  ET     <- LE/lambda
  
  return(ET)
}


#' @rdname LE.to.ET
#' @export
ET.to.LE <- function(ET,Tair){
  
  lambda <- latent.heat.vaporization(Tair)
  LE     <- ET*lambda
  
  return(LE)
}



#' Conversion between Conductance Units
#' 
#' @description Converts conductances from mass (m s-1)
#'              to molar units (mol m-2 s-1), or vice versa.
#' 
#' @aliases ms.to.mol mol.to.ms
#' 
#' @param G_ms       Conductance (m s-1)
#' @param G_mol      Conductance (mol m-2 s-1)
#' @param Tair       Air temperature (deg C)
#' @param pressure   Atmospheric pressure (kPa)
#' @param constants  Kelvin - conversion degree Celsius to Kelvin \cr
#'                   Rgas - universal gas constant (J mol-1 K-1)
#' 
#' @details 
#' The conversions are given by:
#' 
#' \deqn{G_mol = G_ms * pressure / (Rgas * Tair)}
#' 
#' \deqn{G_ms = G_mol * (Rgas * Tair) / pressure}
#' 
#' where Tair is in Kelvin and pressure in Pa (converted from kPa internally)
#' 
#' @references Jones, H.G. 1992. Plants and microclimate: a quantitative approach to environmental plant physiology.
#'             2nd Edition., Cambridge University Press, Cambridge. 428 p
#'             
#' @examples 
#' ms.to.mol(0.005,25,100)
#'             
#' @export
ms.to.mol <- function(G_ms,Tair,pressure,constants=bigleaf.constants()){
  
  Tair     <- Tair + constants$Kelvin
  pressure <- pressure * 1000
  
  G_mol  <- G_ms * pressure / (constants$Rgas * Tair)
  
  return(G_mol)
}


#' @rdname ms.to.mol
#' @export
mol.to.ms <- function(G_mol,Tair,pressure,constants=bigleaf.constants()){
  
  Tair     <- Tair + constants$Kelvin
  pressure <- pressure * 1000
  
  G_ms  <- G_mol * (constants$Rgas * Tair) / (pressure)
  
  return(G_ms)
}



#' Conversions between Humidity Measures
#' 
#' @description Conversion between vapor pressure (e), vapor pressure deficit (VPD), 
#'              specific humidity (q), and relative humidity (rH).
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param e         Vapor pressure (kPa)
#' @param q         Specific humidity (kg kg-1)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param rH        Relative humidity (-)
#' @param constants eps - ratio of the molecular weight of water vapor to dry air (-)
#' 
#' @family humidity conversion
#' 
#' @references Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
#' 
#' @export
VPD.to.rH <- function(VPD,Tair){
  esat <- Esat.slope(Tair)[,"Esat"]
  rH   <- 1 - VPD/esat
  return(rH)
} 


#' @rdname VPD.to.rH
#' @family humidity conversion
#' @export
rH.to.VPD <- function(rH,Tair){
  if(rH > 1){
    warning("relative humidity (rH) has to be between 0 and 1.")
  }
  esat <- Esat.slope(Tair)[,"Esat"]
  VPD  <- esat - rH*esat
  return(VPD)
} 


#' @rdname VPD.to.rH
#' @family humidity conversion
#' @export
VPD.to.e <- function(VPD,Tair){
  esat <- Esat.slope(Tair)[,"Esat"]
  e    <- esat - VPD
  return(e)
}


#' @rdname VPD.to.rH
#' @family humidity conversion
#' @export
e.to.VPD <- function(e,Tair){
  esat <- Esat.slope(Tair)[,"Esat"]
  VPD  <- esat - e 
  return(VPD)
}


#' @rdname VPD.to.rH
#' @family humidity conversion
#' @export
e.to.q <- function(e,pressure,constants=bigleaf.constants()){
  q <- constants$eps * e / (pressure - (1-constants$eps) * e) 
  return(q)
}


#' @rdname VPD.to.rH
#' @family humidity conversion
#' @export
q.to.e <- function(q,pressure,constants=bigleaf.constants()){
  e <- q * pressure / ((1-constants$eps) * q + constants$eps)
  return(e)
}


#' @rdname VPD.to.rH
#' @family humidity conversion
#' @export
q.to.VPD <- function(q,Tair,pressure,constants=bigleaf.constants()){
  esat <- Esat.slope(Tair)[,"Esat"]
  e    <- q.to.e(q,pressure,constants)
  VPD  <- esat - e
  return(VPD)
} 


#' @rdname VPD.to.rH
#' @family humidity conversion
#' @export
VPD.to.q <- function(VPD,Tair,pressure,constants=bigleaf.constants()){
  esat <- Esat.slope(Tair)[,"Esat"]
  e    <- esat - VPD
  q    <- e.to.q(e,pressure,constants)
  return(q) 
} 





#' Conversions between Global Radiation and Photosynthetic Photon Flux Density
#' 
#' @description Converts radiation from W m-2 to umol m-2 s-1 and vice versa.
#' 
#' @param Rg       Global radiation = incoming short-wave radiation at the surface (W m-2)
#' @param PPFD     Photosynthetic photon flux density (umol m-2 s-1)
#' @param J_to_mol Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
#' @param frac_PAR Fraction of incoming solar irradiance that is photosynthetically 
#'                 active radiation (PAR); defaults to 0.5
#'                 
#' @details 
#' The conversion is given by:
#' 
#'  \deqn{PPFD = Rg * frac_PAR * J_to_mol}
#'  
#' by default, the combined conversion factor (\code{frac_PAR * J_to_mol}) is 2.3
#'
#' @examples 
#' # convert a measured incoming short-wave radiation of 500 Wm-2 to 
#' # PPFD in umol m-2 s-1 and backwards
#' Rg.to.PPFD(500)
#' PPFD.to.Rg(1150)
#' 
#' @export
Rg.to.PPFD <- function(Rg,J_to_mol=4.6,frac_PAR=0.5){
  PPFD <- Rg * frac_PAR * J_to_mol
  return(PPFD)
}


#' @rdname Rg.to.PPFD
#' @export
PPFD.to.Rg <- function(PPFD,J_to_mol=4.6,frac_PAR=0.5){
  Rg <- PPFD / frac_PAR / J_to_mol
  return(Rg)
} 




#' Conversion between Mass and Molar Units of Carbon and CO2
#' 
#' @description Converts CO2 quantities from umol CO2 m-2 s-1 to g C m-2 d-1 and vice versa.
#' 
#' @param CO2_flux  CO2 flux (umol CO2 m-2 s-1)
#' @param C_flux    Carbon (C) flux (gC m-2 d-1)
#' @param constants Cmol - molar mass of carbon (kg mol-1)
#' 
#' @examples 
#' umolCO2.to.gC(20)  # gC m-2 d-1
#' 
#' @export
umolCO2.to.gC <- function(CO2_flux,constants=bigleaf.constants()){
  
  C_flux <- CO2_flux * 1e-06 * constants$Cmol * 1000 * 86400
  
  return(C_flux)
}


#' @rdname umolCO2.to.gC
#' @export
gC.to.umolCO2 <- function(C_flux,constants=bigleaf.constants()){
  
  CO2_flux <- (C_flux / 1000 / 86400) / constants$Cmol * 1e06
  
  return(CO2_flux)
}