########################
#### Energy balance ####
########################

#' Biochemical Energy
#' 
#' @description Radiant energy absorbed in photosynthesis or heat release by respiration calculated
#'              from net ecosystem exchange of CO2 (NEE).  
#' 
#' @param NEE     Net ecosystem exchange (umol CO2 m-2 s-1)
#' @param alpha   Energy taken up/released by photosynthesis/respiration per mol CO2 fixed/respired (J umol-1)
#' 
#' @details The following sign convention is employed: NEE is negative when carbon is taken up by the ecosystem.
#'          Positive values of the resulting biochemical energy mean that energy (heat) is taken up by the ecosystem, 
#'          negative ones that heat is released.
#'          The value of alpha is taken from Nobel 1974 (see Meyers & Hollinger 2004), but other values
#'          have been used (e.g. Blanken et al., 1997)
#' 
#' @return \item{Sp -}{biochemical energy (W m-2)}
#' 
#' @references Meyers, T.P., Hollinger, S.E. 2004: An assessment of storage terms in the surface energy
#'             balance of maize and soybean. Agricultural and Forest Meteorology 125, 105-115.
#'             
#'             Nobel, P.S., 1974: Introduction to Biophysical Plant Physiology.
#'             Freeman, New York.
#'             
#'             Blanken, P.D. et al., 1997: Energy balance and canopy conductance of a boreal aspen
#'             forest: Partitioning overstory and understory components. 
#'             Journal of Geophysical Research 102, 28915-28927. 
#'             
#' @examples 
#' # Calculate biochemical energy taken up by the ecosystem with 
#' # a measured NEE of -30umol CO2 m-2 s-1             
#' biochemical.energy(NEE=-30)            
#'            
#' @export 
biochemical.energy <- function(NEE,alpha=0.422){
  Sp <- alpha*-NEE
  return(Sp)
}




#' Energy-Use Efficiency (EUE)
#' 
#' @description Fraction of net radiation fixed by primary productivity.
#' 
#' @param GPP     Gross primary productivity exchange (umol CO2 m-2 s-1)
#' @param alpha   Energy taken up/released by photosynthesis/respiration (J umol-1)
#' @param Rn      Net radiation (W m-2)
#' 
#' @details Energy use efficiency is calculated as:
#' 
#'            \deqn{EUE = sum(GPP)/sum(Rn)}
#'          
#'          where the sums are calculated for complete cases of GPP and Rn over
#'          the entire time period.
#' 
#' @return \item{EUE -}{Energy use efficiency (-)}
#' 
#' @seealso \code{\link{light.use.efficiency}}
#' 
#' @examples
#' energy.use.efficiency(GPP=20,Rn=500)
#' 
#' @export 
energy.use.efficiency <- function(GPP,alpha=0.422,Rn){
  
  Sp <- biochemical.energy(-GPP,alpha)
  
  comp  <- complete.cases(Sp,Rn) 
  
  Sp_sum <- sum(Sp[comp],na.rm=T)
  Rn_sum <- sum(Rn[comp],na.rm=T)
  
  EUE <- Sp_sum/Rn_sum
  
  return(c("EUE"=EUE))
}





#' Energy Balance Closure
#' 
#' @description Calculates the degree of the energy balance non-closure for the entire time span
#'              based on the ratio of two sums (energy balance ratio), and ordinary least squares (OLS).
#' 
#' @param data  Data.frame or matrix containing all required variables.
#' @param Rn    Net radiation (W m-2)
#' @param G     Ground heat flux (W m-2); optional
#' @param S     Sum of all storage fluxes (W m-2); optional
#' @param LE    Latent heat flux (W m-2)
#' @param H     Sensible heat flux (W m-2)
#' @param instantaneous    should the energy balance be calculated at the time step 
#'                         of the observations (\code{TRUE}), or over the entire time period
#'                         provided as input (\code{FALSE})
#' @param missing.G.as.NA  if \code{TRUE}, missing G are treated as \code{NA}s ,otherwise set to 0. 
#' @param missing.S.as.NA  if \code{TRUE}, missing S are treated as \code{NA}s, otherwise set to 0.
#' 
#' 
#' @details The energy balance ratio (EBR) is calculated as:
#'          
#'            \deqn{EBR = sum(LE + H)/sum(Rn - G - S)}
#'          
#'          the sum is taken for all time steps with complete observations (i.e. where
#'          all energy balance terms are available).
#' 
#' @return a named vector containing:
#'         \item{n}{number of complete (all energy balance terms available) observations}
#'         \item{intercept}{intercept of the OLS regression}
#'         \item{slope}{slope of the OLS regression}
#'         \item{r_squared}{r^2 of the OLS regression}
#'         \item{EBR}{energy balance ratio}
#'         
#'         if \code{instantaneous = TRUE}, only \code{EBR} is returned.
#' 
#' @references Wilson K., et al. 2002: Energy balance closure at FLUXNET sites.
#'             Agricultural and Forest Meteorology 113, 223-243.
#'
#' @examples 
#' ## characterize energy balance closure for DE-Tha in June 2014
#' energy.closure(DE_Tha_Jun_2014,instantaneous=FALSE)
#' 
#' ## look at half-hourly closure 
#' EBR_inst <- energy.closure(DE_Tha_Jun_2014,instantaneous=TRUE)
#' summary(EBR_inst)
#' 
#' @importFrom stats complete.cases lm
#' @export
energy.closure <- function(data,Rn="Rn",G=NULL,S=NULL,LE="LE",H="H",instantaneous=FALSE,
                           missing.G.as.NA=FALSE,missing.S.as.NA=FALSE){
  
  check.input(data,list(Rn,LE,H,G,S))
  
  if(!is.null(G)){
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    cat("Ground heat flux G is not provided and set to 0.",fill=TRUE)
    G <- rep(0,nrow(data))
  }
  
  if(!is.null(S)){
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    cat("Energy storage fluxes S are not provided and set to 0.",fill=TRUE)
    S <- rep(0,nrow(data))
  }
  
  if (!instantaneous){
    comp <- complete.cases(Rn,G,S,LE,H)
    n    <- sum(comp)
    
    EBR <- sum(LE[comp] + H[comp]) / sum(Rn[comp] - G[comp] - S[comp])
    
    emod <- lm(c(LE + H) ~ c(Rn - G - S))
    intercept <- summary(emod)$coef[1,1]
    slope     <- summary(emod)$coef[2,1]
    r_squared <- summary(emod)$r.squared
    
    return(c("n"=n,"intercept"=round(intercept,3),"slope"=round(slope,3),"r^2"=round(r_squared,3),"EBR"=round(EBR,3)))
    
  } else {
    
    EBR <- (LE + H) /(Rn - G - S)
    
    return(EBR)
  }
}   



#' Isothermal Net Radiation
#' 
#' @description Calculates the isothermal net radiation, i.e. the net radiation 
#'              that the surface would receive if it had the same temperature than
#'              the air.
#'              
#' @param data       Data.frame or matrix containing all required variables
#' @param Rn         Net radiation (W m-2)
#' @param Tair       Air temperature (degC)
#' @param Tsurf      Surface temperature (degC)
#' @param emissivity Emissivity of the surface (-)
#' @param constants  sigma - Stefan-Boltzmann constant (W m-2 K-4) \cr
#'                   Kelvin - conversion degree Celsius to Kelvin 
#'
#' @details The isothermal net radiation (Rni) is given by:
#'          
#'            \deqn{Rni = Rn + \epsilon * \sigma * (Tsurf^4 - Tair^4)}
#'          
#'          where \eqn{\epsilon} is the emissivity of the surface. Tsurf and Tair
#'          are in Kelvin.
#'          
#' @return \item{Rni -}{isothermal net radiation (W m-2)}
#' 
#' @references Jones, H. 2014: Plants and Microclimate. 3rd edition, Cambridge
#'             University Press.
#' 
#' @examples 
#' # calculate isothermal net radiation of a surface that is 2degC warmer than the air.
#' isothermal.Rn(Rn=400,Tair=25,Tsurf=27,emissivity=0.98) 
#' 
#' @export
isothermal.Rn <- function(data,Rn="Rn",Tair="Tair",Tsurf="Tsurf",emissivity,
                          constants=bigleaf.constants()){
  
  check.input(data,list(Rn,Tair,Tsurf))
  
  Tair  <- Tair + constants$Kelvin
  Tsurf <- Tsurf + constants$Kelvin
  
  Rni <- Rn + emissivity * constants$sigma * (Tsurf^4 - Tair^4)
  
  return(Rni)
  
}