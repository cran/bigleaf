###############################
### WUE metrics calculation ###
###############################

#' Water-Use Efficiency Metrics
#' 
#' @description Calculation of various water use efficiency (WUE) metrics.
#' 
#' @param data      Data.frame or matrix containing all required variables
#' @param GPP       Gross primary productivity (umol CO2 m-2 s-1)
#' @param NEE       Net ecosystem exchange (umol CO2 m-2 s-1)
#' @param LE        Latent heat flux (W m-2)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Tair      Air temperature (deg C)
#' @param constants Cmol - molar mass of carbon (kg mol-1) \cr
#'                  umol2mol - conversion micromole (umol) to mole (mol) \cr
#'                  kg2g - conversion kilogram (kg) to gram (g)
#'
#' @details the following metrics are calculated:
#' 
#'          Water-use efficiency (WUE):
#'          
#'            \deqn{WUE = GPP / ET}
#'          
#'          Water-use efficiency based on NEE (WUE_NEE):
#'          
#'            \deqn{WUE_NEE = NEE / ET}
#'          
#'          Inherent water-use efficiency (IWUE; Beer et al. 2009):
#'          
#'            \deqn{IWUE = (GPP * VPD) / ET}
#'          
#'          Underlying water-use efficiency (uWUE; Zhou et al. 2014):
#'          
#'            \deqn{uWUE= (GPP * sqrt(VPD)) / ET}
#'          
#'          All metrics are calculated based on the median of all values. E.g.
#'          WUE = median(GPP/ET,na.rm=TRUE)
#' 
#' @return a named vector with the following elements:
#'         \item{WUE}{Water-use efficiency (gC (kg H20)-1)}
#'         \item{WUE_NEE}{Water-use efficiency based on NEE (gC (kg H20)-1)}
#'         \item{IWUE}{Inherent water-use efficiency (gC kPa (kg H20)-1)}
#'         \item{uWUE}{Underlying water-use efficiency (gC kPa^0.5 (kg H20)-1)}
#' 
#' @note Units for VPD can also be hPa. Units change accordingly.
#'       WUE_NEE is calculated based on the absolute value of NEE (the sign convention does not matter here).
#' 
#' @references Beer, C., et al., 2009: Temporal and among-site variability of inherent
#'             water use efficiency at the ecosystem level. Global Biogeochemical Cycles 23, GB2018.
#'             
#'             Zhou, S., et al., 2014: The effect of vapor pressure deficit on water
#'             use efficiency at the sub-daily time scale. Geophysical Research Letters 41.
#'     
#' @seealso \code{\link{stomatal.slope}} for a measure of intrinsic WUE          
#'                             
#' @examples 
#' ## filter data for dry periods and daytime at DE-Tha in June 2014
#' DE_Tha_Jun_2014_2 <- filter.data(DE_Tha_Jun_2014,quality.control=FALSE,
#'                                  vars.qc=c("Tair","precip","VPD","H","LE"),
#'                                  filter.growseas=FALSE,filter.precip=TRUE,
#'                                  filter.vars=c("Tair","PPFD","ustar"),
#'                                  filter.vals.min=c(5,200,0.2),
#'                                  filter.vals.max=c(NA,NA,NA),NA.as.invalid=TRUE,
#'                                  quality.ext="_qc",good.quality=c(0,1),
#'                                  missing.qc.as.bad=TRUE,GPP="GPP",doy="doy",
#'                                  year="year",tGPP=0.5,ws=15,min.int=5,precip="precip",
#'                                  tprecip=0.1,precip.hours=24,records.per.hour=2)
#' 
#' ## calculate WUE metrics in the filtered periods
#' WUE.metrics(DE_Tha_Jun_2014_2)
#'                         
#' @importFrom stats median                                     
#' @export
WUE.metrics <- function(data,GPP="GPP",NEE="NEE",LE="LE",VPD="VPD",Tair="Tair",
                        constants=bigleaf.constants()){
  
  check.input(data,list(GPP,NEE,LE,VPD,Tair))
  
  ET  <- LE.to.ET(LE,Tair)                 # kg H2O m-2 s-1
  GPP <- (GPP * constants$umol2mol * constants$Cmol) * constants$kg2g  # gC m-2 s-1
  NEE <- (NEE * constants$umol2mol * constants$Cmol) * constants$kg2g  # gC m-2 s-1
  
  WUE     <- median(GPP/ET,na.rm=TRUE)
  WUE_NEE <- median(abs(NEE)/ET,na.rm=TRUE)
  IWUE    <- median((GPP*VPD)/ET,na.rm=TRUE)
  uWUE    <- median((GPP*sqrt(VPD))/ET,na.rm=TRUE)
  
  return(c(WUE=WUE,WUE_NEE=WUE_NEE,IWUE=IWUE,uWUE=uWUE))
}
