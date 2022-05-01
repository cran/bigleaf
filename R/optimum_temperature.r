############################
#### Optimum Temperature ###-------------------------------------------------------------------------------
############################

#' Optimum temperature of Gross Primary Productivity
#' 
#' @description  Calculates the relationship between Gross Primary Productivity (GPP) and Air Temperature (Tair) 
#'               using boundary line analysis and derives the thermal optima.  This function can also be used to find the 
#'               boundary line relationship and optima of other variables such as NPP and NEP.
#'               
#' @param data          Dataframe containing the Gross Primary Productivity and Air Temperature observations  
#' @param GPP           Name of column (in quotations, eg. "GPP") containing the Gross Primary Productivity observations (umol CO2 m-2 s-1).
#' @param Tair          Name of column (in quotations, eg. "Tair") containing the air temperature (degrees celcius) observations.
#' @param BLine         Quantile at which to place the boundary line in format "0.XX".  Defaults to 0.90.
#' @param Obs_filter    Filter to remove air temperature bins with an insufficient number of observations. Defaults to 30.
#'
#' @details  This function works by first binning GPP and air temperature observations to 1 degree temperature bins and then
#'           deriving the relationship between GPP and air temperature at a defined quantile using boundary line analysis. 
#'           Observations are binned using a rounding function, so that each bin is centered on the degree integer value 
#'           (eg. bin 18 contains values between 17.5 and 18.49). The boundary line is usually placed at the upper boundary 
#'           of the distribution (see Webb 1972) however this functional allows the user to select any quantile,
#'           with the default of 0.9 selected for use with eddy covariance flux observations due to the high 
#'           level of noise in these data (see Bennett et al, 2021). After binning observations, the function removes 
#'           temperature bins with fewer observations than the default of 30 (this value can also be user defined).  
#'           It then calculates the smoothed curve between GPP and air temperature using the loess function and derives 
#'           the thermal optima of GPP (Topt). Topt is defined as the temperature bin at which GPP reaches its maximum
#'           along the smoothed boundary line. 
#'  
#' @return A list containing the following objects:
#' \enumerate{
#'          \item{df.bl: A four column dataframe:}
#'          \itemize{
#'                \item{Tair_bin: air temperature bins in 1 degree increments}
#'                \item{GPP_Bline: Value of GPP at the BLine}
#'                \item{n_obs: number of observations in the air temperature bin}
#'                \item{GPP_Bline_smooth: Value of GPP at the smoothed Bline}
#'                }
#'          \item{opt.temp: A named vector with two elements:}
#'          \itemize{
#'                \item{Topt: Thermal optima of GPP - the air temperature bin with maximum GPP along the smoothed Bline}
#'                \item{GPP_bl: The boundary line GPP observation at Topt}
#'                }
#'                }
#'
#' @examples 
#'  # Locate the relationship between GPP and air temperature using default values 
#'  # for BLine and observation filter.
#'
#'  Gpp_ta <- optimum.temperature(data=AT_Neu_Jul_2010, GPP="GPP", Tair="Tair")
#'  
#'  # Locate the relationship between GPP and air temperature at the 50th percentile, 
#'  # filtering temperature bins with fewer than 10 observations
#'  
#'  Gpp_ta <- optimum.temperature(data=AT_Neu_Jul_2010, 
#'                                GPP="GPP", Tair="Tair", BLine=0.50, Obs_filter=10)
#'  
#' @references Bennett A. et al., 2021: Thermal optima of gross primary productivity are closely aligned with mean air temperatures 
#'             across Australian wooded ecosystems. Global Change Biology 32(3), 280-293
#'              
#'             Webb, R. A. 1972. Use of the Boundary Line in the analysis of biological data. Journal of Horticultural Science 47, 309-319
#' 
#' @importFrom stats loess predict
#'              
#' @export
optimum.temperature <- function(data, GPP="GPP", Tair="Tair", BLine=0.9, Obs_filter=30){
 
  check.input(data, list(GPP, Tair))

  #round to 1degC temp bins
  Tair_bin <- trunc(Tair+sign(Tair)*0.5)
  
  #get boundary line
  df.bl <- aggregate(x=GPP,
                     by = list(Tair_bin = Tair_bin),
                     data = data,
                     FUN = quantile,
                     probs = BLine,
                     type = 8)
  

  n_obs <-aggregate(GPP ~ Tair_bin,
                    FUN = length)
  
  df.bl <-merge(df.bl, n_obs, by= c("Tair_bin"))
  
  colnames(df.bl) <- c("Tair_bin", "GPP_Bline", "n_obs")
  
  #Remove Tair bins with n_obs below filter
  df.bl <- subset(df.bl, n_obs >= Obs_filter)
  
  #get the smoothed boundary line
  df.bl$smooth_bl <- predict(loess(GPP_Bline~Tair_bin, df.bl), df.bl$Tair_bin)
  colnames(df.bl) <- c("Tair_bin", "GPP_Bline", "n_obs", "GPP_Bline_smooth")
  
  #get topt
  Topt.df <- df.bl[order(df.bl$GPP_Bline_smooth, decreasing = TRUE), ]
  Topt <- Topt.df[1,1]
  GPP_bl <- Topt.df[1,2]
  opt.temp <- c("Topt" = Topt, "GPP_bl" = GPP_bl)
  
  #output
  optimum.temp <- list(df.bl, opt.temp)
  names(optimum.temp) <- c("df.bl", "opt.temp")
 return(optimum.temp)
 
}
