#########################
### global constants ####
#########################

#' Constants Used in the bigleaf Package
#'
#' @description This function defines the following constants:
#' 
#' @param cp           Specific heat of air for constant pressure (J K-1 kg-1)
#' @param Rgas         Universal gas constant (J mol-1 K-1)
#' @param Rv           Gas constant of water vapor (J kg-1 K-1) (Stull 1988 p.641)
#' @param Rd           Gas constant of dry air (J kg-1 K-1) (Foken p. 245)
#' @param Md           Molar mass of dry air (kg mol-1)
#' @param Mw           Molar mass of water vapor (kg mol-1)
#' @param eps          Ratio of the molecular weight of water vapor to dry air (=Mw/Md)
#' @param g            Gravitational acceleration (m s-2)
#' @param pressure0    Reference atmospheric pressure at sea level (Pa)
#' @param Tair0        Reference air temperature (K)
#' @param k            von Karman constant
#' @param Cmol         Molar mass of carbon (kg mol-1)
#' @param Omol         Molar mass of oxygen (kg mol-1)
#' @param sigma        Stefan-Boltzmann constant (W m-2 K-4)
#' @param Pr           Prandtl number
#' @param Sc_CO2       Schmidt number for CO2
#' @param Kelvin       Conversion degree Celsius to Kelvin
#' @param DwDc         Ratio of the molecular diffusivities for water vapor and CO2
#' @param days2seconds Seconds per day
#' @param kPa2Pa       Conversion kilopascal (kPa) to pascal (Pa)
#' @param Pa2kPa       Conversion pascal (Pa) to kilopascal (kPa)
#' @param umol2mol     Conversion micromole (umol) to mole (mol)
#' @param mol2umol     Conversion mole (mol) to micromole (umol) 
#' @param kg2g         Conversion kilogram (kg) to gram (g)
#' @param g2kg         Conversion gram (g) to kilogram (kg)
#' @param kJ2J         Conversion kilojoule (kJ) to joule (J)
#' @param J2kJ         Conversion joule (J) to kilojoule (kJ)
#' @param se_median    Conversion standard error (SE) of the mean to SE of the median
#' @param frac2percent Conversion between fraction and percent
#'
#' @details This function is passed as an argument to every function that uses one 
#'          or more constants. Individual constants passed to a function can be 
#'          easily altered. E.g. the following command will change the value of 
#'          the von Karman constant from 0.41 to 0.4:
#'       
#'          \code{bigleaf.constants(k=0.4)}
#'       
#'          the value of a constant can be returned by calling:
#'       
#'          \code{bigleaf.constants()$*name_of_constant*}
#' 
#'          To permanently change the constants contained within this function (which
#'          makes sense for some of them, e.g. for the von Karman constant), 
#'          the command \code{\link[utils]{fixInNamespace}} can be used. E.g.
#'       
#'          \code{fixInNamespace(bigleaf.constants,ns="bigleaf")}
#'       
#'          Note that this has to be repeated every time the package is newly installed/loaded.
#'       
#' @export
bigleaf.constants <- function(
  ## Physical constants
  cp         = 1004.834,        # specific heat of air for constant pressure (J K-1 kg-1)
  Rgas       = 8.31451,         # universal gas constant (J mol-1 K-1)
  Rv         = 461.5,           # gas constant of water vapor (J kg-1 K-1) (Stull 1988 p.641)
  Rd         = 287.0586,        # gas constant of dry air (J kg-1 K-1) (Foken 2008 p. 245)
  Md         = 0.0289645,       # molar mass of dry air (kg mol-1)
  Mw         = 0.0180153,       # molar mass of water vapor (kg mol-1) 
  eps        = 0.622,           # ratio of the molecular weight of water vapor to dry air (=Mw/Md)
  g          = 9.81,            # gravitational acceleration (m s-2)
  pressure0  = 101325,          # reference atmospheric pressure at sea level (Pa)
  Tair0      = 273.15,          # reference air temperature (K)
  k          = 0.41,            # von Karman constant
  Cmol       = 0.012011,        # molar mass of carbon (kg mol-1)
  Omol       = 0.0159994,       # molar mass of oxygen (kg mol-1)
  sigma      = 5.670367e-08,    # Stefan-Boltzmann constant (W m-2 K-4)
  Pr         = 0.71,            # Prandtl number
  Sc_CO2     = 1.07,            # Schmidt number for CO2 (Hicks et al. 1987)
  
  ## Conversion constants
  Kelvin       = 273.15,         # conversion degree Celsius to Kelvin
  DwDc         = 1.6,            # Ratio of the molecular diffusivities for water vapor and CO2
  days2seconds = 86400,          # seconds per day
  kPa2Pa       = 1000,           # conversion kilopascal (kPa) to pascal (Pa)
  Pa2kPa       = 0.001,          # conversion pascal (Pa) to kilopascal (kPa)
  umol2mol     = 1e-06,          # conversion micromole (umol) to mole (mol)
  mol2umol     = 1e06,           # conversion mole (mol) to micromole (umol) 
  kg2g         = 1000,           # conversion kilogram (kg) to gram (g)
  g2kg         = 0.001,          # conversion gram (g) to kilogram (kg)
  kJ2J         = 1000,           # conversion kilojoule (kJ) to joule (J)
  J2kJ         = 0.001,          # conversion joule (J) to kilojoule (kJ)
  se_median    = 1.253,          # conversion standard error (SE) of the mean to SE of the median (http://influentialpoints.com/Training/standard_error_of_median.htm)
  frac2percent = 100             # conversion between fraction and percent
){
  
  list(
    cp = cp, Rgas = Rgas, Rv = Rv, Rd = Rd, Md = Md, Mw = Mw, eps = eps, g = g,
    pressure0 = pressure0, Tair0 = Tair0, k = k, Cmol = Cmol, Omol = Omol,
    sigma = sigma, Pr = Pr, Sc_CO2 = Sc_CO2, Kelvin = Kelvin, DwDc = DwDc,
    days2seconds = days2seconds, kPa2Pa = kPa2Pa, Pa2kPa = Pa2kPa, umol2mol = umol2mol,
    mol2umol = mol2umol, kg2g = kg2g, g2kg = g2kg, kJ2J = kJ2J, J2kJ = J2kJ,
    se_median = se_median, frac2percent = frac2percent
  )
  
}
