#########################
### global constants ####
#########################

#' Constants Used in the bigleaf Package
#'
#' @return This function defines the following constants:
#' \item{cp}{specific heat of air for constant pressure (J K-1 kg-1)}
#' \item{Rgas}{universal gas constant (J mol-1 K-1)}
#' \item{Rv}{gas constant of water vapor (J kg-1 K-1) (Stull 1988 p.641)}
#' \item{Rd}{gas constant of dry air (J kg-1 K-1) (Foken p. 245)}
#' \item{Md}{molar mass of dry air [kg mol-1]}
#' \item{Mw}{molar mass of water vapor (kg mol-1)}
#' \item{eps}{ratio of the molecular weight of water vapor to dry air (-) (=Mw/Md)}
#' \item{Kelvin}{conversion degree Celsius to Kelvin}
#' \item{g}{gravitational acceleration (m s-2)}
#' \item{pressure0}{reference atmospheric pressure at sea level (Pa)}
#' \item{Tair0}{reference air temperature (K)}
#' \item{k}{von Karman constant (-)}
#' \item{Cmol}{molar mass of carbon (kg mol-1)}
#' \item{Omol}{molar mass of oxygen (kg mol-1)}
#' \item{sigma}{Stefan-Boltzmann constant (W m-2 K-4)}
#' \item{DwDc}{Ratio of the molecular diffusivities for water vapor and CO2 (-)}
#' \item{Rbwc}{Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-)}
#'
#' @note Constants contained within this function can be changed permanently (which
#'       makes sense for some of them, e.g. for the von Karman constant), by using
#'       the command \code{\link[utils]{fixInNamespace}}. E.g.
#'       
#'       \code{fixInNamespace(bigleaf.constants,ns="bigleaf")}
#'       
#'       Note that this has to be repeated every time the package is newly installed/loaded.
#'       It might thus be easier to change it directly in the source files.
#'
#' @export
bigleaf.constants <- function(){
  
  list(
    cp         = 1004.834,        # specific heat of air for constant pressure (J K-1 kg-1)
    Rgas       = 8.31451,         # universal gas constant (J mol-1 K-1)
    Rv         = 461.5,           # gas constant of water vapor (J kg-1 K-1) (Stull_1988 p.641)
    Rd         = 287.0586,        # gas constant of dry air (J kg-1 K-1) (Foken p. 245)
    Md         = 0.0289645,       # molar mass of dry air (kg mol-1)
    Mw         = 0.0180153,       # molar mass of water vapor (kg mol-1) 
    eps        = 0.622,           # ratio of the molecular weight of water vapor to dry air (-) (=Mw/Md)
    Kelvin     = 273.15,          # conversion degree Celsius to Kelvin
    g          = 9.81,            # gravitational acceleration (m s-2)
    pressure0  = 101325,          # reference atmospheric pressure at sea level (Pa)
    Tair0      = 273.15,          # reference air temperature (K)
    k          = 0.41,            # von Karman constant (-)
    Cmol       = 0.012011,        # molar mass of carbon (kg mol-1)
    Omol       = 0.0159994,       # molar mass of oxygen (kg mol-1)
    sigma      = 5.670367e-08,    # Stefan-Boltzmann constant (W m-2 K-4)
    DwDc       = 1.6,             # Ratio of the molecular diffusivities for water vapor and CO2 (-)
    Rbwc       = 1.37             # Ratio of the transfer efficiency through the boundary layer for water vapor and CO2 (-)
  )
  
}