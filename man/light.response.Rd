% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigleaf_physiology.r
\name{light.response}
\alias{light.response}
\title{Ecosystem Light Response}
\usage{
light.response(
  data,
  NEE = "NEE",
  Reco = "Reco",
  PPFD = "PPFD",
  PPFD_ref = 2000,
  ...
)
}
\arguments{
\item{data}{Data.frame or matrix containing all required columns}

\item{NEE}{Net ecosystem exchange (umol CO2 m-2 s-1)}

\item{Reco}{Ecosystem respiration (umol CO2 m-2 s-1)}

\item{PPFD}{Photosynthetic photon flux density (umol m-2 s-1)}

\item{PPFD_ref}{Reference PPFD (umol m-2 s-1) for which GPP_ref is estimated.
Default is 2000 umol m-2 s-1.}

\item{...}{Additional arguments to \code{\link[stats]{nls}}}
}
\value{
A \code{nls} model object containing estimates (+/- SE) for alpha and GPP_ref.
}
\description{
calculates GPP_ref at a reference (usually saturating) PPFD and 
             ecosystem quantum yield (alpha) using a rectangular light response curve.
}
\details{
A rectangular light response curve is fitted to NEE data. The curve
         takes the form as described in Falge et al. 2001:
         
            \deqn{-NEE = \alpha PPFD / (1 - (PPFD / PPFD_ref) + \alpha 
                         PPFD / GPP_ref) - Reco}
                      
         where \eqn{\alpha} is the ecosystem quantum yield (umol CO2 m-2 s-1) (umol quanta m-2 s-1)-1, 
         and GPP_ref is the GPP at the reference PPFD (usually at saturating light). \eqn{\alpha} 
         represents the slope of the light response curve, and is a measure for the light use
         efficiency of the canopy. 
         
         The advantage of this equation over the standard rectangular light response
         curve is that GPP_ref at PPFD_ref is more readily interpretable
         as it constitutes a value observed in the ecosystem, in contrast to 
         GPP_ref (mostly named 'beta') in the standard model that occurs at infinite light.
         \code{PPFD_ref} defaults to 2000 umol m-2 s-1, but other values can be used. For 
         further details refer to Falge et al. 2001.
}
\note{
Note the sign convention. Negative NEE indicates that carbon is taken up
        by the ecosystem. Reco has to be 0 or larger.
}
\references{
Falge E., et al. 2001: Gap filling strategies for defensible annual
            sums of net ecosystem exchange. Agricultural and Forest Meteorology 107,
            43-69.
            
            Gilmanov T.G., et al. 2003: Gross primary production and light response
            parameters of four Southern Plains ecosystems estimated using long-term
            CO2-flux tower measurements. Global Biogeochemical Cycles 17, 1071.
            
            Reichstein M., Stoy P.C., Desai A.R., Lasslop G., Richardson A. 2012: 
            Partitioning of net fluxes. In: Eddy Covariance. A practical guide to
            measurement and data analysis. Aubinet M., Vesala T., Papale D. (Eds.).
            Springer.
}
