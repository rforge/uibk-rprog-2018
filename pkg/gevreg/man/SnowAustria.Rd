\name{SnowAustria}
\alias{SnowAustria}
\alias{SnowAustriaGEV}
\alias{SnowAustriaMap}


\title{Random snow depth data for smooth spatial GEV fitting}

\description{
  36 stations with yearly maxima of snow depth for 68 years containing NA's and coordinates longitude (lon), latitude (lat) and elevation (alt), roughly corresponding to the domain of Austia. 
}
\usage{
data("SnowAustria")
}

\format{
Three files named "SnowAustria.rda", "SnowAustriaGEV.rda" and "SnowAustriaMap.rda" are loaded with the command \code{data("SnowAustria", package = "gevreg")}. "SnowAustria.rda" loads the data.frame "SnowAustria" for snow depth data including coordinates at stations. "SnowAustriaGEV.rda" provides the data.frame "SnowAustriaGEV" with GEV parameters valid at the stations, and "SnowAustriaMap.rda" loads the data.frame "SnowAustriaMap".

The first column (\code{station}) in "SnowAustria" is a factor variable containg the levels 1 to 36. Column 2 (\code{hs}) contains snow depth maxima for 68 years per station, including NA's. Columns 3 to 5 (\code{lon}, \code{lat}, \code{alt}) are randomly created coordinates of the 36 stations. For "SnowAustriaGEV" the three columns \code{loc},  \code{scale} and \code{shape} correspond to the location, scale and shape parameters locally estimated at the stations using function \code{\link{gevmle}} of package \code{\link{SpatialExtremes}}. "SnowAustriaMap" provides a gridded set of coordinates of the western part of Austria with Tyrol and Innsbruck in the center, with a horizontal resolution of roughly 1~km. The data.frame comprises three columns "lon", "lat" and "alt" comprising. 
}


\examples{
## Load example data 
data("SnowAustria")
data("SnowAustriaGEV")

## explore dataset
with(SnowAustria, plot(lon,lat))
with(SnowAustria, hist(hs))
df <- cbind(SnowAustriaGEV,data.frame(lon = unique(SnowAustria$lon), alt = unique(SnowAustria$alt)))
plot(df)


}

