\name{gevregdata}
\alias{gevregdata}


\title{Random snow depth data for smooth spatial GEV fitting}

\description{
  36 stations with yearly maxima of snow depth for 68 years conatinag NA's and coordinates longitude (lon), latitude (lat) and elevation (alt), roughly corresponding to the domain of Austia. 
}
\usage{
data("gevregdata")
}

\format{
A list containg the data frames "hs" for snow depth data including coordinates, and "gev_params" with GEV parameters valid at the stations. The first column (\code{station}) in "hs" is a factor variable containg the levels 1 to 36. Column 2 (\code{hs}) contains snow depth maxima for 68 years per station, including NA's. Columns 3 to 5 (\code{lon}, \code{lat}, \code{alt}) are randomly created coordinates. For "gev_params" the three columns \code{loc},  \code{scale} and \code{shape} correspond to the location, scale and shape parameters locally estimated at the stations using function gevmle of package \code{\link[SpatialExtremes]{SpatialExtremes}}.
}


\examples{
## Load example data
data("gevregdata")
hs <- gevregdata[["hs"]]
gev_params <- gevregdata[["gev_params"]]

## explore dataset
with(hs, plot(lon,lat))
with(hs, hist(hs))
df <- cbind(gev_params,data.frame(lon = unique(hs$lon), alt = unique(hs$alt)))
plot(df)


}

