\name{gevregdata}
\alias{gevregdata}


\title{Random snow depth data for GEV fitting}

\description{
  50 Stations with yearly maxima of snow depth for 100 years and coordinates longitude (lon), latitude (lat) and elevation (alt), roughly corresponding to the domain of Austia. 
}
\usage{
data("gevregdata")
}

\format{
A data frame containg 5000 rows and 5 columns. The first column (\code{station}) is a factor variable containg the levels 1 to 3. Column 2 (\code{hs}) contains snow depth maxima for 100 years per station, created with \code{\link[extRemes]{revd}}. Columns 3 to 5 (\code{lon}, \code{lat}, \code{alt}) are randomly created coordinates.  
}



\examples{
## load data
data("gevregdata")

## histogram and summary of snow all depths
hist(gevregdata$hs)
summary(gevregdata$hs)

## plot of locations
plot(gevregdata$lon,gevregdata$lat)



}

