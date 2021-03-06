% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_general.R
\name{toRegularTSStack}
\alias{toRegularTSStack}
\title{Convert a raster stack with irregular time observations to regular time series}
\usage{
toRegularTSStack(x, dts, fun, resol)
}
\arguments{
\item{x}{stack of observations, the first raster is a mask (with values 0 and 1). The pixels of the mask that have a value 1 are processed, the pixels with a 0 value are not processed.}

\item{dts}{dates of the observations}

\item{fun}{unction used to aggregate observations to monthly observations. This should be 'max' or 'mean'.}

\item{resol}{desired temporal resolution of the output. This could be 'monthly' or 'daily'}
}
\value{
stack with regular observations
}
\description{
Convert a raster stack with irregular time observations to regular time series
}
