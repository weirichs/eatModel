\name{plotDevianceMirt}
\alias{plotDevianceMirt}
\title{
Deviance Plot for \pkg{mirt} objects
}
\description{
Plots the deviance change in every iteration for IRT models estimated with \code{\link[mirt]{mirt}}
}
\usage{
plotDevianceMirt(mirt.obj, omitUntil=1)
}
\arguments{
  \item{mirt.obj}{
The object returned by \code{\link[mirt]{mirt}}
}
  \item{omitUntil}{
An optional value indicating number of iterations to be omitted for
plotting.
}}
