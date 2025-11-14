\name{plotDevianceMirt}
\alias{plotDevianceMirt}
\title{
Deviance Plot for \pkg{mirt} objects
}
\description{
Plots the deviance change in every iteration for IRT models estimated with \code{\link[mirt]{mirt}}
}
\usage{
plotDevianceMirt(mirt.obj, omitUntil=1, adaptWindow = TRUE)
}
\arguments{
  \item{mirt.obj}{
The object returned by \code{\link[mirt]{mirt}}
}
  \item{omitUntil}{
An optional value indicating number of iterations to be omitted for
plotting.
}
  \item{adaptWindow}{
Logical: Should the size of the output window be adjusted to accommodate the
additional information about the R session and available memory? This option
should not be used if the plots are not to be displayed directly but written
to a PDF file.
}
}
