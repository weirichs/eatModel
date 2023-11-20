\name{correlationFromRes}
\alias{correlationFromRes}
\title{Gives correlation estimates from the object returned by getResults()}
\description{Gives correlation estimates for multidimensional models. }
\usage{
correlationFromRes(resultsObj, digits = NULL)}
\arguments{
  \item{resultsObj}{
The object returned by \code{\link{getResults}}.
}
  \item{digits}{
Optional: integer value if coefficients should be rounded.
}
}
\value{
A list of data frame(s) with contain(s) correlation matrices.
}
