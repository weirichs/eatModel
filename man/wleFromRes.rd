\name{wleFromRes}
\alias{wleFromRes}
\title{Collects all WLEs from the object returned by getResults()}
\description{Collects WLEs in a wide data frame. }
\usage{
wleFromRes(resultsObj, idVarName = NULL, verbose=TRUE)}
\arguments{
  \item{resultsObj}{
The object returned by \code{\link{getResults}}.
}
  \item{idVarName}{
Optional: character vector of individual student id. This is only to provide compatibility
with older package versions. Specification of this argument is only necessary if the function
gives an error.
}
  \item{verbose}{
Print messages to console?
}
}
\value{
A data frame in the wide format with ten columns.
}
\author{
Sebastian Weirich
}
\examples{
# see examples in the help file of defineModel()
}
