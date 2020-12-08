\name{wleFromRes}
\alias{wleFromRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Collects all WLEs from the object returned by getResults()}
\description{Collects WLEs in a wide data frame. }
\usage{
wleFromRes(resultsObj, idVarName = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultsObj}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}}.
}
  \item{idVarName}{
%%     ~~Describe \code{file} here~~
Optional: character vector of individual student id. This is only to provide compatibility
with older package versions. Specification of this argument is only necessary if the function
gives an error.
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
