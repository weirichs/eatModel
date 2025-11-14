\name{wleRelFromRes}
\alias{wleRelFromRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gives WLE reliability from the object returned by getResults()}
\description{WLE reliability in a data frame. }
\usage{
wleRelFromRes(resultsObj)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultsObj}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}}.
}
}
\value{
A data frame with three columns, containing model name, dimension name, and
the corresponding WLE reliability.
}
\examples{
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
wleRel <- wleRelFromRes(res)
}
