\name{itemFromRes}
\alias{itemFromRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Collect all item results from the object returned by getResults()}
\description{Collect items results in a wide data frame. }
\usage{
itemFromRes(resultsObj)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultsObj}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}}.
}
}
\value{
A data frame in the wide format with several columns, containing model name,
item name, dimension, number of valid responses per item, percentage of correct
responses, item discrimination, item difficulty parameter, optional 2pl discrimination
parameters, standard errors, infit and outfit.
}
\examples{
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
items<- itemFromRes(res)
}
