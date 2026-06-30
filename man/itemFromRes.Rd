\name{itemFromRes}
\alias{itemFromRes}
\title{Collect all item results from the object returned by getResults()}
\description{Collect items results in a wide data frame. }
\usage{
itemFromRes(resultsObj)}
\arguments{
  \item{resultsObj}{
The object returned by \code{\link{getResults}}.
}
%  \item{longFormat}{
%Only applies for partial credit models: Should the output apply in wide format (i.e,
%parameters of an item's categories appear side by side in multiple columns), or
%long format (parameters of an item's categories appear below each other in several lines)
%}
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
