\name{eapRelFromRes}
\alias{eapRelFromRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gives EAP reliability from the object returned by getResults()}
\description{EAP reliability in a data frame. }
\usage{
eapRelFromRes(resultsObj)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultsObj}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}}.
}
}
\value{
A data frame with three columns, model name, domain name, and the
EAP reliability.
}
\examples{
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
eapRel <- eapRelFromRes(res)
}
