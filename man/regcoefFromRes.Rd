\name{regcoefFromRes}
\alias{regcoefFromRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gives regression coefficients (betas) from the object returned by getResults()}
\description{Gives coefficients from latent regression model (conditioning model) in a data frame. }
\usage{
regcoefFromRes(resultsObj, digits = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultsObj}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}}.
}
  \item{digits}{
%%     ~~Describe \code{file} here~~
Optional: integer value if coefficients should be rounded.
}
}
\value{
A list of data frame(s) with three columns, containing model name, dimension name, and
the corresponding regression coefficients.
}
\examples{
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
regs <- regcoefFromRes(res)
}
