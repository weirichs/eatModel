\name{q3FromRes}
\alias{q3FromRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extracts the Q3 matrices from the object returned by getResults()}
\description{Collects Q3 matrices in a list of data frames. }
\usage{
q3FromRes(resultsObj, out = c("wide", "long" ), triangular = FALSE)}
\arguments{
  \item{resultsObj}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}}.
}
  \item{out}{
Specifies the output format. \code{"wide"} gives a triangular matrix,
\code{"long"} gives a long format table with the residual correlation in
a separate column.
}
  \item{triangular}{
Logical: should the wide-format matrix be arranged in triangular shape?
}
}
\value{
A list of data frames, one data.frame per model.
}
\examples{
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
q3   <- q3FromRes(res)
}
