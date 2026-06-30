\name{pvFromRes}
\alias{pvFromRes}
\title{Collects all plausible values from the object returned by getResults()}
\description{Collects plausible values in a wide or long format data frame. }
\usage{
pvFromRes(resultsObj, toWideFormat = TRUE, idVarName = NULL, verbose=TRUE)}
\arguments{
  \item{resultsObj}{
The object returned by \code{\link{getResults}}.
}
  \item{toWideFormat}{
Logical: Should the plausible values be arranged in the wide format? Alternatively, they
will appear in the long format. 
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
A data frame in the wide or long format with several columns, containing model name,
model dimension(s), student ID, and several columns with the plausible values.
}
\examples{
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
pv   <- pvFromRes(res)
}
