\name{eapFromRes}
\alias{eapFromRes}
\title{Collects all EAP results from the object returned by getResults()}
\description{Collects EAP results in a wide data frame. }
\usage{
eapFromRes(resultsObj, idVarName = NULL, verbose = TRUE)}
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
A data frame with several columns, model name, student ID, dimension name,
group name, EAP and the corresponding standard error.
}
\examples{
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
eap  <- eapFromRes(res)
}
