\name{pvFromRes}
\alias{pvFromRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Collects all plausible values from the object returned by getResults()}
\description{Collects plausible values in a wide or long format data frame. }
\usage{
pvFromRes(resultsObj, toWideFormat = TRUE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultsObj}{
The object returned by \code{\link{getResults}}.
}
  \item{toWideFormat}{
Logical: Should the plausible values be arranged in the wide format? Alternatively, they
will appear in the long format. 
}
}
\details{
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
