\name{eapFromRes}
\alias{eapFromRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Collects all EAP results from the object returned by getResults()}
\description{Collects EAP results in a wide data frame. }
\usage{
eapFromRes(resultsObj, idVarName = NULL)}
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data frame in the wide format with ten columns.
}
\author{
Sebastian Weirich
}
\examples{
# see examples in the help file of defineModel()
}
