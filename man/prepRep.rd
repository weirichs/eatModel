\name{prepRep}
\alias{prepRep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prepares 'eatModel' results for jackknifing and trend estimation via 'eatRep'.}
\description{Function uses output of \code{transformToBista} to provide one data.frame
with results from T1 and T2 along with original and transformed linking errors.}
\usage{
prepRep(calibT2, bistaTransfT1, bistaTransfT2, makeIdsUnique = TRUE) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{calibT2}{
%%     ~~Describe \code{file} here~~
The object returned by \code{transformToBista} for the calibration model at T2. 
}
  \item{bistaTransfT1}{
%%     ~~Describe \code{file} here~~
The object returned by \code{transformToBista} for the conditioning model at T1. 
}
  \item{bistaTransfT2}{
%%     ~~Describe \code{file} here~~
The object returned by \code{transformToBista} for the conditioning model at T2. 
}
  \item{makeIdsUnique}{
%%     ~~Describe \code{file} here~~
Logical: Guarantee that person IDs do not overlap between t1 and t2?
}
}
\value{
A single data.frame with person parameter estimates and linking errors for the two
measurement points according to trend estimation. 
}
\author{
Sebastian Weirich
}
\note{
This version is beta. Please use with care!
%%  ~~further notes~~
}
\examples{
# see example 6a in the help file of defineModel() for a detailed demonstration of 
# trend estimation. 
}
