\name{replaceLinkingError}
\alias{replaceLinkingError}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Output interface for \code{\link{multiEquatError}} }
\description{\code{\link{multiEquatError} computes linking errors for three
measurement occasions. \code{replaceLinkingError} serves to integrate this
linking error into the object created by \code{\link{equat1pl}}. The linking
error then can be transformed with subsequent function \code{\link{transformToBista}}.
For an illustration of the workflow, see the examples included in the help page
of \code{\link{multiEquatError}}.}
\usage{
replaceLinkingError (equatingList, multiEquatError_output, verbose = TRUE, digits = 4)}
\arguments{
  \item{equatingList}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{equat1pl}}, for measurement occasion 1 vs. 3.
}
  \item{multiEquatError_output}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{multiEquatError}} which contains the linking error
for measurement occasion 1 vs. 3.
}
  \item{verbose}{
Logical: Print information about old (original) and new (replaced) linking error to console?
}
  \item{digits}{
Number of decimals for printing information about old (original) and new (replaced) linking error
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The 'equatingList' object with replaced linking errors. This object can be processed
subsequently by \code{\link{transformToBista}}
}
\examples{
# see the examples in the help page of 'multiEquatError'
}
