\name{runModel}
\alias{runModel}
\title{Run IRT model specified by 'defineModel' using Conquest, TAM, ort mirt}
\description{Start the estimation of an IRT model defined by \code{\link{defineModel}}. \code{runModel}
should be called with the argument returned by \code{\link{defineModel}}.}
\usage{
runModel(defineModelObj, show.output.on.console = FALSE, 
    show.dos.console = TRUE, wait = TRUE, onlySkeleton = FALSE) }
\arguments{
  \item{defineModelObj}{
The object returned by \code{\link{defineModel}}.
}
  \item{show.output.on.console}{
Applies only if \code{\link{defineModel}} previously was called with \code{software = "conquest"}.
Logical: Should the output of the conquest console be printed on the R console during estimation?
}
  \item{show.dos.console}{
Applies only if \code{\link{defineModel}} previously was called with \code{software = "conquest"}.
Logical: Should the output of the conquest console be printed on screen?
}
  \item{wait}{
Applies only if \code{\link{defineModel}} previously was called with \code{software = "conquest"}.
A logical (not NA) indicating whether the R interpreter should wait for the command to finish,
or run it asynchronously.
}
  \item{onlySkeleton}{
Applies only if \code{\link{defineModel}} previously was called with \code{software = "mirt"}.
If \code{TRUE}, the IRT model is not estimated. Instead, only the model specification is returned
as a data.frame, allowing the user to check whether the model has been specified as desired.
}
}
\value{
If \code{\link{defineModel}} previously was called with \code{software = "tam"}, the returned value
is nearly identically to the corresponding TAM output. Accordingly, if \code{\link{defineModel}}
previously was called with \code{software = "mirt"}, the returned value is nearly identically to the
corresponding mirt output. If \code{\link{defineModel}} previously was called with
\code{software = "conquest"}, the returned value contains only internally used information useful for
\code{\link{getResults}}.
}
\author{
Sebastian Weirich
}
\examples{
# see examples in the help file of defineModel()
}
