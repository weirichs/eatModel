\name{plotDevianceConquest}
\alias{plotDevianceConquest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deviance Plot for \pkg{Conquest} log files
}
\description{
Plots the deviance change in every iteration.
}
\usage{
plotDevianceConquest(logFile, omitUntil=1, reverse=TRUE, change=TRUE)
}
\arguments{
  \item{logFile}{
Character string of the Conquest log file
}
  \item{omitUntil}{
An optional value indicating number of iterations to be omitted for
plotting.
}
  \item{reverse}{
A logical indicating whether the deviance change should be
multiplied by minus 1. The default is \code{TRUE}.
}
\item{change}{An optional logical indicating whether deviance change
or the deviance should be plotted.}
}
\author{
Martin Hecht, Sebastian Weirich
}
