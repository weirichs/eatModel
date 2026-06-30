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
plotDevianceConquest(logFile, omitUntil=1, reverse=TRUE, change=TRUE, adaptWindow = TRUE)
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
or the deviance should be plotted.
}
  \item{adaptWindow}{
Logical: Should the size of the output window be adjusted to accommodate the
additional information about the R session and available memory? This option
should not be used if the plots are not to be displayed directly but written
to a PDF file.
}
}
\author{
Martin Hecht, Sebastian Weirich
}
