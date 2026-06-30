\name{get.prm}
\alias{get.prm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read Conquest prm (parameter) Output Files}
\description{Read Conquest files comprising item parameters (prm). }
\usage{
get.prm(file)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
Character string with the name of the ConQuest prm file (This file is requested using the \code{"export par"}
statement in the conquest syntax file).
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data frame with three columns:
\item{Case}{Case number}
\item{item}{Identifier for this item }
\item{parameter}{parameter estimate for this item}
}
\examples{
file <- system.file("extdata", "twodim.prm", package = "eatModel")
prm  <- get.prm(file)
}
