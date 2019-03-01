\name{get.wle}
\alias{get.wle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read Conquest WLE or MLE Output Files}
\description{Read Conquest files comprising maximum likelihood estimates (MLE) or weighted likelihood estimates
(WLE). }
\usage{
get.wle(file)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
Character string with the name of the ConQuest MLE or WLE file.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data frame with one row per person containing the following columns:
\item{case}{Case number}
\item{ID}{Identifier for this case }
\item{n.solved}{Number of items this person answered correctly}
\item{n.total}{Number of total items presented to this person}
\item{wle}{WLE or MLE estimate. The last number of the columns name indicates the
dimension the WLE or MLE estimate belongs to.}
\item{std.wle}{Standard error of WLE or MLE estimate. The last number of the columns name
indicates the dimension the WLE or MLE estimate belongs to.}
}
\references{
See pp. 230 of Wu, M.L., Adams, R.J., Wilson, M.R., & Haldane, S.A. (2007). \emph{ACER ConQuest
Version 2.0. Generalised Item Response Modeling Software.} Camberwell, Victoria: ACER Press.
}
\author{
Sebastian Weirich
}
