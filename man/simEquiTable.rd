\name{simEquiTable}
\alias{simEquiTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes equivalence table for the Rasch model}
\description{Function provides the equivalence table for unidimensional 1pl models,
specifying the individual competence level for each possible total score of the test.}
\usage{simEquiTable  ( anchor, mRef, sdRef, addConst = 500, multConst = 100, 
cutScores )}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{anchor}{A data frame with anchor parameters on the logit scale, transformed 
  to the metric of the reference population. The first column contains the names of 
  all anchored items. The second column contains anchor parameters.
}
  \item{mRef}{
%%     ~~Describe \code{dif.term} here~~
Scalar: mean of the reference population. 
}
  \item{sdRef}{
%%     ~~Describe \code{split.dif} here~~
Scalar: Standard deviation of the reference population. 
}
  \item{addConst}{
%%     ~~Describe \code{dif.term} here~~
Additive constant for parameter transformation.
}
  \item{multConst}{
%%     ~~Describe \code{abs.dif.bound} here~~
Multiplicative constant for parameter transformation.
}
  \item{cutScores}{
%%     ~~Describe \code{abs.dif.bound} here~~
Named list of one or two elements. "values" is a numeric vector of cut scores (increasing),
"labels" is an optional character vector of cut score labels. Note that "labels" (if specified)
has to be of one more length than "values". 
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of two data frames, including the complete table and the reduced table with
the following 5 columns. 
  \describe{
    \item{score}{Students raw score}
    \item{est}{Estimated individual WLE according to the raw score.}
    \item{bista}{Transformed WLE}
    \item{ks}{competence level}
 }
}
\author{
Johannes Schult
}
\examples{
# read anchor parameter
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)

# use domain 'reading'
prm  <- subset(itemFromRes(res), model == "komplesen")

# use bista cut scores
cuts   <- list ( values = 390+0:3*75, labels = c("I", "II", "III", "IV", "V") )

# create the equivalence table
ret <- simEquiTable( anchor = prm[,c("item", "est")], cutScores = cuts , mRef = 0.039, sdRef = 1.071)
}
