\name{simEquiTable}
\alias{simEquiTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes equivalence table based on simulated data}
\description{Function provides the equivalence table for unidimensional models, 
specifying the individual competence level for each total score of the test.}
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
\dontrun{
# create arbitrary anchor parameter
anchor <- data.frame ( item = paste("i",1:20,sep=""), 
          par = rnorm(20, mean = -.1, sd = 1.5))

# create arbitrary cut scores
# note that the number of labels (if specified) must equal the number of cuts + 1
cuts   <- list ( values = 330+0:4*75, labels = c("1a", "1b", 2:5) )

# create the equivalence table
ret <- simEquiTable( anchor = anchor, cutScores = cuts , mRef = -0.05, sdRef = 0.9)
View(ret$short)       
}
}
