\name{checkDesign}
\alias{checkDesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{check for linked design}
\description{Function checks whether all blocks in a complete or incomplete block design are linked to each other.}
\usage{
checkDesign ( design, bookletColumn) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
%%     ~~Describe \code{file} here~~
A data frame with the test design. All columns (except the \code{booklet} column) are
expected to contain blocks.
}
  \item{bookletColumn}{
%%     ~~Describe \code{file} here~~
Number or name of the booklet column in the design data frame
}
}
\details{
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A logical value (TRUE/FALSE)
}
\author{
Sebastian Weirich
}
\examples{
des1 <- data.frame ( booklet = paste0("B", 1:6), matrix(data = c(1:2,rep(NA,4), 3, 4, 1,rep(NA,6), 2, 1, 3, rep(NA,4), 4, 2) ,
        nrow = 6, ncol = 4, byrow = TRUE), stringsAsFactors = FALSE)
test1<- checkDesign(design = des1, bookletColumn = "booklet")
}

