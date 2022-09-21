\name{checkLinking}
\alias{checkLinking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{check for linked design}
\description{Function checks whether all blocks in a complete or incomplete block
design are linked to each other.}
\usage{
checkLinking(design, blocks=NULL, bookletColumn=NULL, verbose=FALSE) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
%%     ~~Describe \code{file} here~~
A data frame with the test design. All columns (except for the \code{booklet} identifier column) are
expected to contain blocks.
}
  \item{blocks}{
%%     ~~Describe \code{file} here~~
Optional. To check whether a subdomain is completely linked add here all the blocks that belong to this subdomain in a character vector.
}
  \item{bookletColumn}{
%%     ~~Describe \code{file} here~~
Optional. Number or name of the booklet identifier column in the design data frame.
}
  \item{verbose}{
%%     ~~Describe \code{file} here~~
Optional. If \code{TRUE} the function gives more messages.
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
Sebastian Weirich, Karoline Sachse
}
\examples{
# 1. first examples
# a) design linked
des1a <- data.frame(booklet = paste0("B", 1:4),
       Pos1 = c("blockA", "blockB", "blockC", "blockD"),
       Pos2 = c("blockB", "blockC", "blockD", "blockE"),
       Pos3 = c("blockC", "blockD", "blockE", "blockF"))
test1 <- checkLinking(design = des1a, bookletColumn = "booklet")

# b) design not linked:
des1b <- data.frame(Pos1 = c("blockA", "blockH", "blockC", "blockF"),
       Pos2 = c("blockB", "blockC", "blockD", "blockA"),
       Pos3 = c("blockF", "blockG", "blockH", "blockB"))
test2 <- checkLinking(design = des1b)

# 2. second examples: use a 'IQB Bildungstrend 2016'-like design
data(des2)

# design contains three dimensions -- reading, listening and orthography
# linking check must be conducted for each dimension separately.

# a) domain reading (blocks contain "-L")
readblocks <- grep("-L", unique(unlist(des2[,-1])), value=TRUE)
test3 <- checkLinking(design = des2, blocks =readblocks, bookletColumn = "TH")

# b) domain listening (blocks contain "-H")
listenblocks <- grep("-H", unique(unlist(des2[,-1])), value=TRUE)
test4 <- checkLinking(design = des2[,-1], blocks =listenblocks)

# c) domain orthography (blocks contain "-R")
orthoblocks <- grep("-R", unique(unlist(des2[,-1])), value=TRUE)
test5 <- checkLinking(design = des2[,-1], blocks =orthoblocks)
}