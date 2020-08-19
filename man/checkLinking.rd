\name{checkLinking}
\alias{checkLinking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{check for linked design}
\description{Function checks whether all blocks in a complete or incomplete block design are linked to each other.}
\usage{
checkLinking ( design, bookletColumn) }
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
# first and very simple example
des1 <- data.frame ( booklet = paste0("B", 1:6), matrix(data = c(1:2,rep(NA,4), 3, 4, 1,rep(NA,6), 2, 1, 3, rep(NA,4), 4, 2) ,
        nrow = 6, ncol = 4, byrow = TRUE), stringsAsFactors = FALSE)
test1<- checkLinking(design = des1, bookletColumn = "booklet")

# second example: use the original design of the 'IQB Bildungstrend 2016'
# use the 'read_excel' function to import excel design table
library(readxl)
file <- system.file("extdata", "LV16_DE_Testdesign_FINAL.xlsx", package = "eatModel")
des2 <- data.frame(read_excel(file, sheet = 1, col_types = "text", range = "B5:F44"), stringsAsFactors=FALSE)
# design contains three dimensions -- reading, listening and orthography
# linking check must be conducted for each dimension separately. If linking should be checked
# for listening, reading and orthography cells must be set to NA
listening <- des2
for ( i in 2:ncol(listening)) { listening[grep("^D-R|^D-L", listening[,i]),i] <- NA}
test2<- checkLinking(design = listening, bookletColumn = "TH")
# for reading, listening and orthography cells must be set to NA
reading <- des2
for ( i in 2:ncol(reading)) { reading[grep("^D-R|^D-H", reading[,i]),i] <- NA}
test3<- checkLinking(design = reading, bookletColumn = "TH")

# third example: use the original design of the 'IQB Bildungstrend 2018'
file <- system.file("extdata", "Design_2018_final.xlsx", package = "eatModel")
des3 <- data.frame(read_excel(file, sheet = 1, col_types = "text", range = "A6:J86"), stringsAsFactors=FALSE)
cols <- paste0("X", 1:6)
mathe<- des3[which(des3[,3] == "Mathe"),]
test4<- checkLinking(design = mathe[,c("Label", cols)], bookletColumn = "Label")

# sciences
nawi <- des3[which(des3[,3] == "Nawi"),]
# choose 'physik fachwissen'
phyfw<- nawi
for ( i in cols) { phyfw[which(substr(phyfw[,i], 1, 5) != "PhyFw"),i] <- NA}
test5<- checkLinking(design = phyfw[,c("Label", cols)], bookletColumn = "Label")
# choose 'biologie Erkenntnisgewinnung'
bioeg<- nawi
for ( i in cols) { bioeg[which(substr(bioeg[,i], 1, 5) != "BioEg"),i] <- NA}
test6<- checkLinking(design = bioeg[,c("Label", cols)], bookletColumn = "Label")
}

