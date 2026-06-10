\name{checkContextVars}
\alias{checkContextVars}
\title{check for consistency of context variables and items}
\description{Function checks whether some types of context vars, i.e. group,
DIF and weighting variables, are consistent with item variables. The function is
mainly used for internal consistency checks.}
\usage{
checkContextVars (x, varname, type = c("weight", "DIF", "group", "HG"), itemdata,
                  suppressAbort = FALSE, internal = FALSE, renam)}
\arguments{
  \item{x}{
A vector with values of the context variable (e.g., DIF variable)
}
  \item{varname}{
Optional: name of the context variable
}
  \item{type}{
Type of the context variable with following entries allowed: DIF, group, HG, or weight.
}
  \item{itemdata}{
data.frame with item responses
}
  \item{suppressAbort}{
Logical: should the function suppress abort if inconsistencies occur?
}
  \item{internal}{
Logical: is only used for internal use. Recommend to set to FALSE.
}
  \item{renam}{
Optional: A data.frame with two columns containing the renamed item names. This is necessary 
because, when using Conquest, item names cannot exceed ten characters. If there are more than 
ten characters, the items are temporarily renamed internally, and after the Conquest analysis 
is complete, the renaming is reversed. To ensure that the console output contains the original 
item names, a renaming back to the original is also performed for each console output. The 
specification of the \code{renam} argument is only necessary when the function is called
internally; otherwise, the value \code{NULL} should always be used here. 
}
}
\value{
A list
}
\examples{
data(trends)
# first reshape the data for the first time of measurement set into wide format
datW <- reshape2::dcast(trends[which(trends[,"year"] == 2010),],
                        idstud+sex+ses+language~item, value.var="value")
chk1 <- checkContextVars(datW[,"language"], "language", type="DIF",
                         itemdata = datW[,-c(1:4)], renam=NULL)
chk1$info
}

