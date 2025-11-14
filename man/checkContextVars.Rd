\name{checkContextVars}
\alias{checkContextVars}
\title{check for consistency of context variables and items}
\description{Function checks whether some types of context vars, i.e. group,
DIF and weighting variables, are consistent with item variables. The function is
mainly used for internal consistency checks.}
\usage{
checkContextVars (x, varname, type = c("weight", "DIF", "group", "HG"), itemdata,
                  suppressAbort = FALSE, internal = FALSE)}
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
                         itemdata = datW[,-c(1:4)])
chk1$info
}

