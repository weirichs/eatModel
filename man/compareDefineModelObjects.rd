\name{compareDefineModelObjects}
\alias{compareDefineModelObjects}
\title{Compares two \code{\link{defineModel}} objects}
\description{Function facilitates quality checks by comparing two objects returned
by \code{\link{defineModel}} whether there are functionally equivalent. }
\usage{
compareDefineModelObjects (ref, tar, round = TRUE, digits = 3)}
\arguments{
  \item{ref}{
The reference object returned by \code{\link{defineModel}}.
}
  \item{tar}{
The target object returned by \code{\link{defineModel}}.
}
  \item{round}{
Logical. Round decimal values before the comparison?
}
  \item{digits}{
Only necessary if \code{round = TRUE}. Number of digits for rounding
}
}
\value{
No value is returned. Function prints messages to the console if inequalities occur.
}
\examples{
data(trends)

# first reshape the data for the first time of measurement set into wide format
datW <- reshape2::dcast(trends[which(trends[,"year"] == 2010),],
                        idstud+sex+ses+language~item, value.var="value")
                        
# second, create the q matrix from the long format data frame
qMat <- unique(trends[ ,c("item","domain")])
qMat <- data.frame ( qMat[,"item", drop=FALSE], model.matrix(~domain-1, data = qMat))

# model 1
def1 <- defineModel(datW, items = -c(1:4), id="idstud", qMatrix = qMat, software="tam",
        HG.var = c("sex", "ses"), method="gauss")

# model 2 (misspecified)
def2 <- defineModel(datW, items = -c(1:6), id="idstud", qMatrix = qMat[-c(20,35,40),],
        software="tam", HG.var = c("sex", "ses", "language"), method="quadrature", nodes = 30)

# compare both outputs
compareDefineModelObjects (def1, def2)
}
