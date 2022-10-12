\name{trends}

\docType{data}

\alias{trends}

\title{Reading and listening achievement test data obtained from large-scale assessment at three times of measurement}

\description{
This data set contains fictional achievement scores of 13524 students in two domains (reading, listening) in the long format.
}

\usage{data(trends)}

\format{'data.frame':   404281 obs. of 17 variables:
  \describe{
    \item{year}{Year of evaluation}
    \item{idstud}{Student identifier}
    \item{idclass}{Class identifier}
    \item{wgt}{individual student weight}
    \item{jkzone}{jackknife zone (primary sampling unit)}
    \item{jkrep}{jackknife replicate}
    \item{country}{The country an examinee stems from}
    \item{language}{spoken language at home}
    \item{ses}{student's socio economical status}
    \item{sex}{student's sex}
  	\item{domain}{The domain the variable belongs to}
  	\item{booklet}{booklet identifier. equal booklet identifiers indicate equal booklets across years (assessment cycles)}
  	\item{block}{block identifier}
  	\item{task}{task identifier}
  	\item{item}{item identifier}
  	\item{pos}{position of the block within the booklet}
  	\item{value}{The response of the student to the item (0=incorrect; 1=correct)}
 }
}

\source{Simulated data}

\keyword{datasets}

\examples{
data(trends)
# number of students per year, country and domain
by(data=trends, INDICES = trends[,"year"], FUN = function(x) { tapply(x[,"idstud"], x[,c("country", "domain")], FUN = function(y){length(unique(y))})})
# number of items per year, country and domain
by(data=trends, INDICES = trends[,"year"], FUN = function(x) { tapply(x[,"item"], x[,c("country", "domain")], FUN = function(y){length(unique(y))})})

# students nested in classes?
lme4::isNested(trends[,"idstud"], trends[,"idclass"])
# items nested in tasks?
lme4::isNested(trends[,"item"], trends[,"task"])
# tasks nested in blocks? no, few tasks occur in more than one block
lme4::isNested(trends[,"task"], trends[,"block"])
# tasks nested in blocks for specific years?
by(data=trends, INDICES = trends[,"year"], FUN = function (y) {lme4::isNested(y[,"task"], y[,"block"]) })
# blocks nested in booklets?
lme4::isNested(trends[,"block"], trends[,"booklet"])
}



