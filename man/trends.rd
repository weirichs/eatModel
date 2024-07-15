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
  	\item{format}{item format}
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

# no overlapping student IDs between assessment cycles
ids <- by(data=trends, INDICES = trends[,"year"], FUN = function (x) {unique(x[,"idstud"])})
length(intersect(ids[[1]], ids[[2]]))
length(intersect(ids[[1]], ids[[3]]))
length(intersect(ids[[2]], ids[[3]]))

# sampling weights substantially differ between countries due to stratified sampling
eatTools::roundDF(do.call("rbind",  by(data=trends, INDICES = trends[,c("year", "country")], FUN = function (x) {data.frame ( trends[1,c("year", "country")], eatTools::descr(x[!duplicated(x[,"idstud"]),"wgt"])[,c("Minimum", "Maximum", "Mean", "Median", "SD")], stringsAsFactors = FALSE)})), digits = 3)

# which booklets occur in which assessment cycles?
# see, for example: Bo01 only occurs 2010; Bo02 occurs 2010, 2015, and 2022; Bo83 occurs 2015 and 2020
reshape2::dcast(do.call("rbind", by(data=trends, INDICES = trends[,"year"], FUN = function (x) {data.frame ( x[1,"year", drop=FALSE], table(x[!duplicated(x[,"idstud"]),"booklet"]), stringsAsFactors = FALSE)})), year~Var1, value.var = "Freq")

# which reading tasks occur in which assessment cycles?
# see, for example: T01 occurs 2010, 2015, and 2022; T27 only occurs 2020
reshape2::dcast(do.call("rbind", by(data=subset(trends,domain=="reading"), INDICES = subset(trends,domain=="reading")[,"year"], FUN = function (x) {data.frame ( x[1,"year", drop=FALSE], table(x[!duplicated(x[,"idstud"]),"task"]), stringsAsFactors = FALSE)})), year~Var1, value.var = "Freq")

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



