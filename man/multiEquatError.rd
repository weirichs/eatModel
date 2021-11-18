\name{multiEquatError}
\alias{multiEquatError}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{1pl linking errors for three measurements}
\description{Function computes linking errors based on the 1pl linking according
to \code{\link[sirt]{equating.rasch}} from the \code{sirt} package for three time
points of measurement.}
\usage{
multiEquatError (x1, x2, x3, difBound = 1 )}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}} for the first measurement. Alternatively,
a data.frame with two columns. First column: item identifier. Second column: item parameter.
}
  \item{x2}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}} for the second measurement. Alternatively,
a data.frame with two columns. First column: item identifier. Second column: item parameter.
}
  \item{x3}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}} for the third measurement. Alternatively,
a data.frame with two columns. First column: item identifier. Second column: item parameter.
}
  \item{difBound}{
%%     ~~Describe \code{file} here~~
Defines the boundary. Items with absolute linking DIF greater than the boundary will
be removed from the linking procedure. 
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
description needed
}
\author{
Karoline Sachse
}
\examples{
data(trends)
# calibrate all three measurements, using a unidimensional model each (only reading)
results <- by(data = trends, INDICES = trends[,"year"], FUN = function (y){
           dat <- reshape2::dcast(subset ( y, domain == "reading"), idstud~item, value.var="value")
           def <- defineModel(dat=dat, items= -1, id="idstud", software="tam")
           run <- runModel(def)
           res <- getResults(run)
           return(res)})
lErrors <- multiEquatError(results[[1]], results[[2]], results[[3]], difBound = 0.64)

# second example: calibrate all three measurements, using a unidimensional model each
# reading and listening are handled subsequently, using the model split function
results2<- by(data = trends, INDICES = trends[,"year"], FUN = function (y){
           qmat<- unique(y[ ,c("item","domain")])
           qmat<- data.frame ( qmat[,"item", drop=FALSE], model.matrix(~domain-1, data = qmat))
           spl <- splitModels ( qMatrix = qmat,  nCores = 1)
           dat <- reshape2::dcast(y, idstud~item, value.var="value")
           def <- defineModel(dat=dat, splittedModels = spl, id="idstud", software="tam")
           run <- runModel(def)
           res <- getResults(run)
           return(res)})
lErrors2<- multiEquatError(results2[[1]], results2[[2]], results2[[3]], difBound = 0.64)

# third example: use three item parameter lists for equating
# borrow item parameter data from the 'sirt' package
data(data.pars1.rasch, package="sirt")

# use the first three item parameter lists
e1 <- subset ( data.pars1.rasch, study == "study1")[,c("item", "b")]
e2 <- subset ( data.pars1.rasch, study == "study2")[,c("item", "b")]
e3 <- subset ( data.pars1.rasch, study == "study3")[,c("item", "b")]
multiEquatError(e1, e2, e3)
}
