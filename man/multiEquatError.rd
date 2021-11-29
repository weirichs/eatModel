\name{multiEquatError}
\alias{multiEquatError}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{1pl linking errors for three measurement occasions}
\description{Function computes linking errors based on the 1pl linking procedure according
to \code{\link[sirt]{equating.rasch}} from the \code{sirt} package for three measurement
occasions.}
\usage{
multiEquatError (x1, x2, x3, difBound = 1, dependentDIF =FALSE, testletStr = NULL, verbose = TRUE )}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}} for the first measurement occasion. Alternatively,
a data.frame with two columns. First column: item identifier. Second column: item difficulty parameter.
}
  \item{x2}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}} for the second measurement occasion. Alternatively,
a data.frame with two columns. First column: item identifier. Second column: item difficulty parameter.
}
  \item{x3}{
%%     ~~Describe \code{file} here~~
The object returned by \code{\link{getResults}} for the third measurement occasion. Alternatively,
a data.frame with two columns. First column: item identifier. Second column: item difficulty parameter.
}
  \item{difBound}{
%%     ~~Describe \code{file} here~~
Defines the boundary. Items with absolute linking DIF greater than the boundary will
be removed from the linking procedure. 
}
  \item{dependentDIF}{
Logical. If DIF is assumed to be correlated between measurement occasions one and two
versus two and three, this covariance will be taken into account. If the assumption of
linear dependency does not seem reasonable, this option possibly results in overfitting
and underestimation of the true linking error.
}
  \item{testletStr}{
%%     ~~Describe \code{file} here~~
Optional: A data.frame with two columns. First column: testlet identifier, should be named \code{unit}.
Second column: item identifier, should be named \code{item}. Only necessary if linking errors should be
estimated via jackknife procedure. If items are nested in testlets, Monseur & Berezner, 2007 recommend
jackknife proceudres for linking error estimation in order to avoid bias.
}
  \item{verbose}{
%%     ~~Describe \code{file} here~~
Logical: print linking informations to console?
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
Karoline Sachse, with code borrowed from the sirt package
}
\references{
Monseur, C., & Berezner, A. (2007). The Computation of Equating Errors in International Surveys
in Education. \emph{Journal of Applied Measurement 8}(3): 323-335.
}
\examples{
data(trends)

################################################################################
###    Example 1: Linking for three measurement occasions and one domain     ###
################################################################################

# calibrate all three measurements, using a unidimensional model each (only reading)
results <- by(data = trends, INDICES = trends[,"year"], FUN = function (y){
           dat <- reshape2::dcast(subset ( y, domain == "reading"), idstud~item, value.var="value")
           def <- defineModel(dat=dat, items= -1, id="idstud", software="tam")
           run <- runModel(def)
           res <- getResults(run)
           return(res)})
lErrors <- multiEquatError(results[[1]], results[[2]], results[[3]], difBound = 0.64)

# dependent DIF
lErrDif <- multiEquatError(results[[1]], results[[2]], results[[3]], difBound = 0.64, dependentDIF =TRUE)

# Following code demonstrates the replacement of old linking error in the 'equatingList' object
# for further processing
# direct linking 1 to 3 (1 is reference)
it1     <- itemFromRes(results[[1]])
eq1.vs.3<- equat1pl(results[[3]], prmNorm = it1[,c("item", "est")], difBound = 0.64, iterativ = TRUE)

# replace 'direct' linking error with 'indirect' linking error from 'multiEquatError()'
eq1.vs.3<- replaceLinkingError (equatingList=eq1.vs.3, multiEquatError_output=lErrors)


################################################################################
###    Example 2: Linking for three measurement occasions and two domains    ###
################################################################################

# calibrate all three measurements, using a unidimensional model each
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

# dependent DIF
lErrDif2<- multiEquatError(results2[[1]], results2[[2]], results2[[3]], difBound = 0.64, dependentDIF =TRUE)

# direct linking 1 to 3 (1 is reference)
it1     <- itemFromRes(results2[[1]])
eq1.vs.3<- equat1pl(results2[[3]], prmNorm = it1[,c("item", "est")], difBound = 0.64, iterativ = TRUE)

# replace 'direct' linking error with 'indirect' linking error from 'multiEquatError()'
eq1.vs.3<- replaceLinkingError (equatingList=eq1.vs.3, multiEquatError_output=lErrors2)


################################################################################
###    Example 3: Jackknife linking for three measurement and two domains    ###
################################################################################

# we reuse the previously created object 'results2' here
allItems<- rbind(itemFromRes(results2[[1]]), itemFromRes(results2[[2]]), itemFromRes(results2[[3]]))
allItems<- unique(allItems[,"item"])
tl.frame<- data.frame (testlet = substr(allItems,1,3), item = allItems, stringsAsFactors = FALSE)

# repeat the linking from example 2
jkErr2  <- multiEquatError(results2[[1]], results2[[2]], results2[[3]], difBound = 0.64, testletStr = tl.frame)


################################################################################
###            Example 4: Linking for three item parameter lists             ###
################################################################################

# use three item parameter lists for equating
# borrow item parameter data from the 'sirt' package
data(data.pars1.rasch, package="sirt")

# use the first three item parameter lists
e1 <- subset ( data.pars1.rasch, study == "study1")[,c("item", "b")]
e2 <- subset ( data.pars1.rasch, study == "study2")[,c("item", "b")]
e3 <- subset ( data.pars1.rasch, study == "study3")[,c("item", "b")]
multiEquatError(e1, e2, e3)
}
