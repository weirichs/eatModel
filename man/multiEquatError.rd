\name{multiEquatError}
\alias{multiEquatError}
\title{1pl linking errors for three measurement occasions}
\description{Function combines linking errors based on the 1pl linking procedure according
to \code{\link{equat1pl}} for three measurement occasions.}
\usage{
multiEquatError (eq.1_2, eq.2_3, eq.1_3, dependentDIF =FALSE, verbose=TRUE)}
\arguments{
  \item{eq.1_2}{
The object returned by \code{\link{equat1pl}} for linking the first and second measurement occasion.
}
  \item{eq.2_3}{
The object returned by \code{\link{equat1pl}} for linking the second and third measurement occasion.
}
  \item{eq.1_3}{
The object returned by \code{\link{equat1pl}} for linking the first and third measurement occasion.
}
  \item{dependentDIF}{
Logical. If DIF is assumed to be correlated between measurement occasions one and two
versus two and three, this covariance will be taken into account. If the assumption of
linear dependency does not seem reasonable, this option possibly results in overfitting
and underestimation of the true linking error.
}
  \item{verbose}{
Logical: print linking informations to console?
}
}
\value{
A list of data.frames with linking errors, according to the number of domains. The output is
intended to be used as argument of the replaceLinkingError function.
}
\author{
Sebastian Weirich
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

# calibrate all three measurements, using a unidimensional model each
# (only reading, no testlets, no jackknife)
results <- by(data = trends, INDICES = trends[,"year"], FUN = function (y){
           dat <- reshape2::dcast(subset ( y, domain == "reading"), idstud~item, value.var="value")
           def <- defineModel(dat=dat, items= -1, id="idstud", software="tam")
           run <- runModel(def)
           res <- getResults(run)
           return(res)})
items   <- lapply(results, itemFromRes)
eq.1_2  <- equat1pl(items[[1]][,c("item", "est")], items[[2]][,c("item", "est")],difBound = 0.64, iterativ = TRUE)
eq.2_3  <- equat1pl(items[[2]][,c("item", "est")], items[[3]][,c("item", "est")],difBound = 0.64, iterativ = TRUE)
eq.1_3  <- equat1pl(items[[1]][,c("item", "est")], items[[3]][,c("item", "est")],difBound = 0.64, iterativ = TRUE)
lErrors <- multiEquatError(eq.1_2, eq.2_3, eq.1_3)

# dependent DIF
lErrDif <- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)

# Following code demonstrates the replacement of old linking error in the 'equatingList' object
# for further processing
# direct linking 1 to 3 (1 is reference)
eq1.vs.3<- equat1pl(results[[3]], prmNorm = items[[1]][,c("item", "est")], difBound = 0.64, iterativ = TRUE)

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
eq.1_2  <- equat1pl(results2[[1]], itemFromRes(results2[[2]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE)
eq.2_3  <- equat1pl(results2[[2]], itemFromRes(results2[[3]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE)
eq.1_3  <- equat1pl(results2[[1]], itemFromRes(results2[[3]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE)
lErrors2<- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)

# direct linking 1 to 3 (1 is reference)
it1     <- itemFromRes(results2[[1]])
eq1.vs.3<- equat1pl(results2[[3]], prmNorm = it1[,c("item", "est")], difBound = 0.64, iterativ = TRUE)

# replace 'direct' linking error with 'indirect' linking error from 'multiEquatError()'
eq1.vs.3<- replaceLinkingError (equatingList=eq1.vs.3, multiEquatError_output=lErrors2)


################################################################################
###    Example 3: Jackknife linking for three measurement and two domains    ###
################################################################################

# we reuse the previously created object 'results2' here
eq.1_2  <- equat1pl(results2[[1]], itemFromRes(results2[[2]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE)
eq.2_3  <- equat1pl(results2[[2]], itemFromRes(results2[[3]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE)
eq.1_3  <- equat1pl(results2[[1]], itemFromRes(results2[[3]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE)
lErrors3<- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)


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

# pairwise linking
eq.1_2 <- equat1pl(e1, e2,difBound = 0.64, iterativ = TRUE)
eq.2_3 <- equat1pl(e2, e3,difBound = 0.64, iterativ = TRUE)
eq.1_3 <- equat1pl(e1, e3,difBound = 0.64, iterativ = TRUE)
multiEquatError(eq.1_2, eq.2_3, eq.1_3)
}
