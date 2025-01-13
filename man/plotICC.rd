\name{plotICC}
\alias{plotICC}
\title{Plots item characteristic curves.}
\description{Function provides item characteristic plots for each item. To date, 
only dichotomous 1pl and 2pl models are supported.}
\usage{plotICC  ( resultsObj, defineModelObj, items = NULL, personPar = c("WLE", "EAP", "PV"),
       personsPerGroup = 30, pdfFolder = NULL, smooth = 7 )}
\arguments{
  \item{resultsObj}{
The object returned by \code{getResults}.
}
  \item{defineModelObj}{
The object returned by \code{defineModel}.
}
  \item{items}{
Optional: A vector of items for which the ICC should be plotted. If NULL, ICCs of all items
will be collected in a common pdf. The \code{pdfFolder} argument must not be NULL if the ICC of more
than one item should be plotted, i.e. if items is not NULL or a vector of length > 1.
}
  \item{personPar}{
Which person parameter should be used for plotting? To mimic the
behavior of the S3 plot method of \code{TAM}, use \code{"WLE"}.
}
  \item{personsPerGroup}{
%%     ~~Describe \code{file} here~~
Specifies the number of persons in each interval of the theta scale for dividing the 
persons in various groups according to mean EAP score. 
}
  \item{pdfFolder}{
%%     ~~Describe \code{file} here~~
Optional: A folder with writing access for the pdf file. Necessary only if ICCs for more
than one item should be plotted.
}
  \item{smooth}{
%%     ~~Describe \code{file} here~~
Optional: A parameter (integer value) for smoothing the plot. If the number of examinees is high, the
icc plot may become scratchy. \code{smooth} defines the maximum number of discrete nodes
across the theta scale for evaluating the icc. Higher values result in a less smooth icc. To mimic the
behavior of the S3 plot method of \code{TAM}, use the value 7.
}
}
\author{
Sebastian Weirich
}
\note{
This function is beta! Use with care...
}
\examples{
data(trends)
# choose only 2010
dat <- trends[which(trends[,"year"] == 2010),]
# choose reading
dat <- dat[which(dat[,"domain"] == "reading"),]

# first reshape the data set into wide format
datW <- reshape2::dcast(dat, idstud~item, value.var="value")

# defining the model: specifying q matrix is not necessary
mod1 <- defineModel(dat=datW, items= -1, id="idstud", software = "tam")

# run the model
run1 <- runModel(mod1)

# get the results
res1 <- getResults(run1)

# plot for one item 
plotICC  ( resultsObj = res1, defineModelObj = mod1, items = "T04_04")
}
