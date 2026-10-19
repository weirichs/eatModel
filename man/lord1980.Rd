\name{lord1980}
\alias{lord1980}
\title{Applies the DIF rules of Lord (1980)}
\description{Given a data.frame containing the item parameters from a differential
item functioning (DIF) analysis, the rules for A-, B-, and C-DIF implemented by the
Educational Testing Service (ETS) can be applied to the items here. For each item,
the program indicates whether its DIF parameter meets the A, B, or C criterion.
Users generally do not need to use this function directly, as it is implemented
within \code{\link{getResults}}. The exception is DIF models performed as part
of partial credit analyses.}
\usage{
lord1980(dat, absDifCol, seDifCol, lbCol, ubCol, adb=NULL, sdb=NULL) }
\arguments{
  \item{dat}{
A data.frame with DIF item parameters. The most common use case here would be to
use the \code{"item*step"} sheet from the Conquest showfile, which is re-imported
using \code{\link{get.shw}} or \code{\link{get.shw.v5}}.
}
  \item{absDifCol}{
The column or column number in \code{dat} that contains the absolute DIF value.
}
  \item{seDifCol}{
The column or column number in \code{dat} that contains the standard error of the
DIF parameter.
}
  \item{lbCol}{
The column or column number in \code{dat} that contains the lower bound of the
confidence interval.
}
  \item{ubCol}{
The column or column number in \code{dat} that contains the upper bound of the
confidence interval.
}
  \item{adb}{
The absolute value of the DIF, which is to be used as a threshold. This is only
necessary if the user wishes to use a custom criterion for the DIF.
}
  \item{sdb}{
The significance value of the DIF, which is to be used as a threshold. This is only
necessary if the user wishes to use a custom criterion for the DIF.
}
}
\details{
The user-specific DIF indicator is 1 if the absolute DIF value is greater than
adb at a significance level of sdb.
}
\value{
The function returns the original data.frame, which is expanded by two columns:
an indicator column that, if desired, contains the DIF criterion specified by
the user; and a column containing the DIF according to the ETS criteria (A, C, and C).
}
\examples{
# see example 11 in the help file of defineModel()
}
         
