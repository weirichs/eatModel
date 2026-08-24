\name{transformToBista}
\alias{transformToBista}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transformation of item and person parameters to the Bista metric.}
\description{Function uses output of \code{\link{equat1pl}} to provide two data.frames
captured in a list, one for the item parameters on the metric defined for the evaluation
of educational standards in Germany (``bista metric''), and one for the person parameters
(PVs) on the bista metric .}
\usage{
transformToBista ( equatingList, refPop, cuts, weights = NULL, defaultM = 500, 
                   defaultSD = 100, q3bound = .2, roman = FALSE, vera = TRUE,
                   idVarName = NULL, years = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{equatingList}{
The object returned by \code{\link{equat1pl}}.
}
  \item{refPop}{
Optional: Data frame with at least three columns. First column indicates the domain name.
Note that this name must match the domain names in the output of \code{\link{getResults}}.
Second column contains the mean of the reference population. Third column contains
the standard deviation of the reference population. Fourth column optionally contains 
the transformed mean on the Bista metric of the reference population. Fifth column 
optionally contains the transformed standard deviation on the Bista metric of the 
reference population. If the fourth and fifth columns are missing, values will be 
defaulted to 500/100. If \code{refPop} is not specified, mean and SD will be computed
from the data, optionally using weights (if the \code{weights} argument is specified).
}
  \item{cuts}{
A named list with cut scores. Names of the list must match the domain names in the
output of \code{\link{getResults}}. Each element of the list is a named list with one or two
elements---the cut scores (in ascending order) and (optionally) the labels of the 
stages. The first element (cut scores) must be named \code{values}. The second element
(labels of the stages) must be named \code{labels}. See the examples of \code{\link{defineModel}}
for further details.
}
  \item{weights}{
Optional: a data.frame with two columns, first column is person identifier, second
columns is individual caseweight. Necessary for the transformation of linking 
error for (ordered) factors and/or if descriptives of the reference population 
should be computed directly from the data. See the examples of \code{\link{defineModel}}
for further details. 
}
  \item{defaultM}{
Mean of the reference population on the ``bista'' metric.
}
  \item{defaultSD}{
Standard deviation of the reference population in the ``bista'' metric.
}
  \item{q3bound}{
Define the absolute boundary of Q3 values which should be captured in the item parameter
list according to the guidelines of the ``Vergleichsarbeiten''.
}
  \item{roman}{
Logical: Use roman numbers for competence level column in the shortened item parameter
table dedicated for the ``Vergleichsarbeiten''? This only has a meaning if the argument
\code{vera} is set to \code{TRUE}.
}
  \item{vera}{
Logical: Prepare item parameter list according to the guidelines of the ``Vergleichsarbeiten''?
}
  \item{idVarName}{
Optional: character vector of individual student id. This is only to provide compatibility
with older package versions. Specification of this argument is only necessary if the function
gives an error.
}
  \item{years}{
Optional: numeric vector with two elements, indicating the both years of assessment.
Only necessary if an additional linking error object should be created.
}
}
\value{
A list with several data frames: the first one contains original and transformed
item parameters and competence levels. If multiple models were analyzed sequentially
via a model split, the first column of the data.frame indicates which model is being
considered. The following columns contain the item identifiers, information on which
dimension the item belongs to, response frequencies, and untransformed and transformed
item parameters. The ``estTransf'' column contains the parameter transformed to the
logit metric. ``estTransf625'' contains the parameter for a response probability of
62.5 percent instead of 50 percent. ``estTransfBista'' is the item parameter transformed
to the Bista metric. For dichotomous 1pl items, the transformation to ``estTransf625''
is as follows:
\deqn{
    \beta_{625}=\beta + log(\frac{.625}{1-.625})
}{%
  \beta_{625}=\beta + log(\frac{.625}{1-.625})
    }
For dichotomous 2pl models estimated using tam, the transformation is as follows:
\deqn{
    \beta_{625}=(\beta + log(\frac{.625}{1-.625})) \alpha^{-1}
}{%
  \beta_{625}=(\beta + log(\frac{.625}{1-.625})) \alpha^{-1}
    }
where \eqn{\alpha} is the item discrimination parameter. For dichotomous 2pl models
estimated using mirt, the transformation is as follows:
\deqn{
    \beta_{625}=\beta + log(\frac{.625}{1-.625}) \alpha^{-1}
}{%
  \beta_{625}=\beta + log(\frac{.625}{1-.625}) \alpha^{-1}
    }
The transformation to the metric for educational standards (``estTransfBista'') is
performed according to the following formula:
\deqn{
    \beta_{bista}=\frac{\beta_{625} - m_{ref}}{sd_{ref}} sd_{refTF} + m_{refTF}
}{%
  \beta_{bista}=\frac{\beta_{625} - m_{ref}}{sd_{ref}} sd_{refTF} + m_{refTF}
    }
where \eqn{m_{ref}} is the empirical mean of the reference population, \eqn{sd_{ref}}
is the empirical standard deviation of the reference population, \eqn{m_{refTF}} is
the transformed mean on the Bista metric of the reference population (usually 500),
and \eqn{sd_{ref}} is the transformed standard deviation on the Bista metric of the
reference population (usually 100). In the partial credit model, \eqn{\beta_{625}}
equals the thurstonian threshold of the step parameter according to a probabiliy
level of .625.
The second data.frame in the list contains the individual parameters (PVs). The
third data.frame contains the means and standard deviations for all dimensions
of the reference population. In terms of structure, this data.frame corresponds
to the object that the \code{transformToBista} function expects as the \code{refPop}
argument.
}
\author{
Sebastian Weirich
}
\examples{
# see example 5, 6, and 6a in the help file of defineModel()
}
