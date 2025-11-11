\name{reading}

\docType{data}

\alias{reading}

\title{Partial credit responses for reading comprehension data obtained from large-scale assessment}

\description{
This data set contains fictional achievement scores of 2595 students in the long format.
}

\usage{data(trends)}

\format{'data.frame':    64830 obs. of  22 variables:
  \describe{
    \item{taskID}{task ID}
    \item{bookletID}{booklet ID}
    \item{uniqueID}{unique student identifier}
    \item{ID}{student identifier}
    \item{sex}{student's sex}
    \item{language}{student's language spoken at home}
    \item{age}{student's age in years}
    \item{country}{the country the student stems from (anonymized)}
    \item{domain}{domain}
    \item{sex}{student's sex}
  	\item{domain}{The domain the variable belongs to}
  	\item{itemType}{response format of the item}
  	\item{type}{item type}
  	\item{taskLabel}{short description of the task}
  	\item{position}{item position (i.e., block position)}
  	\item{item}{item identifier}
    \item{pos}{position of the block within the booklet}
  	\item{valueSum}{The sum score for each student of all variables belonging to an item. Use this variable for partial credit model, i.e., polytomous responses.}
  	\item{valueAgg}{aggregated item score. Use this variable for dichotomous responses.}
  	\item{valueMax}{maximum score which was possible for this item}
 }
}

\source{Simulated data}

\keyword{datasets}

