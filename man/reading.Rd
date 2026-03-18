\name{reading}

\docType{data}

\alias{reading}

\title{Dichotomous and partial credit responses for reading comprehension data obtained from large-scale assessment}

\description{
This data set contains fictional achievement scores of 2595 students in the long format along with item and person covariates. The data set is based on an actual pilot survey from 2012, but has been anonymized for use in the package.
}

\usage{data(reading)}

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
    \item{domain}{domain. As the data set only contains reading comprehension items and tasks, this variable is a constant wioth value 'reading'.}
  	\item{itemType}{response format of the item}
  	\item{type}{item type: The items that were part of the pilot study are labeled 'pilot'. The link items that provide the connection to the educational standards metric are called 'link'. Items from the IGLU study are called 'iglu'.}
  	\item{taskLabel}{short description of the task}
  	\item{textType}{brief description of the stimulus text on which the test tasks were based}
  	\item{position}{item position (i.e., block position)}
  	\item{item}{item identifier}
  	\item{valueSum}{The sum score for each student of all variables belonging to an item. Use this variable for partial credit model, i.e., polytomous responses.}
  	\item{valueAgg}{Aggregated item score. Use this variable for dichotomous responses.}
  	\item{valueMax}{maximum score which was possible for this item}
 }
}

\source{pilot study of a reading comprehension test of the Institute of Educational Progress}

\keyword{datasets}

