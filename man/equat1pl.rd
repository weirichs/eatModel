\name{equat1pl}
\alias{equat1pl}
\title{1pl equating with optional elimination of linking DIF items}
\description{Function does the 1pl linking according to \code{\link[sirt]{equating.rasch}} from the \code{sirt} package.
Moreover, optional elimination of items with linking DIF is allowed and linking error may be estimated
via the jackknife method if testlets are specified.}
\usage{
equat1pl(results, prmNorm, item = NULL, domain = NULL, testlet = NULL, value = NULL, 
         excludeLinkingDif = TRUE, difBound = 1, iterativ = FALSE,
         method = c("Mean.Mean", "Haebara", "Stocking.Lord", "robust", "Haberman"), itemF = NULL,
         domainF = NULL, testletF = NULL, valueF = NULL, estimation=c("OLS", "BSQ", "HUB", "MED", "LTS", "L1", "L0"),
         b_trim=Inf, lts_prop=.5)}
\arguments{
  \item{results}{
The object returned by \code{\link{getResults}}. Alternatively, a data.frame with item
parameters of the focus group. In this case, additional arguments (\code{itemF}, 
\code{domainF}, and \code{valueF}) have to be defined.
}
  \item{prmNorm}{
Data frame with normed anchor item parameters. Data frame must have at least two
columns: items and item difficulties. Use the further arguments \code{item}, \code{domain},
\code{testlet} and \code{value} to define the column in which the corresponding parameter
can be found. If \code{item}, \code{domain}, \code{testlet} and \code{value} is NULL, \code{prmNorm} must
have only two columns: First column items, second column item difficulties. 
Column names then are irrelevant. 
}
  \item{item}{
Optional: Give the number or name of the item identifier column in prmNorm.
}
  \item{domain}{
Optional: Give the number or name of the domain name in prmNorm. Only necessary if
item identifiers are not unique in \code{prmNorm} (for example, if one item occurs several
times, with a global item parameter and a domain-specific item parameter. Domain names 
in prmNorm must match dimension names in the object returned by \code{\link{getResults}}.
}
  \item{testlet}{
Optional: Give the number or name of the testlet name in prmNorm. Only necessary if
linking errors should be estimated via the jackknife method. 
}
  \item{value}{
Optional: Give the number or name of the parameter column in \code{prmNorm}.
}
  \item{excludeLinkingDif}{
Logical. Should items with linking DIF excluded?
}
  \item{difBound}{
Defines the boundary. Items with absolute linking DIF greater than the boundary will
be removed from the linking procedure. 
}
  \item{iterativ}{
Logical. Should the exclusion of linking DIF items executed in an iterative loop?
(i.e. start with all items and compute linking constant, than remove the item with
the most pronounced DIF and compute linking constant, and so on, until no item is
left with |DIF|> \code{difBound}.
}
  \item{method}{
Linking method. If \code{"Mean.Mean"}, \code{"Haebara"}, or \code{"Stocking.Lord"},
the function \code{\link[sirt]{equating.rasch}} from the \code{sirt} package is called. If
\code{"robust"}, the function \code{\link[sirt]{linking.robust}} from the \code{sirt} package
is called. If \code{"Haberman"}, the function \code{\link[sirt]{linking.haberman}} from the
\code{sirt} package is called.
}
  \item{itemF}{
Optional: Give the number or name of the item column in results.
}
  \item{domainF}{
Optional: Give the number or name of the domain column in results.
}
  \item{testletF}{
Optional: Give the number or name of the testlet column in teh results object. Only necessary if
linking errors should be estimated via the jackknife method.
}
  \item{valueF}{
Optional: Give the number or name of the parameter column in results.
}
  \item{estimation}{
Applies only if \code{method} equals \code{"Haberman"}. Estimation method. See the
help page of \code{linking.haberman} from the \code{sirt} package for further details.
}
  \item{b_trim}{
Applies only if \code{method} equals \code{"Haberman"}. Trimming parameter for item slopes.
See the help page of \code{linking.haberman} from the \code{sirt} package for further details.
}
  \item{lts_prop}{
Applies only if \code{method} equals \code{"Haberman"}. Proportion of retained observations
in \code{"LTS"} regression estimation. See the help page of \code{linking.haberman} from the
\code{sirt} package for further details.
}
}
\value{
A list with equating information intended for further transformation by the \code{\link{transformToBista}}
function. 
}
\examples{
# see example 5, 6, and 6a in the help file of defineModel()
}
