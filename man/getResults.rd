\name{getResults}
\alias{getResults}
\title{Collect all results from Conquest/TAM analysis into a common data frame}
\description{First the IRT model should be defined using \code{\link{defineModel}}. Afterwards,
call \code{\link{runModel}} with the argument returned by \code{\link{defineModel}} to start the estimation.
The last step then is to create a results frame using \code{getResults}. }
\usage{
getResults( runModelObj, overwrite = FALSE, Q3 = TRUE, q3theta = c("pv", "wle", "eap"), 
            q3MinObs = 0, q3MinType = c("singleObs", "marginalSum"), omitFit = FALSE, 
            omitRegr = FALSE, omitWle = FALSE, omitPV = FALSE, abs.dif.bound = 0.6, 
            sig.dif.bound = 0.3, p.value = 0.9, nplausible = NULL, ntheta = 2000,
            normal.approx = FALSE, samp.regr = FALSE, theta.model=FALSE, np.adj=8,
            group = NULL, beta_groups = TRUE, level = .95, n.iter = 1000,
            n.burnin = 500, adj_MH = .5, adj_change_MH = .05, refresh_MH = 50, 
            accrate_bound_MH = c(.45, .55),	sample_integers=FALSE, theta_init=NULL,
            print_iter = 20, verbose = TRUE, calc_ic=TRUE, omitUntil = 1, seed=NA) }
\arguments{
  \item{runModelObj}{
The object returned by \code{\link{runModel}}.
}
  \item{overwrite}{
Logical. Should result files be overwritten if exist?
}
  \item{Q3}{
Logical. Estimate the Q3 statistic according to Yen (1984)? Note: this is only
possible for uni-dimensional models. If \code{software == "tam"}, Q3 statistic
is estimated using the \code{\link[TAM]{tam.modelfit}} function. If \code{software == "Conquest"},
Q3 statistic is estimated using the \code{\link[sirt]{Q3}} function from the \code{sirt}
package.
}
  \item{q3theta}{
Specify whether the Q3 statistic should be estimated using PVs, WLEs or EAPs
as the theta variable.
}
  \item{q3MinObs}{
Q3 statistic might be untrustworthy if item covariance estimation is based on very
few observations. Define the minimum number of observation which should be fulfilled
for Q3 estimation. 
}
  \item{q3MinType}{
If \code{"singleObs"}, \code{q3MinObs} argument is based on the least number of observations in
the \eqn{2\times 2} 0/1 frequency table of item pairs. If \code{"marginalSum"}, \code{q3MinObs} argument is based on
the sum of marginals in the \eqn{2\times 2} 0/1 frequency table of item pairs.
}
  \item{omitFit}{
Logical. Should item fit values be included into the results?
}
  \item{omitRegr}{
Logical. Should regression parameters and their standard errors be included into the results?
}
  \item{omitWle}{
Logical. Should WLE estimates be included into the results?
}
  \item{omitPV}{
Logical. Should plausible values be included into the results?
}
  \item{abs.dif.bound}{
Applies only if DIF analyses are performed before. When DIF-Parameter are evaluated,
this specifies the critical value for absolute DIF. See the details section for further
information.
}
  \item{sig.dif.bound}{
Applies only if DIF analyses are performed before. When DIF-Parameter are evaluated,
this specifies the critical value for confidence interval DIF. See the details section for further
information.
}
  \item{p.value}{
Applies only if DIF analyses are performed before. When DIF-Parameter are evaluated,
this specifies the critical p-value for confidence interval DIF. See the details section for further
information.
}
 \item{nplausible}{
Applies only if \code{software = "tam"}: Number of plausible values to be drawn. Note: 
number of plausible values were already defined in \code{\link{defineModel}}, because
Conquest needs to know the number of PVs prior to estimation. In \code{TAM}, it 
is possible to redefine the number of plausible values and overwrite the definition
that was given in \code{\link{defineModel}}.
}
\item{ntheta}{
Applies only if \code{software = "tam"}. Following description is borrowed from the help 
file of \code{\link[TAM]{tam.pv}} from the \code{TAM} package: Number of ability nodes for
plausible value imputation. Note that in this function ability nodes are simulated 
for the whole sample, not for every person (contrary to the software Conquest).
}
\item{normal.approx}{
Applies only if \code{software = "tam"}. Following description is borrowed from the help 
file of \code{\link[TAM]{tam.pv}} from the \code{TAM} package: An optional logical indicating
whether the individual posterior distributions should be approximated by a normal 
distribution? The default is \code{FALSE}. In the case \code{normal.approx=TRUE}
(normal distribution approximation), the number of ability nodes \code{ntheta} can 
be substantially smaller than 2000, say 200 or 500. The normal approximation is 
implemented for unidimensional and multidimensional models.
}
\item{samp.regr}{
Applies only if \code{software = "tam"}. Following description is borrowed from the help 
file of \code{\link[TAM]{tam.pv}} from the \code{TAM} package: An optional logical indicating
whether regression coefficients should be fixed in the plausible value imputation or
also sampled from their posterior distribution? The default is \code{FALSE}. Sampled 
regression coefficients are obtained by nonparametric bootstrap.
}
\item{theta.model}{
Applies only if \code{software = "tam"}. Following description is borrowed from the help 
file of \code{\link[TAM]{tam.pv}} from the \code{TAM} package: Logical indicating whether the
theta grid from the \code{tamobj} object should be used for plausible value
imputation. In case of \code{normal.approx=TRUE}, this should be sufficient in many 
applications.
}
\item{np.adj}{
Applies only if \code{software = "tam"}. Following description is borrowed from the help 
file of \code{\link[TAM]{tam.pv}} from the \code{TAM} package: This parameter defines the
``spread'' of the random theta values for drawing plausible values when 
\code{normal.approx=FALSE}. If \eqn{s_{EAP}} denotes the standard deviation 
of the posterior distribution of theta (in the one-dimensional case), then theta
is simulated from a normal distribution with standard deviation \code{np.adj} 
times \eqn{s_{EAP}}.
}
\item{group}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. Optional 
vector of group identifiers. See the help page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{beta_groups}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. See the help 
page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{level}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. Confidence level
in bayesian approach. See the help page of\code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{n.iter}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. Number of 
iterations in the bayesian approach. See the help page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{n.burnin}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. Number of 
burn-in iterations in the bayesian approach. See the help page of \code{\link[TAM]{tam.pv.mcmc}} for
further details.
}
\item{adj_MH}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. See the help 
page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{adj_change_MH}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. See the help 
page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{refresh_MH}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. See the help 
page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{accrate_bound_MH}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. See the help 
page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{sample_integers}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. Logical
indicating whether weights for complete cases should be sampled in bootstrap. See the help
page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{theta_init}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. Optional matrix
with initial theta values. See the help page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{print_iter}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. See the help 
page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{verbose}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. See the help 
page of \code{\link[TAM]{tam.pv.mcmc}} for further details.
}
\item{calc_ic}{
Applies only if \code{software = "tam"} and \code{pvMethod = "bayesian"}. Logical
indicating whether information criteria should be computed. See the help
page of \code{tam.pv.mcmc} for further details.
}
  \item{omitUntil}{
Argument is passed to \code{\link{plotDevianceConquest}}: An optional value indicating
number of iterations to be omitted for plotting.
}
  \item{seed}{
Fixed simulation seed. This value is directly passed on to the \code{\link[TAM]{tam.fit}} function.
}
}
\details{
If \code{\link{defineModel}} was run with software Conquest, a path argument (\code{'dir'})
is necessary. The path argument is optional for software TAM. If \code{'dir'} was
specified, \code{getResults} additionally writes its output into the specified folder, 
using the \code{analysis.name} argument for file naming. Otherwise, \code{getResults}
only returnes the result data frame. 

If DIF analyses were performed before, the user can specify the criteria according to
which DIF should be interpreted or evaluated. By default, the ETS criteria (Zieky, 1993)
are used which classify DIF into three distinct categories, "A", "B", or "C". Small DIF ("A")
corresponds to absolute DIF values below .43 (no significance test is performed here); medium
DIF ("B") corresponds to absolute DIF values between .43 and .64 which are significantly
different from zero. High DIF ("C") corresponds to absolute DIF values above .64 which are
significantly different from .43 (DeMars, 2011; Monahan et al. 2007). Alternatively, the
three arguments \code{abs.dif.bound}, \code{sig.dif.bound}, and \code{p.value} allow to
specify user-defined dichotomous criteria. If items should be flagged as DIF, if the
absolute value increases 0.5 and significantly exceeds 0.1 at a alpha level of 0.05, use
\code{abs.dif.bound = 0.5} and \code{sig.dif.bound = 0.1} and \code{p.value = 0.95}.
}
\value{
A data frame in the long format with ten columns.
\item{model}{The name of the model (as specified by the user in \code{analysis.name}.}
\item{source}{The estimation software (i.e, conquest or TAM) }
\item{var1}{The variable name for which the corresponding value is given, i.e. its indicator. }
\item{var2}{Additional variable information if necessary.}
\item{type}{Type of coefficient (for example, random or fixed).}
\item{indicator.group}{The type of the group the corresponding variable belongs to.}
\item{group}{The group the corresponding variable belongs to. Note: group is nested within \code{indicator.group}.}
\item{par}{The type of the parameter.}
\item{derived.par}{Optionally: The derived parameter.}
\item{value}{The value of the corresponding estimate.}
}
\examples{
# see examples in the help file of defineModel()
}
\references{
DeMars, C. E. (2011). An analytic comparison of effect sizes for differential item functioning.
Applied Measurement in Education, 24 (3), 189-209. https://doi.org/10.1080/08957347.2011.580255

Monahan, P. O., McHorney, C. A., Stump, T. E. & Perkins, A. J. (2007). Odds ratio,
delta, ETS classification, and standardization measures of DIF magnitude for binary
logistic regression. Journal of Educational and Behavioral Statistics, 32 (1), 92-109.
https://doi.org/10.3102/1076998606298035

Zieky, M. (1993). Practical questions in the use of DIF statistics in item development. In P.
W. Holland & H. Wainer (Eds.), \emph{Differential item functioning} (pp. 337–347). Hillsdale, NJ:
Lawrence Erlbaum.
}