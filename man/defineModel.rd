\name{defineModel}
\alias{defineModel}
\title{Prepares IRT analysis for Conquest and TAM}
\description{Facilitates data analysis using the software Conquest and/or TAM. It automatically
checks data for IRT consistency, generates Conquest syntax, label, anchor and data files or
corresponding TAM call for a single model specified by several arguments in R. Finally, an
R object is created which contain the required input for Conquest or TAM. To start the estimation,
call \code{\link{runModel}} with the argument returned by \code{defineModel}.}
\usage{
defineModel (dat, items, id, splittedModels = NULL,
   irtmodel = c("1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", "2PL.groups",
   "GPCM.design", "3PL"), qMatrix=NULL, DIF.var=NULL, HG.var=NULL, group.var=NULL,
   weight.var=NULL, anchor = NULL, domainCol=NULL, itemCol=NULL, valueCol=NULL,
   check.for.linking = TRUE, minNperItem = 50, removeMinNperItem = FALSE,
   boundary = 6, remove.boundary = FALSE, remove.no.answers = TRUE,
   remove.no.answersHG = TRUE, remove.missing.items = TRUE,
   remove.constant.items = TRUE, remove.failures = FALSE,
   remove.vars.DIF.missing = TRUE, remove.vars.DIF.constant = TRUE,
   verbose=TRUE, software = c("conquest","tam"), dir = NULL, analysis.name,
   schooltype.var = NULL, model.statement = "item",  compute.fit = TRUE,
   pvMethod = c("regular", "bayesian"), fitTamMmlForBayesian = TRUE,
   n.plausible=5, seed = NULL,
   conquest.folder= NULL,
   constraints=c("cases","none","items"), std.err=c("quick","full","none"),
   distribution=c("normal","discrete"),
   method=c("gauss", "quadrature", "montecarlo", "quasiMontecarlo"),
   n.iterations=2000, nodes=NULL, p.nodes=2000, f.nodes=2000,converge=0.001,
   deviancechange=0.0001, equivalence.table=c("wle","mle","NULL"), use.letters=FALSE,
   allowAllScoresEverywhere = TRUE, guessMat = NULL, est.slopegroups = NULL,
   fixSlopeMat = NULL, slopeMatDomainCol=NULL, slopeMatItemCol=NULL,
   slopeMatValueCol=NULL, progress = FALSE, Msteps = NULL, increment.factor=1 ,
   fac.oldxsi=0, export = list(logfile = TRUE, systemfile = FALSE, history = TRUE,
   covariance = TRUE, reg_coefficients = TRUE, designmatrix = FALSE))}
\arguments{
  \item{dat}{
A data frame containing all variables necessary for analysis.
}
  \item{items}{
Names or column numbers of variables with item responses. Item response values must
be numeric (i.e. 0, 1, 2, 3 ... ). Character values (i.e. A, B, C ... or a, b, c, ...)
are not allowed. Class of item columns are expected to be numeric or integer.
Columns of class \code{character} will be transformed.
}
  \item{id}{
Name or column number of the identifier (ID) variable.
}
  \item{splittedModels}{
Optional: Object returned by \code{\link{splitModels}}. Definition for multiple model handling.
}
  \item{irtmodel}{
Specification of the IRT model. The argument corresponds to the \code{irtmodel}
argument of \code{\link[TAM]{tam.mml}}. See the help page of \code{\link[TAM]{tam.mml}} for further details.
}
  \item{qMatrix}{
Optional: A named data frame indicating how items should be grouped to dimensions. The
first column contains the unique names of all items and should be named ``item''. The other
columns contain dimension definitions and should be named with the respective
dimension names. A positive value (e.g., 1 or 2 or 1.4) indicates the loading weight
with which an item loads on the dimension, a value of 0 indicates that the respective
item does not load on this dimension. If no q matrix is specified by the user, an
unidimensional structure is assumed.
}
  \item{DIF.var}{
Name or column number of one grouping variable for which differential item
functioning analysis is to be done.
}
  \item{HG.var}{
Optional: Names or column numbers of one or more context variables (e.g., sex, school).
These variables will be used for latent regression model in Conquest or TAM.
}
  \item{group.var}{
Optional: Names or column numbers of one or more grouping variables. Descriptive
statistics for WLEs and Plausible Values will be computed separately for each
group in Conquest.
}
  \item{weight.var}{
Optional: Name or column number of one weighting variable.
}
  \item{anchor}{
Optional: A named data frame with anchor parameters. The first column contains the
names of all items which are used to be anchor items and may be named item.
The second column contains anchor parameters. Anchor items can be a subset of the
items in the dataset and vice versa. If the data frame contains more than two
columns, columns must be named explicitly using the following arguments
\code{domainCol}, \code{itemCol}, and \code{valueCol}.
}
  \item{domainCol}{
Optional: Only necessary if the \code{anchor} argument was used to define
anchor parameters. Moreover, specifying \code{domainCol} is only necessary, if the
item identifiers in \code{anchor} are not unique---for example, if a specific item
occurs with two parameters, one domain-specific item parameter and one additional
``global'' item parameter. The domain column than must specify which parameter
belongs to which domain.
}
  \item{itemCol}{
Optional: Only necessary if the \code{anchor} argument was used to define
anchor parameters. Moreover, specifying \code{itemCol} is only necessary, if the
\code{anchor} data frame has more than two columns. The \code{itemCol} column than
must specify which column contains the item identifier.
}
  \item{valueCol}{
Optional: Only necessary if the \code{anchor} argument was used to define
anchor parameters. Moreover, specifying \code{valueCol} is only necessary, if the
\code{anchor} data frame has more than two columns. The \code{valueCol} column than
must specify which column contains the item parameter values.
}
  \item{check.for.linking}{
A logical value indicating whether the items in dataset are checked for being
connected with each other via design.
}
  \item{minNperItem}{
Numerical: A message is printed on console if an item has less valid values than the number
defined in \code{minNperItem}.
}
  \item{removeMinNperItem}{
Logical: Remove items with less valid responses than defined in \code{minNperItem}?
}
  \item{boundary}{
Numerical: A message is printed on console if a subject has answered less than the number of items
defined in boundary.
}
  \item{remove.boundary}{
Logical: Remove subjects who have answered less items than defined in the \code{boundary} argument?
}
  \item{remove.no.answers}{
Logical: Should persons without any item responses being removed prior to analysis?
}
  \item{remove.no.answersHG}{
Logical: Should persons without any responses on any background variable being removed prior to analysis?
}
  \item{remove.missing.items}{
Logical: Should items without any item responses being removed prior to analysis?
}
  \item{remove.constant.items}{
Logical: Should items without variance being removed prior to analysis?
}
  \item{remove.failures}{
Logical: Should persons without any correct item response (i.e., only \dQuote{0} responses) being removed prior to analysis?
}
  \item{remove.vars.DIF.missing}{
Logical: Applies only in DIF analyses. Should items without any responses in at least one
DIF group being removed prior to analyses? Note: Conquest may crash if these items
remain in the data.
}
  \item{remove.vars.DIF.constant}{
Logical: Applies only in DIF analyses. Should items without variance in at least one
DIF group being removed prior to analyses? Note: Conquest may crash if these items
remain in the data.
}
  \item{verbose}{
A logical value indicating whether messages are printed on the R console.
}
  \item{software}{
The desired estimation software for the analysis.
}
  \item{dir}{
The directory in which the output will be written to. If \code{software = "conquest"},
\code{dir} must be specified. If \code{software = "tam"}, \code{dir} is not mandatory.
}
  \item{analysis.name}{
A character string specifying the analysis name. If \code{software = "conquest"},
\code{analysis.name} must be specified. All Conquest input and output files will
named \code{analysis.name} with their corresponding extensions. If \code{software = "tam"},
\code{analysis.name} is not mandatory. In the case of multiple models estimation,
\code{split.models} automatically defines \code{analysis.name} for each model.
}
  \item{schooltype.var}{
Optional: Name or column number of the variable indicating the school type (e.g.
academic track, non-academic track). Only necessary if \emph{p} values should be
computed for each school type separately.
}
  \item{model.statement}{
Optional: Applies only if \code{software = "conquest"}. A character string given the model
statement in the Conquest syntax. If omitted, the statement is generated automatically
with respect to the defined model.
}
  \item{compute.fit}{
Applies only if \code{software = "conquest"}. Compute item fit statistics?
}
  \item{pvMethod}{
Applies only if \code{software = "tam"}: Specifies whether PVs should be drawn regularly
or using a Bayesian algorithm.
}
  \item{fitTamMmlForBayesian}{
Logical, applies only if \code{software = "tam"}: If PVs are drawn using a Bayesian
algorithm, it is not necessary to fit the model via \code{\link[TAM]{tam.mml}} before. \code{fitTamMmlForBayesian}
specifies whether the model should be fitted before though. See the help page of \code{tam.pv.mcmc}
for further details.
}
  \item{n.plausible}{
The number of plausible values which are to be drawn from the conditioning model.
}
  \item{seed}{
Optional: Set seed value for analysis. The value will be used in Conquest syntax file ('set seed'-statement,
see conquest manual, p. 225) or in TAM (control$seed). Note that seed only occurs for stochastic integration.
}
  \item{conquest.folder}{
Applies only if \code{software = "conquest"}. A character string with path and name
of the Conquest console, for example \code{"c:/programme/conquest/console_Feb2007.exe"}.
If \code{NULL}, function tries to retrieve Conquest console from
\code{"i:/Methoden/00_conquest_console/console_Feb2007.exe"}. If the path is'nt valid,
the user is requested to specify the path of the console.
}
  \item{constraints}{
A character string specifying how the scale should be constrained. Possible options
are \code{"cases"} (default), \code{"items"} and \code{"none"}. When anchor parameter are specified in
anchor, constraints will be set to \code{"none"} automatically. In \code{TAM} the option
\code{"none"} is not allowed. (See the help file of \code{\link[TAM]{tam.mml}} for further details.)
}
  \item{std.err}{
Applies only if \code{software = "conquest"}. A character string specifying which
type of standard error should be estimated. Possible options are \code{"full"}, \code{"quick"}
(default) and \code{"none"}. See Conquest manual pp.167 for details on standard error estimation.
}
  \item{distribution}{
Applies only if \code{software = "conquest"}. A character string indicating the
a priori trait distribution. Possible options are \code{"normal"} (default) and \code{"discrete"}.
See Conquest manual pp.167 for details on population distributions.
}
  \item{method}{
A character string indicating which method should be used for numerical or stochastic
integration. Possible options are \code{"gauss"} (Gauss-Hermmite quadrature: default),
\code{"quadrature"} (Bock/Aitken quadrature) and \code{"montecarlo"}. See Conquest manual
pp.167 for details on these methods. When using \code{software = "tam"}, \code{"gauss"} and
\code{"quadrature"} essentially leads to numerical integration, i.e TAM is called with
\code{control$snodes = 0} and with \code{control$nodes = seq(-6,6,len=nn)}, where
\code{nn} equals the number of nodes specified in the \code{nodes} argument of
\code{defineModel} (see below). When using \code{software = "tam"}, \code{"montecarlo"}
leads to calling TAM with \code{control$QMC = FALSE} and \code{snodes = nn}, where
\code{nn} equals the number of nodes specified in the \code{nodes} argument of
\code{defineModel}. When using \code{software = "tam"}, \code{"quasiMontecarlo"}
leads to calling TAM with \code{control$QMC = TRUE} and \code{snodes = nn}, where
\code{nn} equals the number of nodes specified in the \code{nodes} argument of
\code{defineModel}. To met the \code{software = "tam"} default (numerical integration), use
\code{software="tam", nodes = 21}.
}
  \item{n.iterations}{
%%     ~~Describe \code{dif.term} here~~
An integer value specifying the maximum number of iterations for which estimation
will proceed without improvement in the deviance.
}
  \item{nodes}{
%%     ~~Describe \code{dif.term} here~~
An integer value specifying the number of nodes to be used in the analysis. The
default value is 20. When using \code{software = "tam"}, the value specified here
leads to calling TAM with \code{nodes = 20} AND \code{snodes = 0} if \code{"gauss"} or
\code{"quadrature"} was used in the \code{method} argument. If \code{"montecarlo"} or \code{"quasiMontecarlo"}
was used in the \code{method} argument, the value specified here leads to calling TAM with
\code{control$snodes = 20}. For numerical integration, for example,
\code{method = "gauss"} and \code{nodes = 21} (TAM default) may be appropriate.
For quasi monte carlo integration, \code{method = "quasiMontecarlo"} and \code{nodes = 1000}
may be appropriate (TAM authors recommend to use at least 1000 nodes).
}
  \item{p.nodes}{
%%     ~~Describe \code{dif.term} here~~
Applies only if \code{software = "conquest"}. An integer value specifying the
number of nodes that are used in the approximation of the posterior distributions,
which are used in the drawing of plausible values and in the calculation of EAP
estimates. The default value is 2000.
}
  \item{f.nodes}{
%%     ~~Describe \code{dif.term} here~~
Applies only if \code{software = "conquest"}. An integer value specifying the
number of nodes that are used in the approximation of the posterior distributions
in the calculation of fit statistics. The default value is 2000.
}
  \item{converge}{
%%     ~~Describe \code{dif.term} here~~
An integer value specifying the convergence criterion for parameter estimates.
The estimation will terminate when the largest change in any parameter estimate
between successive iterations of the EM algorithm is less than converge. The
default value is 0.001.
}
  \item{deviancechange}{
An integer value specifiying the convergence criterion for the deviance. The
estimation will terminate when the change in the deviance between successive
iterations of the EM algorithm is less than deviancechange. The default value
is 0.0001.
}
  \item{equivalence.table}{
Applies only if \code{software = "conquest"}. A character string specifying the
type of equivalence table to print. To be more precise, the person estimator that
is to be used to create the table must be specified here. Possible options are
\code{"wle"} (default), \code{"mle"} and NULL. If NULL, no equivalence table is
computed at all.
}
  \item{use.letters}{
Applies only if \code{software = "conquest"}. A logical value indicating whether
item response values should be coded as letters. This option can be used in partial
credit models comprising items with more than 10 categories to avoid response columns
with width 2 in Conquest.
}
  \item{allowAllScoresEverywhere}{
%%     ~~Describe \code{dif.term} here~~
Applies only if \code{software = "Conquest"}. Defines score statement generation
in multidimensional polytomous models. Consider two dimensions, `reading' and `listening'.
In `reading', values 0, 1, 2, 3 occur. In `listening', values 1, 2, 3, 4 occur. If \code{TRUE},
values 0, 1, 2, 3, 4 are defined for both dimensions. Otherwise, values 0, 1, 2, 3
are defined for `reading', values 1, 2, 3, 4 are defined for `listening'.
}
  \item{guessMat}{
%%     ~~Describe \code{dif.term} here~~
Applies only if \code{software = "tam"} for 3PL models. A named data frame
with two columns indicating for which items a common guessing parameter should
be estimated. The first column contains the names of all items in the analysis
and should be named \code{"item"}. The second column is numerical (integer values
recommended) and allocates the items to groups. For each group of items, a
separate guessing parameter is estimated. If the value in the second columns equals
zero, the guessing parameter is fixed to zero.
}
  \item{est.slopegroups}{
%%     ~~Describe \code{dif.term} here~~
Applies only if \code{software = "tam"} for 2PL models. Optionally, a named data frame
with two columns indicating for which items a common discrimination parameter should
be estimated. The first column contains the names of all items in the analysis
and should be named \code{"item"}. The second column is numerical (integer values
recommended) and allocates the items to groups. For each group of items, a
separate discrimination parameter is estimated. Without specifying \code{est.slopegroups},
a discrimination parameter for each item is estimated.
}
  \item{fixSlopeMat}{
%%     ~~Describe \code{dif.term} here~~
Applies only if \code{software = "tam"} for 2PL models. Optionally, a named data frame
with two columns indicating for which items a fixed discrimination should be assumed.
The first column contains the names of the items which discrimination should be fixed.
Note that item indicators should be unique---if not, use further arguments \code{slopeMatDomainCol},
\code{slopeMatItemCol} and \code{slopeMatValueCol}. The second column is numerical and contains
the discrimination value. Note: To date, this works only for between item dimensionality models.
Within item dimensionality models must be specified directly in TAM, using the \code{B.fixed}
argument of \code{\link[TAM]{tam.mml}}. Items which discrimation should be estimated should not occur in this data frame.
}
  \item{slopeMatDomainCol}{
%%     ~~Describe \code{dif.term} here~~
Optional: Only necessary if the \code{fixSlopeMat} argument was used to define fixed slope
parameters. Moreover, specifying \code{slopeMatDomainCol} is only necessary, if the item
identifiers in \code{fixSlopeMat} are not unique---for example, if a specific item occurs
with two slope parameters, one domain-specific item slope parameter and one additional
``global'' item parameter. The domain column than must specify which parameter belongs to which domain.
}
  \item{slopeMatItemCol}{
%%     ~~Describe \code{dif.term} here~~
Optional: Only necessary if the \code{fixSlopeMat} argument was used to define fixed slope
parameters. Moreover, specifying \code{itemCol} is only necessary, if the \code{fixSlopeMat}
data frame has more than two columns. The \code{itemCol} column than must specify which column
contains the item identifier.
}
  \item{slopeMatValueCol}{
%%     ~~Describe \code{dif.term} here~~
Optional: Only necessary if the \code{fixSlopeMat} argument was used to define slope parameters.
Moreover, specifying \code{valueCol} is only necessary, if the \code{fixSlopeMat} data frame has
more than two columns. The \code{valueCol} column than must specify which column contains the item parameter values.
}
  \item{progress}{
%%     ~~Describe \code{dif.term} here~~
Applies only if \code{software = "tam"}. Print estimation progress messages on console?
}
  \item{Msteps}{
Number of M steps for item parameter estimation. A high value of M steps could be helpful in cases of non-convergence.
The default value is 4; the default for 3pl models is set to 10.
}
  \item{increment.factor}{
%%     ~~Describe \code{dif.term} here~~
Applies only if \code{software = "tam"}. Should only be varied if the model does not converge.
See help page of \code{\link[TAM]{tam.mml}} for further details.
}
  \item{fac.oldxsi}{
%%     ~~Describe \code{dif.term} here~~
Applies only if \code{software = "tam"}. Should only be varied if the model does not converge.
See help page of \code{\link[TAM]{tam.mml}} for further details.
}
  \item{export}{
%%     ~~Describe \code{dif.term} here~~
Applies only if \code{software = "conquest"}. Specifies which additional files should be written
on hard disk.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list which contains information about the desired estimation. The list is intended for
further processing via \code{\link{runModel}}. Structure of the list varies depending on
whether multiple models were called using \code{\link{splitModels}} or not. If
\code{\link{splitModels}} was called, the number of elements in the list equals
the number of models defined via \code{\link{splitModels}}. Each element in the list is
a list with various elements:
  \item{software}{
Character string of the software which is intended to use for the further estimation,
i.e. \code{"conquest"} or \code{"tam"}
}
  \item{qMatrix}{
The Q matrix allocating items to dimensions.
}
  \item{all.Names}{
Named list of all relevant variables of the data set.
}
  \item{dir}{
Character string of the directory the results are to be saved.
}
  \item{analysis.name}{
Character string of the analysis' name.
}
  \item{deskRes}{
Data frame with descriptives (e.g., p values) of the test items.
}
  \item{discrim}{
Data frame with item discrimination values.
}
  \item{perNA}{
The person identifiers of examinees which are excluded from the analysis due to solely missing values.
}
  \item{per0}{
The person identifiers of examinees which have solely false responses. if \code{remove.failues} was set
to be TRUE, these persons are excluded from the data set.
}
  \item{perA}{
The person identifiers of examinees which have solely correct responses.
}
  \item{perExHG}{
The person identifiers of examinees which are excluded from the analysis due to missing values on explicit variables.
}
  \item{itemsExcluded}{
Character string of items which were excluded, for example due to zero variance or solely missing values.
}
If \code{software == "conquest"}, the output additionally includes the following elements:
  \item{input}{
Character string of the path with Conquest input (cqc) file.
}
  \item{conquest.folder}{
Character string of the path of the conquest executable file.
}
  \item{model.name}{
Character string of the model name.
}
If \code{software == "tam"}, the output additionally includes the following elements:
  \item{anchor}{
Optional: data frame of anchor parameters (if anchor parameters were defined).
}
  \item{daten}{
The prepared data for TAM analysis.
}
  \item{irtmodel}{
Character string of the used IRT model.
}
  \item{est.slopegroups}{
Applies for 2pl modeling. Information about which items share a common slope parameter.
}
  \item{guessMat}{
Applies for 3pl modeling. Information about which items share a common guessing parameter.
}
  \item{control}{
List of control parameters for TAM estimation.
}
  \item{n.plausible}{
Desired number of plausible values.
}
}
\author{
Sebastian Weirich
}
\examples{
################################################################################
###               Preparation: necessary for all examples                    ###
################################################################################

# load example data
# data set 'trends' contains item response data for three measurement occasions
# in a single data.frame
data(trends)

# first reshape the data for the first time of measurement set into wide format
datW <- reshape2::dcast(trends[which(trends[,"year"] == 2010),],
                        idstud+sex+ses+language~item, value.var="value")

# second, create the q matrix from the long format data frame
qMat <- unique(trends[ ,c("item","domain")])
qMat <- data.frame ( qMat[,"item", drop=FALSE], model.matrix(~domain-1, data = qMat))


################################################################################
###                Example 1: Unidimensional Rasch Model                     ###
################################################################################

# Example 1: define and run a unidimensional Rasch model with all variables in dataset
# (reading and listening together) using "TAM".

# defining the model: specifying q matrix is not necessary
mod1 <- defineModel(dat=datW, items= -c(1:4), id="idstud", analysis.name = "unidim",
        software="tam" )

# run the model
run1 <- runModel(mod1)

# get the results
res1 <- getResults(run1)

# extract the item parameters from the results object
item <- itemFromRes ( res1 )


################################################################################
###       Example 1a: Unidimensional Rasch Model with DIF estimation         ###
################################################################################

# nearly the same procedure as in example 1. Using 'sex' as DIF variable
# note that 'sex' is a factor variable here. Conquest needs all explicit variables
# to be numeric. Variables will be automatically transformed to numeric by
# 'defineModels'. However, it might be the better idea to transform the variable
# manually.
datW[,"sexNum"] <- car::recode ( datW[,"sex"] , "'male'=0; 'female'=1", as.factor = FALSE)

# as we have defined a new variable ('sexNum') in the data, it is a good idea
# to explicitly specify item columns ... instead of saying 'items= -c(1:3)' which
# means: Everything except column 1 to 3 are item columns
items<- grep("^T[[:digit:]]{2}", colnames(datW))
mod1a<- defineModel(dat=datW, items= items, id="idstud", DIF.var = "sexNum",
        analysis.name = "unidimDIF", software="tam",fac.oldxsi=.1, increment.factor=1.07)

# run the model
run1a<- runModel(mod1a)

# get the results
res1a<- getResults(run1a)


################################################################################
###        Example 2a: Multidimensional Rasch Model with anchoring           ###
################################################################################

# Example 2a: running a multidimensional Rasch model on a subset of items with latent
# regression. Use item parameter from the first model as anchor parameters

# read in anchor parameters from the results object of the first example
aPar <- itemFromRes ( res1 )
aPar <- aPar[,c("item", "est")]

# defining the model: specifying q matrix now is necessary.
# Please note that all latent regression variables have to be of class numeric.
# If regression variables are factors, dummy variables automatically will be used.
# (This behavior is equivalent as in lm() for example.)
mod2a<- defineModel(dat=datW, items= grep("^T[[:digit:]]{2}", colnames(datW)),
        id="idstud", analysis.name = "twodim",  qMatrix = qMat,
        HG.var = c("language","sex", "ses"), anchor = aPar, n.plausible = 20,
        software="tam", fac.oldxsi=.1, increment.factor=1.07)

# run the model
run2a<- runModel(mod2a)

# get the results
res2a<- getResults(run2a)


################################################################################
###        Example 2b: Multidimensional Rasch Model with equating            ###
################################################################################

# Example 2b: running a multidimensional Rasch model on a subset of items
# defining the model: specifying q matrix now is necessary.
mod2b<- defineModel(dat=datW, items= grep("^T[[:digit:]]{2}", colnames(datW)),
        id="idstud", analysis.name = "twodim2", qMatrix = qMat,
        n.plausible = 20, software="tam", fac.oldxsi=.1, increment.factor=1.07)

# run the model
run2b<- runModel(mod2b)

# get the results
res2b<- getResults(run2b)

### equating (wenn nicht verankert)
eq2b <- equat1pl( results = res2b, prmNorm = aPar)

### transformation to the 'bista' metric: needs reference population definition
ref  <- data.frame ( domain = c("domainreading", "domainlistening"),
        m = c(0.03890191, 0.03587727), sd= c(1.219, 0.8978912))
cuts <- list ( domainreading = list ( values = 390+0:3*75),
        domainlistening = list ( values = 360+0:3*85))
tf2b <- transformToBista ( equatingList = eq2b, refPop = ref, cuts = cuts)


################################################################################
###            Example 3: Multidimensional Rasch Model in TAM                ###
################################################################################

# Example 3: the same model in TAM
# we use the same anchor parameters from example 1

# estimate model 2 with latent regression and anchored parameters in TAM
# specification of an output folder (via 'dir' argument) no longer necessary
mod2T<- defineModel(dat=datW, items= grep("^T[[:digit:]]{2}", colnames(datW)),
        id="idstud", qMatrix = qMat, HG.var = "sex", anchor = aPar, software = "tam")

# run the model
run2T<- runModel(mod2T)

# Object 'run2T' is of class 'tam.mml'
class(run2T)

# the class of 'run2T' corresponds to the class defined by the TAM package; all
# functions of the TAM package intended for further processing (e.g. drawing
# plausible values, plotting deviance change etc.) work, for example:
wle  <- tam.wle(run2T)

# Finally, the model result are collected in a single data frame
res2T<- getResults(run2T)


################################################################################
###    Example 4: define und run multiple models defined by 'splitModels'    ###
################################################################################

# Example 4: define und run multiple models defined by 'splitModels'
# Model split is possible for different groups of items (i.e. domains) and/or
# different groups of persons (for example, federal states within Germany)

# define person grouping
pers  <- data.frame ( idstud = datW[,"idstud"] , group1 = datW[,"sex"],
         group2 = datW[,"language"], stringsAsFactors = FALSE )

# define 18 models, splitting according to person groups and item groups separately
# by default, multicore processing is applied
l1    <- splitModels ( qMatrix = qMat, person.groups = pers, nCores = 1)

# apply 'defineModel' for each of the 18 models in 'l1'
modMul<- defineModel(dat = datW, items = grep("^T[[:digit:]]{2}", colnames(datW)),
         id = "idstud", check.for.linking = TRUE, splittedModels = l1, software = "tam")

# run all models
runMul<- runModel(modMul)

# get results of all models
resMul<- getResults(runMul)


################################################################################
###          Example 5: Linking and equating for multiple models             ###
################################################################################

# Example 5: define und run multiple models according to different domains (item groups)
# and further linking/equating. This example mimics the routines necessary for the
# 'Vergleichsarbeiten' at the Institute of Educational Progress (IQB)

# specify two models according to the two domains 'reading' and 'listening'
l2    <- splitModels ( qMatrix = qMat,  nCores = 1)

# define 2 models
mods  <- defineModel(dat = datW, id = "idstud", check.for.linking = TRUE,
         splittedModels = l2, software = "tam")

# run 2 models
runs  <- runModel(mods)

# get the results
ress  <- getResults(runs)

# only for illustration, we create arbitrary 'normed' parameters for anchoring
prmNrm<- itemFromRes(ress)[ sample ( 1:56, 31,FALSE) ,c("item", "est")]
prmNrm[,"est"] <- prmNrm[,"est"] - 0.6 + rnorm ( 31, 0, 0.75)

# anchoring without exclusion of linking DIF items (DIF items will only be identified)
anch  <- equat1pl ( results = ress, prmNorm = prmNrm, excludeLinkingDif = FALSE,
         difBound = 0.6)

# anchoring with exclusion of linking DIF items
anch2 <- equat1pl ( results = ress, prmNorm = prmNrm, excludeLinkingDif = TRUE,
         difBound = 0.6, iterativ = FALSE)

# anchoring with iterative exclusion of linking DIF items
anch3 <- equat1pl ( results = ress, prmNorm = prmNrm, excludeLinkingDif = TRUE,
         difBound = 0.6, iterativ = TRUE)

# transformation to the Bista metric
# first we arbitrarily define mean and standard deviation of the reference
# population according to both dimensions (defined in the Q matrix):
# reading and listening
# Note that the the first column of the 'refPop' data frame must include the
# domain names. Domain names must match the names defined in the Q matrix
refPop<- data.frame ( domain = c("domainreading", "domainlistening"),
        m = c(0.03890191, 0.03587727), sd= c(1.219, 0.8978912))

# second, we specify a list with cut scores. Values must be in ascending order.
# Labels of the competence stages are optional. If no labels are specified,
# the will be defaulted to 1, 2, 3 ... etc.
# Note: if labels are specified, there must be one label more than cut scores.
# (i.e. 4 cut scores need 5 labels, etc.)
cuts  <- list ( domainreading = list ( values = 390+0:3*75),
         domainlistening = list ( values = 360+0:3*85))

# transformation
dfr   <- transformToBista ( equatingList = anch3, refPop = refPop, cuts=cuts )
head(dfr$itempars)
head(dfr$personpars)


################################################################################
###      Example 5a: Linking for multiple models, including global domain    ###
################################################################################

# Example 5a: define und run multiple models according to different domains (item groups)
# and further linking/equating. Same as example 5, but extended for the 'global'
# domain.

# add the 'global' domain (reading and listening together) in the Q matrix
qMat2 <- qMat
qMat2[,"global"] <- 1

# specify two models according to the two domains 'reading' and 'listening'
l3    <- splitModels ( qMatrix = qMat2, nCores = 1)

# define 3 models
mods3 <- defineModel(dat = datW, id = "idstud", check.for.linking = TRUE,
         splittedModels = l3, software = "tam")

# run 3 models
runs3 <- runModel(mods3)

# get the results
res3  <- getResults(runs3)

# only for illustration, we create arbitrary 'normed' parameters for anchoring.
# Each item now has to item parameter: one is domain-specific, one is the global
# item parameter. Hence, each item occurs two times in 'prmNrm'
prmNrm<- itemFromRes(ress)[ sample ( 1:56, 31,FALSE) ,c("dimension", "item", "est")]
prmNrm[,"est"] <- prmNrm[,"est"] - 0.6 + rnorm ( 31, 0, 0.75)
prmGlo<- prmNrm
prmGlo[,"dimension"] <- "global"
prmNrm<- rbind ( prmNrm, prmGlo)

# if the item identifier in 'prmNrm' is not unique, 'equat1pl' has to know which
# parameter in 'prmNrm' belongs to which dimension/domain. Hence, the dimension
# is added to 'prmNrm'.

# anchoring: if 'prmNrm' has more than 2 columns, the columns of 'prmNrm' must be
# specified in 'equat1pl'
anch3 <- equat1pl ( results = res3, prmNorm = prmNrm, item = "item", value = "est",
         domain = "dimension", excludeLinkingDif = FALSE, difBound = 0.6)

# transformation to the Bista metric
# first we arbitrarily define mean and standard deviation of the reference
# population according to the three dimensions (defined in the Q matrix):
# reading, listening, and global
# Note that the the first column of the 'refPop' data frame must include the
# domain names. Domain names must match the names defined in the Q matrix
refPop<- data.frame ( domain = c("domainreading", "domainlistening", "global"),
        m = c(0.03890191, 0.03587727, 0.03), sd= c(1.219, 0.8978912, 1.05))

# second, we specify a list with cut scores. Values must be in ascending order.
# Labels of the competence stages are optional. If no labels are specified,
# the will be defaulted to 1, 2, 3 ... etc.
# Note: if labels are specified, there must be one label more than cut scores.
# (i.e. 4 cut scores need 5 labels, etc.)
cuts  <- list ( domainreading = list ( values = 390+0:3*75),
         domainlistening = list ( values = 360+0:3*85),
         global = list ( values = 400+0:2*100, labels = c("A1", "A2", "B1", "B2")))

# transformation
dfr   <- transformToBista ( equatingList = anch3, refPop = refPop, cuts=cuts )
head(dfr$itempars)
head(dfr$personpars)


################################################################################
###    Example 6: Scaling, linking and equating for multiple models (II)     ###
################################################################################

# Example 6 and 6a: define und run multiple models according to different domains
# (item groups) and different person groups. This example mimics the routines
# necessary for the 'Laendervergleich/Bildungstrend' at the Institute for
# Educational Progress (IQB). Example 6 demonstrates routines without trend
# estimation, i.e. for the first (or solely) measurement occasion. Example 6a
# demonstrates routines for the second measurement occasion (t2) which is linked
# to t1. Example 6b demonstrates routines for the third measurement occasion (t3)
# which is linked to the common scale of t1 and t2.

# data set 'trends' contains item response data for three measurement occasions
# in a single data.frame
data(trends)

# Preparation: assume time of measurement 't1' corresponds to the year 2010.
# This is the year of the reference population. We consider both domains,
# reading and listening in a single call. Hence, the data.frame 'datT1' contains
# items of both domains
datT1<- reshape2::dcast(subset ( trends, year == 2010),
        idstud+country+sex+ses+language~item, value.var="value")

# generate Q matrix
qMat <- unique(trends[ ,c("item","domain")])
qMat <- data.frame ( qMat[,"item", drop=FALSE], model.matrix(~domain-1, data = qMat))

# First step: item calibration in separate unidimensional models for each domain
# (the models are short and simple, so we don't need multicore)
modsT1<- splitModels ( qMatrix = qMat, nCores = 1)

# define 2 models. Note: not all items of the Q matrix are present in the data.
# Items which occur only in the Q matrix will be ignored.
defT1 <- defineModel(dat = datT1, id = "idstud", check.for.linking = TRUE,
         splittedModels = modsT1, software = "tam")

# run (calibrate) the 2 models subsequently
runT1 <- runModel(defT1)

# get the results of the two unidimensional models
resT1 <- getResults(runT1, omitPV = TRUE, omitWle = TRUE, Q3 = FALSE)

# extract item parameters from the 'results' object
# t1 is the reference measurement occasion, i.e. no linking/equating is necessary
itemT1<- itemFromRes(resT1)

# Second step: drawing plausible values separately for each country.
# A two-dimensional (reading/listening) model is specified separately for each
# person group (= each country) with item parameters fixed at their calibration
# values. Moreover, a latent regression model is used (in the actual 'Laendervergleich',
# regressors are principal components of previously imputed background variables). We use 'sex',
# 'ses' and 'language' as regressors. For convenience, 'ses' is scaled (mean = 0, sd = 1)
datT1[,"ses_scaled"] <- scale(datT1[,"ses"])[,1]

# have a look at 'sex' and 'language at home' in each country:
table(datT1[,c("country", "sex")])
table(datT1[,c("country", "language")])

# Running second step: split models according to person groups
# ('all.persons' must be FALSE, otherwise the whole group would be treated as
# a separate distinct group.)
modT1P<- splitModels ( person.groups = datT1[,c("idstud", "country")],
         all.persons = FALSE, nCores = 1)

# define the 2 country-specific 2-dimensional models, specifying latent regression
# model and fixed item parameters.
defT1P<- defineModel(dat = datT1, items = itemT1[,"item"], id = "idstud",
         check.for.linking = TRUE, splittedModels = modT1P, qMatrix = qMat,
         anchor = itemT1[,c("item", "est")],
         HG.var = c("sex", "ses_scaled", "language"),  software = "tam")

# run the 2 models (estimation needs approx. 20 seconds)
runT1P<- runModel(defT1P)

# get the results (to save time, item fit estimation is skipped)
resT1P<- getResults(runT1P, omitWle = TRUE, Q3 = FALSE)

# latent regression coefficients for the three countries and two dimensions
regcoefFromRes(resT1P, digits = 3)

# equating is not necessary, as the models run with fixed item parameters
# However, to prepare for the transformation on the 'bista' metric, run
# 'equat1pl' with empty arguments
ankT1P<- equat1pl ( results = resT1P)

# transformation to the 'bista' metric
# Note: if the sample was drawn from the reference population, mean and SD
# are not yet known. So we ignore the 'refPop' argument in 'transformToBista'
# and simply define the cut scores. The function then assumes that the parameter
# stem from the reference population and estimates its mean and sd.
cuts  <- list ( domainreading = list ( values = 390+0:3*75),
         domainlistening = list ( values = 360+0:3*85))

# transformation: omit 'refPop' argument
dfrT1P<- transformToBista ( equatingList = ankT1P, cuts=cuts, vera=FALSE )

# mean and sd of the reference population (population of 't1')
dfrT1P[["refPop"]]


################################################################################
###     Example 6a: Extend example 6 with trend estimation (now for t2)      ###
################################################################################

# Example 6a needs the objects (Q matrix, item parameters, ...) created in example 6
# Preparation: assume time of measurement 't2'.
datT2<- reshape2::dcast(subset ( trends, year == 2015),
        idstud+country+sex+ses+language~item, value.var="value")

# First step: item calibration in separate unidimensional models for each domain
modsT2<- splitModels ( qMatrix = qMat, nCores = 1)

# define 2 models. Items which occur only in the Q matrix but not in the data
# will be ignored.
defT2 <- defineModel(dat = datT2, id = "idstud", check.for.linking = TRUE,
         splittedModels = modsT2, software = "tam")

# run 2 models for calibration
runT2 <- runModel(defT2)

# get the results
resT2 <- getResults(runT2)

# collect item parameters
itemT2<- itemFromRes(resT2)

# Second step: compute linking constant between 't1' and 't2' with the iterative
# exclusion of linking DIF items and computation of linking error. We use the
# 'itemT1' object created in example 6 for reference item parameters. The linking
# procedure is executed consecutively for listening and reading.
L.t1t2<- equat1pl ( results = resT2, prmNorm = itemT1[,c("item", "est")],
         excludeLinkingDif = TRUE, difBound = 0.64, iterativ = TRUE)

# linking constant is negative: students performance at T2 is worse than T1
# see that DIF is non symmetrical for reading: the four items with the highest
# amount of DIF have all positive DIF values. DIF exclusion hence shrinks linking
# constant towards zero.
# Third step: transform item parameters of 't2' to the metric of 't1'
# We now need to specify the 'refPop' argument. We use the values from 't1' which
# serves as the reference. To capture linking errors in a separate data.frame
# within the returned list, we define the years of assessment
ref   <- dfrT1P[["refPop"]]
T.t1t2<- transformToBista ( equatingList = L.t1t2, refPop=ref, cuts = cuts,
         vera=FALSE, years = c(2010,2015))

# The object 'T.t1t2' now contains transformed person and item parameters with
# original and transformed linking errors. See for example person parameter:
head(T.t1t2$personpars)

# Fourth step: drawing plausible values for 't2'. We use the transformed item
# parameters (captured in 'T.t1t2') for anchoring

# Running second step: split models according to person groups (countries)
# ('all.persons' must be FALSE, otherwise the whole group would be treated as
# a separate distinct group.)
modT2P<- splitModels ( person.groups = datT2[,c("idstud", "country")] ,
         all.persons = FALSE, nCores = 1)

# define the 2 country-specific 2-dimensional models, specifying latent regression
# model and fixed item parameters. We used the transformed item parameters (captured
# in 'T.t1t2[["itempars"]]' --- using the 'estTransf' column) for anchoring.
# Again, 'ses' is scaled (mean = 0, sd = 1)
datT2[,"ses_scaled"] <- scale(datT2[,"ses"])[,1]
defT2P<- defineModel(dat = datT2, items = itemT2[,"item"], id = "idstud",
         check.for.linking = TRUE, splittedModels = modT2P, qMatrix = qMat,
         anchor = T.t1t2[["itempars"]][,c("item", "estTransf")],
         HG.var = c("sex", "ses_scaled", "language"), software = "tam")

# run the 2 models (estimation takes approx. 29 seconds)
runT2P<- runModel(defT2P)

# get the results
resT2P<- getResults(runT2P)

# equating is not necessary, as the models run with fixed item parameters
# However, to prepare for the transformation on the 'bista' metric, run
# 'equat1pl' with empty arguments
ankT2P<- equat1pl ( results = resT2P)

# transformation to the 'bista' metric, using the previously defined cut scores
# and the reference population mean and sd from 't1'
dfrT2P<- transformToBista ( equatingList = ankT2P, refPop=ref, cuts=cuts, vera=FALSE)


################################################################################
###     Example 6b: Extend example 6a with trend estimation (now for t3)     ###
################################################################################

# Example 6b needs the objects (Q matrix, item parameters, ...) created in example 6
# and 6a. Preparation: assume time of measurement 't3'.
datT3<- reshape2::dcast(subset ( trends, year == 2020),
        idstud+country+sex+ses+language~item, value.var="value")

# First step: item calibration in separate unidimensional models for each domain
modsT3<- splitModels ( qMatrix = qMat, nCores = 1)

# define 2 models. Items which occur only in the Q matrix but not in the data
# will be ignored.
defT3 <- defineModel(dat = datT3, id = "idstud", check.for.linking = TRUE,
         splittedModels = modsT3, software = "tam")

# run 2 models
runT3 <- runModel(defT3)

# get the results
resT3 <- getResults(runT3)

# collect item parameters
itemT3<- itemFromRes(resT3)

# Second step: compute linking constant. We link the items of 't3' to the items
# of 't2' which are already transformed to the metric of 't1' (we use the item
# parameters which were used for plausible values imputation at 't2' as norm
# parameters)
L.t2t3<- equat1pl ( results = resT3, prmNorm = T.t1t2[["itempars"]][,c("item", "estTransf")],
         excludeLinkingDif = TRUE, difBound = 0.64, iterativ = TRUE)

# linking constant is negative: students performance at T3 is worse than T1
# Third step: transform item parameters of 't3' to the common metric of 't1' and 't2'
# We already know the 'refPop' values.
ref   <- dfrT1P[["refPop"]]
T.t2t3<- transformToBista ( equatingList = L.t2t3, refPop=ref, cuts = cuts,
         vera=FALSE, years = c(2015,2020))

# Fourth step: drawing plausible values for 't3'. We use the transformed item
# parameters (captured in 'T.t2t3') for anchoring

# Running second step: split models according to person groups (countries)
# ('all.persons' must be FALSE, otherwise the whole group would be treated as
# a separate distinct group.)
modT3P<- splitModels ( person.groups = datT3[,c("idstud", "country")] ,
         all.persons = FALSE, nCores = 1)

# define the 2 country-specific 2-dimensional models, specifying latent regression
# model and fixed item parameters. We used the transformed item parameters (captured
# in 'T.t2t3[["itempars"]]' --- using the 'estTransf' column) for anchoring.
# Again, 'ses' is scaled (mean = 0, sd = 1)
datT3[,"ses_scaled"] <- scale(datT3[,"ses"])[,1]
defT3P<- defineModel(dat = datT3, items = itemT3[,"item"], id = "idstud",
         check.for.linking = TRUE, splittedModels = modT3P, qMatrix = qMat,
         anchor = T.t2t3[["itempars"]][,c("item", "estTransf")],
         HG.var = c("sex", "ses_scaled", "language"), software = "tam")

# run the 2 models (estimation takes approx. 20 seconds)
runT3P<- runModel(defT3P)

# get the results
resT3P<- getResults(runT3P)

# equating is not necessary, as the models run with fixed item parameters
# However, to prepare for the transformation on the 'bista' metric, run
# 'equat1pl' with empty arguments
ankT3P<- equat1pl ( results = resT3P)

# transformation to the 'bista' metric, using the previously defined cut scores
# and the reference population mean and sd from 't1'
dfrT3P<- transformToBista ( equatingList = ankT3P, refPop=ref, cuts=cuts, vera=FALSE)


################################################################################
###                 Example 6c: Prepare trend estimation                     ###
################################################################################

# Example 6c needs the objects (Q matrix, item parameters, ...) created in example 6,
# 6a, and 6b.
# We now collect the person parameter estimates (plausible values) from t1, t2, and t3
# in a common data.frame. The person estimates are already collected in the objects
# previously created by 'transformToBista()'. Not all columns are necessary.

# plausibles values of measurement occasion 1 ('t1'): add the year to the data.frame
persT1<- data.frame ( year = 2010,
         dfrT1P[["personpars"]][,c("idstud", "dimension", "imp", "value", "valueTransfBista", "traitLevel")],
         stringsAsFactors = FALSE)

# plausibles values of measurement occasion 2 ('t2'): add the year to the data.frame
persT2<- data.frame ( year = 2015,
         dfrT2P[["personpars"]][,c("idstud", "dimension", "imp", "value", "valueTransfBista", "traitLevel")],
         stringsAsFactors = FALSE)

# plausibles values of measurement occasion 3 ('t3'): add the year to the data.frame
persT3<- data.frame ( year = 2020,
         dfrT3P[["personpars"]][,c("idstud", "dimension", "imp", "value", "valueTransfBista", "traitLevel")],
         stringsAsFactors = FALSE)

# bind together in a common data.frame
pers  <- rbind(persT1, persT2, persT3)

# merge background variables to plausible values data
# first we have to create the 'domain' column in plausible values data
pers[,"domain"] <- car::recode(pers[,"dimension"], "'domainlistening'='listening'; 'domainreading'='reading'")
pers[,"dimension"] <- NULL
pers  <- merge(unique(trends[,c("year", "idclass", "idstud", "domain", "country", "language", "ses", "sex")]),
         pers, by = c("year", "idstud", "domain"), all = FALSE)

# collect linking errors
# t1 vs. t2: linking errors were computed in example 6a.
let1t2<- T.t1t2[["linkingErrors"]]

# t2 vs. t3: linking errors were computed in example 6b.
let2t3<- T.t2t3[["linkingErrors"]]

# t1 vs. t3: linking errors were not yet computed: link t3 to t1 to create linking error template
L.t1t3<- equat1pl ( results = resT3, prmNorm = itemT1[,c("item", "est")],
         excludeLinkingDif = TRUE, difBound = 0.64, iterativ = TRUE)

# indirect linking ('chained' linking)
chain <- multiEquatError (x1=resT1, x2=resT2, x3=resT3, difBound = 0.64, verbose = TRUE )

# replace direct linking errors with indirect linking errors
L.t1t3<- replaceLinkingError (equatingList =L.t1t3, multiEquatError_output=chain)

# transform linking errors
ref   <- dfrT1P[["refPop"]]
tle   <- transformToBista ( equatingList = L.t1t3, refPop=ref, cuts = cuts,
         vera=FALSE, years = c(2010,2020))
let1t3<- tle[["linkingErrors"]]

# bind all linking errors in a common data.frame
lErr  <- rbind(let1t2, let2t3, let1t3)
lErr[,"domain"] <- car::recode(lErr[,"domain"], "'domainlistening'='listening'; 'domainreading'='reading'")


################################################################################
###  Example 6d: eatRep estimation with plausible values and linking errors  ###
################################################################################

# Example 6d needs the objects ('pers', 'lErr') created in example 6c
# load the 'eatRep' package ... note: needs eatRep version 0.14.0 or higher
library(eatRep)

# compute means for both countries with trend, for both domains separately,
# using replications methods (jackknife-1)
means <- by(data = pers, INDICES = pers[,"domain"], FUN = function ( dim ) {
         m <- repMean(datL = dim, ID="idstud", PSU = "idclass", type = "jk1",
              imp = "imp", groups = "country", dependent = "valueTransfBista",
              trend = "year", linkErr = lErr[which(lErr[,"domain"] == dim[1,"domain"]),])
         r <- report(m, add = list(domain = dim[1,"domain"]))
         return(r)})
means <- do.call("rbind", means)

# additionally: differ the sex-specific means in each country from the sex-specific means
# in the whole population? Are the differences (male vs. female) in each country different
# from the difference (male vs. female) in the whole population?
means2<- by(data = pers, INDICES = pers[,"domain"], FUN = function ( dim ) {
         m <- repMean(datL = dim, ID="idstud", PSU = "idclass", type = "jk1",
              imp = "imp", groups = c("country","sex"), group.differences.by = "sex",
              group.splits = 0:1, cross.differences = TRUE,crossDiffSE.engine= "lm",
              dependent = "valueTransfBista", trend = "year",
              linkErr = lErr[which(lErr[,"domain"] == dim[1,"domain"]),])
         r <- report(m, add = list(domain = dim[1,"domain"]))
         return(r)})
means2<- do.call("rbind", means2)


################################################################################
###  Example 6e: eatRep (repTable) with plausible values and linking errors  ###
################################################################################

# Example 6e needs the objects ('pers', 'lErr') created in example 6c
# load the 'eatRep' package ... note: needs eatRep version 0.14.0 or higher
library(eatRep)

# compute frequencies for trait levels, for both domains, with trend
freqs <- by(data = pers, INDICES = pers[,"domain"], FUN = function ( dim ) {
         m <- repTable(datL = dim, ID="idstud", PSU = "idclass", type = "jk1",
              imp = "imp", groups = "country", dependent = "traitLevel",
              trend = "year", linkErr = lErr[which(lErr[,"domain"] == dim[1,"domain"]),])
         r <- report(m, add = list(domain = dim[1,"domain"]))
         return(r)})
freqs <- do.call("rbind", freqs)

# additionally: sex differences in each country, using 'group.differences.by' argument
# Note: for frequency tables group differences may result in a chi square test or in
# a difference of each categories' frequency.
# first: request chi square test
freqs1<- by(data = pers, INDICES = pers[,"domain"], FUN = function ( dim ) {
         m <- repTable(datL = dim, ID="idstud", PSU = "idclass", type = "jk1",
              imp = "imp", groups = c("country","sex"), group.differences.by = "sex",
              chiSquare = TRUE,dependent = "traitLevel", trend = "year",
              linkErr = lErr[which(lErr[,"domain"] == dim[1,"domain"]),])
         r <- report(m, add = list(domain = dim[1,"domain"]))
         return(r)})
freqs1<- do.call("rbind", freqs1)

# differences for each competence level (chiSquare = FALSE)
# (for faster computation, we omit jackknife procedure)
freqs2<- by(data = pers, INDICES = pers[,"domain"], FUN = function ( dim ) {
         m <- repTable(datL = dim, ID="idstud", type = "none", imp = "imp",
              groups = c("country","sex"), group.differences.by = "sex",
              chiSquare = FALSE,dependent = "traitLevel", trend = "year",
              linkErr = lErr[which(lErr[,"domain"] == dim[1,"domain"]),],
              group.splits = 0:2, cross.differences = TRUE)
         r <- report(m, add = list(domain = dim[1,"domain"]))
         return(r)})
freqs2<- do.call("rbind", freqs2)


################################################################################
###        Example 7: Linking and equating for multiple models (III)         ###
################################################################################

# For the purpose of illustration, this example repeats analyses of example 6,
# using a 2pl estimation now. Example 7 demonstrates routines without trend
# estimation.

# Preparation: assume time of measurement 't1' corresponds to the year 2003.
datT1<- reshape2::dcast(subset ( trends, year == 2010),
        formula = idstud+sex+country+language+ses~item, value.var="value")

# First step: item calibration in separate unidimensional models for each domain
# split 2 models. Note: not all items of the Q matrix are present in the data.
# Items which occur only in the Q matrix will be ignored.
modsT1<- splitModels ( qMatrix = qMat, nCores = 1)

# lets specify a 2pl model with constraints: a common discrimination for all
# items belonging to the same task
slopes<- data.frame ( variable = qMat[,"item"],
         slope = as.numeric(as.factor(substr(qMat[,"item"],1,3))))

# prepare 2pl model
defT1 <- defineModel(dat = datT1, id = "idstud", check.for.linking = TRUE,
         splittedModels = modsT1, irtmodel = "2PL.groups", est.slopegroups = slopes,
         software = "tam")

# run 2 models
runT1 <- runModel(defT1)

# get the results
resT1 <- getResults(runT1)

# extract item parameters from the 'results' object
itemT1<- itemFromRes(resT1)

# Second step: drawing plausible values. Two-dimensional model is specified for
# each person group with fixed item parameters. Moreover, a latent regression
# model is used (in the actual 'Laendervergleich', regressors are principal
# components).
# define person grouping
# Running second step: split models according to person groups
# ('all.persons' must be FALSE, otherwise the whole group would be treated as
# a separate distinct group.)
modT1P<- splitModels ( person.groups = datT1[,c("idstud", "country")], all.persons = FALSE, nCores = 1)

# define the 2 country-specific 2-dimensional models, specifying latent regression
# model and fixed item and fixed slope parameters.
defT1P<- defineModel(dat = datT1, items = itemT1[,"item"], id = "idstud", irtmodel = "2PL",
         check.for.linking = TRUE, splittedModels = modT1P, qMatrix = qMat,
         anchor = itemT1[,c("item", "est")], fixSlopeMat = itemT1[,c("item", "estSlope")],
         HG.var = c("ses", "sex", "language"), software = "tam")

# run the 2 models
runT1P<- runModel(defT1P)

# get the results
resT1P<- getResults(runT1P, Q3 = FALSE)
}

