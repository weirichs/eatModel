# Prepares IRT analysis for Conquest, TAM, or mirt

Facilitates data analysis using the software Conquest, TAM, or mirt. It
automatically checks data for IRT consistency, generates Conquest
syntax, label, anchor and data files or corresponding TAM/mirt call for
a single model specified by several arguments in R. Finally, an R object
is created which contain the required input for Conquest or TAM. To
start the estimation, call
[`runModel`](https://weirichs.github.io/eatModel/reference/runModel.md)
with the argument returned by `defineModel`.

## Usage

``` r
defineModel (dat, items, id, splittedModels = NULL,
   irtmodel = c("1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", "GPCM.groups", "2PL.groups", "GPCM.design", "3PL"),
   qMatrix=NULL, DIF.var=NULL, HG.var=NULL, group.var=NULL, weight.var=NULL, anchor = NULL, 
   domainCol=NULL, itemCol=NULL, valueCol=NULL, catCol = NULL, check.for.linking = TRUE, minNperItem = 50, removeMinNperItem = FALSE,
   boundary = 6, remove.boundary = FALSE, remove.no.answers = TRUE, remove.no.answersHG = TRUE, 
   remove.missing.items = TRUE, remove.constant.items = TRUE, remove.failures = FALSE, 
   remove.vars.DIF.missing = TRUE, remove.vars.DIF.constant = TRUE, 
   verbose=TRUE, software = c("conquest","tam", "mirt"), dir = NULL, analysis.name,
   schooltype.var = NULL, model.statement = "item",  compute.fit = TRUE,
   pvMethod = c("regular", "bayesian"), fitTamMmlForBayesian = TRUE,
   n.plausible=5, seed = NULL, conquest.folder= NULL,
   constraints=c("cases","none","items"), std.err=c("quick","full","none"), distribution=c("normal","discrete"),
   method=c("gauss", "quadrature", "montecarlo", "quasiMontecarlo"), n.iterations=2000,
   nodes=NULL, p.nodes=2000, f.nodes=2000,converge=0.001,deviancechange=0.0001,
   equivalence.table=c("wle","mle","NULL"), use.letters=FALSE,
   allowAllScoresEverywhere = TRUE, guessMat = NULL, est.slopegroups = NULL,
   fixSlopeMat = NULL, slopeMatDomainCol=NULL, slopeMatItemCol=NULL, slopeMatValueCol=NULL, 
   progress = NULL, Msteps = NULL, increment.factor=1 , fac.oldxsi=0,
   export = list(logfile = TRUE, systemfile = FALSE, history = TRUE,
   covariance = TRUE, reg_coefficients = TRUE, designmatrix = FALSE))
```

## Arguments

- dat:

  A data frame containing all variables necessary for analysis.

- items:

  Names or column numbers of variables with item responses. Item
  response values must be numeric (i.e. 0, 1, 2, 3 ... ). Character
  values (i.e. A, B, C ... or a, b, c, ...) are not allowed. Class of
  item columns are expected to be numeric or integer. Columns of class
  `character` will be transformed.

- id:

  Name or column number of the identifier (ID) variable.

- splittedModels:

  Optional: Object returned by
  [`splitModels`](https://weirichs.github.io/eatModel/reference/splitModels.md).
  Definition for multiple model handling.

- irtmodel:

  If `software = "conquest"`, the argument is ignored as the IRT model
  is specified with the `model.statement` argument. If
  `software = "tam"`, the argument `irtmodel` corresponds to the
  `irtmodel` argument of
  [`tam.mml`](https://rdrr.io/pkg/TAM/man/tam.mml.html). See the help
  page of [`tam.mml`](https://rdrr.io/pkg/TAM/man/tam.mml.html) for
  further details. If `software = "mirt"`, `irtmodel` is a data.frame
  with two columns. First column: item identifier. Second column: Type
  of response model for this item. See the `itemtype` argument in the
  help file of
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.html) for
  further details. Additionally, see example 10 below for specifying a
  mixed model with Rasch-, 2pl, partial credit, and generalized partial
  credit items.

- qMatrix:

  Optional: A named data frame indicating how items should be grouped to
  dimensions. The first column contains the unique names of all items
  and should be named “item”. The other columns contain dimension
  definitions and should be named with the respective dimension names. A
  positive value (e.g., 1 or 2 or 1.4) indicates the loading weight with
  which an item loads on the dimension, a value of 0 indicates that the
  respective item does not load on this dimension. If no q matrix is
  specified by the user, an unidimensional structure is assumed.

- DIF.var:

  Name or column number of one grouping variable for which differential
  item functioning analysis is to be done.

- HG.var:

  Optional: Names or column numbers of one or more context variables
  (e.g., sex, school). These variables will be used for latent
  regression model in Conquest or TAM.

- group.var:

  Applies only if `software = "tam"`. Optional: Optional: Names or
  column numbers of one or more grouping variables. Descriptive
  statistics for WLEs and Plausible Values will be computed separately
  for each group in Conquest.

- weight.var:

  Optional: Name or column number of one weighting variable.

- anchor:

  Optional: A named data frame with anchor parameters. If the data.frame
  only has two columns, the columns do not need to be named. In this
  case, the item name should be in the first column. The second column
  contains anchor parameters. Anchor items can be a subset of the items
  in the dataset and vice versa. If the data frame contains more than
  two columns, columns must be named explicitly using the following
  arguments `domainCol`, `itemCol`, `valueCol`, and `catCol`.

- domainCol:

  Optional: Only necessary if the `anchor` argument was used to define
  anchor parameters. Moreover, specifying `domainCol` is only necessary,
  if the item identifiers in `anchor` are not unique—for example, if a
  specific item occurs with two parameters, one domain-specific item
  parameter and one additional “global” item parameter. The domain
  column than must specify which parameter belongs to which domain.

- itemCol:

  Optional: Only necessary if the `anchor` argument was used to define
  anchor parameters. Moreover, specifying `itemCol` is only necessary,
  if the `anchor` data frame has more than two columns. The `itemCol`
  column then must specify which column contains the item identifier.

- valueCol:

  Optional: Only necessary if the `anchor` argument was used to define
  anchor parameters. Moreover, specifying `valueCol` is only necessary,
  if the `anchor` data frame has more than two columns. The `valueCol`
  column then must specify which column contains the item parameter
  values.

- catCol:

  Optional: Only necessary if the `anchor` argument was used to define
  anchor parameters in a partial credit model. The `catCol` column then
  must specify which column contains the item parameter category
  indices. The values in this column must be named according to the
  convention `"Cat1"`, `"Cat2"`, `"Cat3"`, and so on.

- check.for.linking:

  A logical value indicating whether the items in dataset are checked
  for being connected with each other via design.

- minNperItem:

  Numerical: A message is printed on console if an item has less valid
  values than the number defined in `minNperItem`.

- removeMinNperItem:

  Logical: Remove items with less valid responses than defined in
  `minNperItem`?

- boundary:

  Numerical: A message is printed on console if a subject has answered
  less than the number of items defined in boundary.

- remove.boundary:

  Logical: Remove subjects who have answered less items than defined in
  the `boundary` argument?

- remove.no.answers:

  Logical: Should persons without any item responses being removed prior
  to analysis?

- remove.no.answersHG:

  Logical: Should persons without any responses on any background
  variable being removed prior to analysis?

- remove.missing.items:

  Logical: Should items without any item responses being removed prior
  to analysis?

- remove.constant.items:

  Logical: Should items without variance being removed prior to
  analysis?

- remove.failures:

  Logical: Should persons without any correct item response (i.e., only
  “0” responses) being removed prior to analysis?

- remove.vars.DIF.missing:

  Logical: Applies only in DIF analyses. Should items without any
  responses in at least one DIF group being removed prior to analyses?
  Note: Conquest may crash if these items remain in the data.

- remove.vars.DIF.constant:

  Logical: Applies only in DIF analyses. Should items without variance
  in at least one DIF group being removed prior to analyses? Note:
  Conquest may crash if these items remain in the data.

- verbose:

  A logical value indicating whether messages are printed on the R
  console.

- software:

  The desired estimation software for the analysis.

- dir:

  The directory in which the output will be written to. If
  `software = "conquest"`, `dir` must be specified. If
  `software = "tam"`, `dir` is not mandatory.

- analysis.name:

  A character string specifying the analysis name. If
  `software = "conquest"`, `analysis.name` must be specified. All
  Conquest input and output files will named `analysis.name` with their
  corresponding extensions. If `software = "tam"`, `analysis.name` is
  not mandatory. In the case of multiple models estimation,
  `split.models` automatically defines `analysis.name` for each model.

- schooltype.var:

  Optional: Name or column number of the variable indicating the school
  type (e.g. academic track, non-academic track). Only necessary if *p*
  values should be computed for each school type separately.

- model.statement:

  Optional: Applies only if `software = "conquest"`. A character string
  given the model statement in the Conquest syntax. If omitted, the
  statement is generated automatically with respect to the defined
  model.

- compute.fit:

  Applies only if `software = "conquest"`. Compute item fit statistics?

- pvMethod:

  Applies only if `software = "tam"`: Specifies whether PVs should be
  drawn regularly or using a Bayesian algorithm.

- fitTamMmlForBayesian:

  Logical, applies only if `software = "tam"`: If PVs are drawn using a
  Bayesian algorithm, it is not necessary to fit the model via
  [`tam.mml`](https://rdrr.io/pkg/TAM/man/tam.mml.html) before.
  `fitTamMmlForBayesian` specifies whether the model should be fitted
  before though. See the help page of `tam.pv.mcmc` for further details.

- n.plausible:

  The number of plausible values which are to be drawn from the
  conditioning model.

- seed:

  Optional: Set seed value for analysis. The value will be used in
  Conquest syntax file ('set seed'-statement, see conquest
  manual, p. 225) or in TAM (control\$seed). Note that seed only occurs
  for stochastic integration.

koennen

- conquest.folder:

  Applies only if `software = "conquest"`. A character string with path
  and name of the Conquest console, for example
  `"c:/programme/conquest/console_Feb2007.exe"`. In package version
  0.7.24 and later, conquest executable file is included in the package,
  so the user needn't to specify this argument.

- constraints:

  A character string specifying how the scale should be constrained.
  Possible options are `"cases"` (default), `"items"` and `"none"`. When
  anchor parameter are specified in anchor, constraints will be set to
  `"none"` automatically. In `TAM` the option `"none"` is not allowed.
  (See the help file of
  [`tam.mml`](https://rdrr.io/pkg/TAM/man/tam.mml.html) for further
  details.)

- std.err:

  Applies only if `software = "conquest"`. A character string specifying
  which type of standard error should be estimated. Possible options are
  `"full"`, `"quick"` (default) and `"none"`. See Conquest manual pp.167
  for details on standard error estimation.

- distribution:

  Applies only if `software = "conquest"`. A character string indicating
  the a priori trait distribution. Possible options are `"normal"`
  (default) and `"discrete"`. See Conquest manual pp.167 for details on
  population distributions.

- method:

  A character string indicating which method should be used for
  numerical or stochastic integration. Possible options are `"gauss"`
  (Gauss-Hermite quadrature: default), `"quadrature"` (Bock/Aitken
  quadrature) and `"montecarlo"`. See Conquest manual pp.167 for details
  on these methods. When using `software = "tam"`, `"gauss"` and
  `"quadrature"` essentially leads to numerical integration, i.e TAM is
  called with `control$snodes = 0` and with
  `control$nodes = seq(-6,6,len=nn)`, where `nn` equals the number of
  nodes specified in the `nodes` argument of `defineModel` (see below).
  Whether or not `control$QMC` is set to TRUE or FALSE, depends on
  whether `"montecarlo"` or `"quasiMontecarlo"` is chosen.
  `"montecarlo"` leads to calling TAM with `control$QMC = FALSE` and
  `snodes = nn`, where `nn` equals the number of nodes specified in the
  `nodes` argument of `defineModel`. `"quasiMontecarlo"` leads to
  calling TAM with `control$QMC = TRUE` and `snodes = nn`, where `nn`
  equals the number of nodes specified in the `nodes` argument of
  `defineModel`. To met the `software = "tam"` default (Quasi Monte
  Carlo integration with `control$QMC = TRUE`), use
  `software="tam", nodes = 21, method = "quasiMontecarlo"`.

- n.iterations:

  An integer value specifying the maximum number of iterations for which
  estimation will proceed without improvement in the deviance.

- nodes:

  An integer value specifying the number of nodes to be used in the
  analysis. The default value is 20. When using `software = "tam"`, the
  value specified here leads to calling TAM with `nodes = 20` AND
  `snodes = 0` if `"gauss"` or `"quadrature"` or `"quasiMontecarlo"` was
  used in the `method` argument. If `"montecarlo"` was used in the
  `method` argument, the value specified here leads to calling TAM with
  `control$snodes = 20`. For numerical integration, for example,
  `method = "gauss"` and `nodes = 21` (TAM default) may be appropriate.
  For quasi monte carlo integration, `method = "quasiMontecarlo"` and
  `nodes = 1000` may be appropriate (TAM authors recommend to use at
  least 1000 nodes).

- p.nodes:

  Applies only if `software = "conquest"`. An integer value specifying
  the number of nodes that are used in the approximation of the
  posterior distributions, which are used in the drawing of plausible
  values and in the calculation of EAP estimates. The default value is
  2000.

- f.nodes:

  Applies only if `software = "conquest"`. An integer value specifying
  the number of nodes that are used in the approximation of the
  posterior distributions in the calculation of fit statistics. The
  default value is 2000.

- converge:

  An integer value specifying the convergence criterion for parameter
  estimates. The estimation will terminate when the largest change in
  any parameter estimate between successive iterations of the EM
  algorithm is less than converge. The default value is 0.001.

- deviancechange:

  An integer value specifiying the convergence criterion for the
  deviance. The estimation will terminate when the change in the
  deviance between successive iterations of the EM algorithm is less
  than deviancechange. The default value is 0.0001.

- equivalence.table:

  Applies only if `software = "conquest"`. A character string specifying
  the type of equivalence table to print. To be more precise, the person
  estimator that is to be used to create the table must be specified
  here. Possible options are `"wle"` (default), `"mle"` and NULL. If
  NULL, no equivalence table is computed at all.

- use.letters:

  Applies only if `software = "conquest"`. A logical value indicating
  whether item response values should be coded as letters. This option
  can be used in partial credit models comprising items with more than
  10 categories to avoid response columns with width 2 in Conquest.

- allowAllScoresEverywhere:

  Applies only if `software = "Conquest"`. Defines score statement
  generation in multidimensional polytomous models. Consider two
  dimensions, \`reading' and \`listening'. In \`reading', values 0, 1,
  2, 3 occur. In \`listening', values 1, 2, 3, 4 occur. If `TRUE`,
  values 0, 1, 2, 3, 4 are defined for both dimensions. Otherwise,
  values 0, 1, 2, 3 are defined for \`reading', values 1, 2, 3, 4 are
  defined for \`listening'.

- guessMat:

  Applies only if `software = "tam"` for 3PL models. A named data frame
  with two columns indicating for which items a common guessing
  parameter should be estimated. The first column contains the names of
  all items in the analysis and should be named `"item"`. The second
  column is numerical (integer values recommended) and allocates the
  items to groups. For each group of items, a separate guessing
  parameter is estimated. If the value in the second columns equals
  zero, the guessing parameter is fixed to zero.

- est.slopegroups:

  Applies only if `software = "tam"` for 2PL models. Optionally, a named
  data frame with two columns indicating for which items a common
  discrimination parameter should be estimated. The first column
  contains the names of all items in the analysis and should be named
  `"item"`. The second column is numerical (integer values recommended)
  and allocates the items to groups. For each group of items, a separate
  discrimination parameter is estimated. Without specifying
  `est.slopegroups`, a discrimination parameter for each item is
  estimated.

- fixSlopeMat:

  Applies only if `software = "tam"` for 2PL models. Optionally, a named
  data frame with two columns indicating for which items a fixed
  discrimination should be assumed. The first column contains the names
  of the items which discrimination should be fixed. Note that item
  indicators should be unique—if not, use further arguments
  `slopeMatDomainCol`, `slopeMatItemCol` and `slopeMatValueCol`. The
  second column is numerical and contains the discrimination value.
  Note: To date, this works only for between item dimensionality models.
  Within item dimensionality models must be specified directly in TAM,
  using the `B.fixed` argument of
  [`tam.mml`](https://rdrr.io/pkg/TAM/man/tam.mml.html). Items which
  discrimation should be estimated should not occur in this data frame.

- slopeMatDomainCol:

  Optional: Only necessary if the `fixSlopeMat` argument was used to
  define fixed slope parameters. Moreover, specifying
  `slopeMatDomainCol` is only necessary, if the item identifiers in
  `fixSlopeMat` are not unique—for example, if a specific item occurs
  with two slope parameters, one domain-specific item slope parameter
  and one additional “global” item parameter. The domain column than
  must specify which parameter belongs to which domain.

- slopeMatItemCol:

  Optional: Only necessary if the `fixSlopeMat` argument was used to
  define fixed slope parameters. Moreover, specifying `itemCol` is only
  necessary, if the `fixSlopeMat` data frame has more than two columns.
  The `itemCol` column then must specify which column contains the item
  identifier.

- slopeMatValueCol:

  Optional: Only necessary if the `fixSlopeMat` argument was used to
  define slope parameters. Moreover, specifying `valueCol` is only
  necessary, if the `fixSlopeMat` data frame has more than two columns.
  The `valueCol` column than must specify which column contains the item
  parameter values.

- progress:

  Logical: Applies only if `software = "tam"` or `software = "mirt"`.
  Print estimation progress messages on console? If `NULL` (the
  default), progress is printed for `software = "mirt"`, but not for
  `software = "tam"`.

- Msteps:

  Number of M steps for item parameter estimation. A high value of M
  steps could be helpful in cases of non-convergence. The default value
  is 4; the default for 3pl models is set to 10.

- increment.factor:

  Applies only if `software = "tam"`. Should only be varied if the model
  does not converge. See help page of
  [`tam.mml`](https://rdrr.io/pkg/TAM/man/tam.mml.html) for further
  details.

- fac.oldxsi:

  Applies only if `software = "tam"`. Should only be varied if the model
  does not converge. See help page of
  [`tam.mml`](https://rdrr.io/pkg/TAM/man/tam.mml.html) for further
  details.

- export:

  Applies only if `software = "conquest"`. Specifies which additional
  files should be written on hard disk.

## Value

A list which contains information about the desired estimation. The list
is intended for further processing via
[`runModel`](https://weirichs.github.io/eatModel/reference/runModel.md).
Structure of the list varies depending on whether multiple models were
called using
[`splitModels`](https://weirichs.github.io/eatModel/reference/splitModels.md)
or not. If
[`splitModels`](https://weirichs.github.io/eatModel/reference/splitModels.md)
was called, the number of elements in the list equals the number of
models defined via
[`splitModels`](https://weirichs.github.io/eatModel/reference/splitModels.md).
Each element in the list is a list with various elements:

- software:

  Character string of the software which is intended to use for the
  further estimation, i.e. `"conquest"` or `"tam"`

- qMatrix:

  The Q matrix allocating items to dimensions.

- all.Names:

  Named list of all relevant variables of the data set.

- dir:

  Character string of the directory the results are to be saved.

- analysis.name:

  Character string of the analysis' name.

- deskRes:

  Data frame with descriptives (e.g., p values) of the test items.

- discrim:

  Data frame with item discrimination values.

- perNA:

  The person identifiers of examinees which are excluded from the
  analysis due to solely missing values.

- per0:

  The person identifiers of examinees which have solely false responses.
  if `remove.failues` was set to be TRUE, these persons are excluded
  from the data set.

- perA:

  The person identifiers of examinees which have solely correct
  responses.

- perExHG:

  The person identifiers of examinees which are excluded from the
  analysis due to missing values on explicit variables.

- itemsExcluded:

  Character string of items which were excluded, for example due to zero
  variance or solely missing values.

If `software == "conquest"`, the output additionally includes the
following elements:

- input:

  Character string of the path with Conquest input (cqc) file.

- conquest.folder:

  Character string of the path of the conquest executable file.

- model.name:

  Character string of the model name.

If `software == "tam"`, the output additionally includes the following
elements:

- anchor:

  Optional: data frame of anchor parameters (if anchor parameters were
  defined).

- daten:

  The prepared data for TAM analysis.

- irtmodel:

  Character string of the used IRT model.

- est.slopegroups:

  Applies for 2pl modeling. Information about which items share a common
  slope parameter.

- guessMat:

  Applies for 3pl modeling. Information about which items share a common
  guessing parameter.

- control:

  List of control parameters for TAM estimation.

- n.plausible:

  Desired number of plausible values.

## Author

Sebastian Weirich

## Examples

``` r
################################################################################
###               Preparation: necessary for all examples                    ###
################################################################################

# load example data
# data set 'trends' contains item response data for three measurement occasions
# in a single data.frame
data(trends)

# second, create the q matrix from long format data
qMat <- unique(trends[ ,c("item","domain")])
qMat <- data.frame ( qMat[,"item", drop=FALSE], model.matrix(~domain-1, data = qMat))


################################################################################
###                Example 1: Unidimensional Rasch Model                     ###
################################################################################

# Example 1: define and run a unidimensional Rasch model for the 2010 cohort
# with all variables from booklet 02 (this booklet contains  reading and listening
# items together). For the estimation, Conquest is used.
# Since the wide-format dataset generated here will be used in later analyses,
# several covariates - such as gender and SES - are generated at the same time.
# These variables are not yet used in this example.
datW <- subset(trends, booklet == "Bo02" & year == 2010) |>
        reshape2::dcast(idstud+language+sex+ses+country~item, value.var = "value")

# defining the unidimensional model: specifying q matrix is not necessary
mod1 <- defineModel(dat=datW, items= -c(1:5), id="idstud", analysis.name = "unidim", dir = tempdir())
#> Error in defineModel(dat = datW, items = -c(1:5), id = "idstud", analysis.name = "unidim",     dir = tempdir()): Assertion on 'conquest.folder' failed: Directory expected, but file in place: 'C:/Users/grewered/AppData/Local/Programs/R/R-4.4.2/library/eatModel/exec/console_Feb2007.exe'.

# run the model
run1 <- runModel(mod1)
#> Error: object 'mod1' not found

# get the results
res1 <- getResults(run1)
#> Error: object 'run1' not found

# extract the item parameters from the results object
item1<- itemFromRes(res1)
#> Error: object 'res1' not found


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
# to explicitly specify item columns ... instead of saying 'items= -c(1:5)' which
# means: Everything except column 1 to 5 are item columns
items<- grep("^T[[:digit:]]{2}", colnames(datW))
mod1a<- defineModel(dat=datW, items= items, id="idstud", DIF.var = "sexNum",
        analysis.name = "unidimDIF", dir = tempdir())
#> Error in defineModel(dat = datW, items = items, id = "idstud", DIF.var = "sexNum",     analysis.name = "unidimDIF", dir = tempdir()): Assertion on 'conquest.folder' failed: Directory expected, but file in place: 'C:/Users/grewered/AppData/Local/Programs/R/R-4.4.2/library/eatModel/exec/console_Feb2007.exe'.

# run the model
run1a<- runModel(mod1a)
#> Error: object 'mod1a' not found

# get the results
res1a<- getResults(run1a)
#> Error: object 'run1a' not found


################################################################################
### Example 1b: Rasch Model with DIF variable in customized model statement  ###
################################################################################

# builing on example 1a, the DIF variable should be used without interaction term
# and without eliminating items which do not vary within DIF groups. We need the
# 'sexNum' variable which was defined previously
# Example based on 'model 6' in 'r:\Nawi\main\10_Studien\Auswertung\Alex Robitzsch R-Skripte\skalierung_conquest_mit_R_04.R'
items<- grep("^T[[:digit:]]{2}", colnames(datW))
mod1b<- defineModel(dat=datW, items= items, id="idstud", model.statement = "item + sexNum",
        analysis.name = "unidimCustom", dir = tempdir())
#> Error in defineModel(dat = datW, items = items, id = "idstud", model.statement = "item + sexNum",     analysis.name = "unidimCustom", dir = tempdir()): Assertion on 'conquest.folder' failed: Directory expected, but file in place: 'C:/Users/grewered/AppData/Local/Programs/R/R-4.4.2/library/eatModel/exec/console_Feb2007.exe'.

# run the model
run1b<- runModel(mod1b)
#> Error: object 'mod1b' not found

# get the results
res1b<- getResults(run1b)
#> Error: object 'run1b' not found

# item parameter
it1b <- itemFromRes(res1b)
#> Error: object 'res1b' not found


################################################################################
###        Example 2a: Multidimensional Rasch Model with anchoring           ###
################################################################################

# Example 2a: running a multidimensional Rasch model on a subset of items with latent
# regression. Use item parameter from the first model as anchor parameters

# read in anchor parameters from the results object of the first example
aPar <- itemFromRes(res1)[,c("item", "est")]
#> Error: object 'res1' not found

# defining the model: specifying q matrix now is necessary.
# Please note that all latent regression variables have to be of class numeric.
# If regression variables are factors, dummy variables automatically will be used.
# (This behavior is equivalent as in lm() for example.)
mod2a<- defineModel(dat=datW, items= grep("^T[[:digit:]]{2}", colnames(datW)),
        id="idstud", analysis.name = "twodim", HG.var = c("language","sex", "ses"),
        qMatrix = qMat, anchor = aPar, n.plausible = 20,dir = tempdir())
#> Error in defineModel(dat = datW, items = grep("^T[[:digit:]]{2}", colnames(datW)),     id = "idstud", analysis.name = "twodim", HG.var = c("language",         "sex", "ses"), qMatrix = qMat, anchor = aPar, n.plausible = 20,     dir = tempdir()): Assertion on 'conquest.folder' failed: Directory expected, but file in place: 'C:/Users/grewered/AppData/Local/Programs/R/R-4.4.2/library/eatModel/exec/console_Feb2007.exe'.

# run the model
run2a<- runModel(mod2a)
#> Error: object 'mod2a' not found

# get the results
res2a<- getResults(run2a)
#> Error: object 'run2a' not found


################################################################################
###        Example 2b: Multidimensional Rasch Model with equating            ###
################################################################################

# Example 2b: running a multidimensional Rasch model on a subset of items
# without anchoring. Defining the model: specifying q matrix now is necessary.
mod2b<- defineModel(dat=datW, items= grep("^T[[:digit:]]{2}", colnames(datW)),
        id="idstud", analysis.name = "twodim2", qMatrix = qMat,
        n.plausible = 20, dir = tempdir())
#> Error in defineModel(dat = datW, items = grep("^T[[:digit:]]{2}", colnames(datW)),     id = "idstud", analysis.name = "twodim2", qMatrix = qMat,     n.plausible = 20, dir = tempdir()): Assertion on 'conquest.folder' failed: Directory expected, but file in place: 'C:/Users/grewered/AppData/Local/Programs/R/R-4.4.2/library/eatModel/exec/console_Feb2007.exe'.

# run the model
run2b <- runModel(mod2b)
#> Error: object 'mod2b' not found

# get the results
res2b<- getResults(run2b)
#> Error: object 'run2b' not found

### equating (wenn nicht verankert)
eq2b <- equat1pl( results = res2b, prmNorm = aPar, difBound=.64, iterativ=TRUE)
#> Error: object 'res2b' not found

### transformation to the 'bista' metric: needs reference population definition
### we use some arbitrary values here
ref  <- data.frame ( domain = c("domainreading", "domainlistening"),
        m = c(0.03890191, 0.03587727), sd= c(1.219, 0.8978912))
cuts <- list ( domainreading = list ( values = 390+0:3*75),
        domainlistening = list ( values = 360+0:3*85))
tf2b <- transformToBista ( equatingList = eq2b, refPop = ref, cuts = cuts)
#> Error: object 'eq2b' not found


################################################################################
###            Example 3: Multidimensional Rasch Model in TAM                ###
################################################################################

# estimate model 2 with latent regression and anchored parameters in TAM
# specification of an output folder (via 'dir' argument) no longer necessary
mod3T<- defineModel(dat=datW, items= grep("^T[[:digit:]]{2}", colnames(datW)),
        id="idstud", qMatrix = qMat, HG.var = "sex", anchor = aPar, software = "tam")
#> Following 250 item(s) missed in data frame will be removed from Q matrix: 
#>     T07_08, T02_02, T07_01, T02_06, T07_05, T07_06, T02_07, T07_02, T07_04, T02_03, T07_09, T02_01, T02_04, T07_03, T02_05, T07_07, T07_10, T12_11, T09_12, T01_08, T13_15, T15_03, T13_04, T13_07, T13_16, T13_01, T15_02, T15_01, T15_15, T15_16, T15_09, T13_12, T13_13, T15_14, T13_08, T13_06, T15_12, T13_10, T13_17, T15_11, T13_05, T13_09, T15_10, T15_13, T03_07, T03_01, T03_08, T03_05, T03_06, T03_03, T10_08, T11_06, T11_03, T11_01, T04_02, T04_01, T04_06, T04_04, T11_04, T10_01, T11_08, T10_02, T10_09, T10_07, T11_02, T04_05, T04_03, T11_05, T10_06, T11_10, T11_09, T04_07, T15_17, T15_18, T13_18, T10_10, T04_09, T04_08, T11_0X, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07, T06_05, T05_03, T05_06, T05_05, T06_02, T06_04, T06_03, T05_04, T05_02, T05_01, T06_01, T08_04, T08_06, T08_03, T08_07, T08_02, T08_05, T08_01, T08_09, T08_08, T03_09, T06_1X, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> Background variable(s) 'sex' of class 
#>     'factor' will be converted to indicator variables.
#>     Variable 'sex' was converted to 1 indicator(s) with name(s) 'sexmale'.
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Error: object 'aPar' not found

# run the model
run3T<- runModel(mod3T)
#> Error: object 'mod3T' not found

# Object 'run2T' is of class 'tam.mml'
class(run3T)
#> Error: object 'run3T' not found

# the class of 'run2T' corresponds to the class defined by the TAM package; all
# functions of the TAM package intended for further processing (e.g. drawing
# plausible values, plotting deviance change etc.) work, for example:
wle  <- tam.wle(run3T)
#> Error: object 'run3T' not found

# Finally, the model result are collected in a single data frame
res3T<- getResults(run3T)
#> Error: object 'run3T' not found

# repeat the same model with mirt
items<- grep("^T[[:digit:]]{2}", colnames(datW), value=TRUE)
irtM <- data.frame(item = items, task = substr(items,1,4), stringsAsFactors = FALSE) |> dplyr::mutate(irtmod = "Rasch")
mod3M<- defineModel(dat=datW, items= grep("^T[[:digit:]]{2}", colnames(datW)),id="idstud",
        irtmodel = irtM[,c("item", "irtmod")],  anchor = aPar, qMatrix = qMat, HG.var = "sex",  software = "mirt")
#> Error in match.arg(arg = irtmodel, choices = c("1PL", "2PL", "PCM", "PCM2",     "RSM", "GPCM", "2PL.groups", "GPCM.design", "3PL")): 'arg' must be NULL or a character vector
run3M<- runModel(mod3M)
#> Error: object 'mod3M' not found
res3M<- getResults(run3M)
#> Error: object 'run3M' not found
itM  <- itemFromRes(res3M)
#> Error: object 'res3M' not found


################################################################################
###    Example 4: define und run multiple models defined by 'splitModels'    ###
################################################################################

# Example 4: define und run multiple models defined by 'splitModels'
# Model split is possible for different groups of items (i.e. domains) and/or
# different groups of persons (for example, federal states within Germany)

# define person grouping
pers  <- data.frame ( idstud = datW[,"idstud"] , group1 = datW[,"sex"],
         group2 = datW[,"country"], stringsAsFactors = FALSE )

# define 24 models, splitting according to person groups and item groups separately
# by default, multicore processing is applied
l1    <- splitModels ( qMatrix = qMat, person.groups = pers, nCores = 1)
#> ---------------------------------
#> splitModels: generating 24 models
#> ........................
#> see <returned>$models
#> number of cores: 1
#> ---------------------------------

# In this example, the sample sizes in the individual subgroups or models are far
# too small. The example is therefore intended solely to illustrate how the
# analysis can be broken down by item and respondent group.
table(datW[,c("sex", "country")])
#>         country
#> sex      countryA countryB countryC
#>   female       42       27       20
#>   male         32       35       16

# apply 'defineModel' for each of the 24 models in 'l1'
modMul<- defineModel(dat = datW, items = grep("^T[[:digit:]]{2}", colnames(datW)),
         id = "idstud", check.for.linking = TRUE, splittedModels = l1, software = "tam")
#> 
#> Specification of 'qMatrix' and 'person.groups' results in 24 model(s).
#> Error in eval(cl1[[u]]): object 'datW' not found

# run all models
runMul<- runModel(modMul)
#> Error: object 'modMul' not found

# get results of all models
resMul<- getResults(runMul)
#> Error: object 'runMul' not found


################################################################################
###          Example 5: Linking and equating for multiple models             ###
################################################################################

# Example 5: define und run multiple models according to different domains (item groups)
# and further linking/equating. This example mimics the routines necessary for the
# 'Vergleichsarbeiten' (VERA) at the Institute of Educational Progress (IQB)
# Let's assume that the 2015 cohort serves as the subsample for the 'Vergleichsarbeiten'.
vera <- subset(trends, booklet == "Bo02" & year == 2015) |>
        reshape2::dcast(idstud+language+sex+ses+country~item, value.var = "value")

# specify two models according to the two domains 'reading' and 'listening'
l2    <- splitModels ( qMatrix = qMat,  nCores = 1)
#> --------------------------------
#> splitModels: generating 2 models
#> ..
#> see <returned>$models
#> number of cores: 1
#> --------------------------------

# define 2 models. free estimation without anchoring
mods  <- defineModel(dat = vera, id = "idstud", check.for.linking = TRUE,
         splittedModels = l2, software = "tam")
#> 
#> Specification of 'qMatrix' and 'person.groups' results in 2 model(s).
#> 
#> 
#> =========================================
#> Model No. 1
#>     Model name:           domainlistening
#>     Number of items:      16
#>     Number of persons:    133
#>     Number of dimensions: 1
#> =========================================
#> 
#> Warning! Model No. 1, model name: 'domainlistening': 118 from 134 items listed the Q matrix not found in data:
#>     T12_02, T13_15, T15_03, T13_04, T13_07, T13_16, T13_01, T15_02, T15_01, T15_15, T15_16, T15_09, T13_12, T13_13, T15_14, T13_08, T13_06, T15_12, T13_10, T13_17, T15_11, T13_05, T13_09, T15_10, T15_13, T15_17, T15_18, T13_18, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> Following 118 item(s) missed in data frame will be removed from Q matrix: 
#>     T12_02, T13_15, T15_03, T13_04, T13_07, T13_16, T13_01, T15_02, T15_01, T15_15, T15_16, T15_09, T13_12, T13_13, T15_14, T13_08, T13_06, T15_12, T13_10, T13_17, T15_11, T13_05, T13_09, T15_10, T15_13, T15_17, T15_18, T13_18, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> 
#> 
#> =======================================
#> Model No. 2
#>     Model name:           domainreading
#>     Number of items:      17
#>     Number of persons:    133
#>     Number of dimensions: 1
#> =======================================
#> 
#> Warning! Model No. 2, model name: 'domainreading': 132 from 149 items listed the Q matrix not found in data:
#>     T07_08, T02_02, T07_01, T02_06, T07_05, T07_06, T02_07, T07_02, T07_04, T02_03, T07_09, T02_01, T02_04, T07_03, T02_05, T07_07, T07_10, T09_11, T01_05, T03_07, T03_01, T03_08, T03_05, T03_06, T03_03, T10_08, T11_06, T11_03, T11_01, T04_02, T04_01, T04_06, T04_04, T11_04, T10_01, T11_08, T10_02, T10_09, T10_07, T11_02, T04_05, T04_03, T11_05, T10_06, T11_10, T11_09, T04_07, T10_10, T04_09, T04_08, T11_0X, T06_05, T05_03, T05_06, T05_05, T06_02, T06_04, T06_03, T05_04, T05_02, T05_01, T06_01, T08_04, T08_06, T08_03, T08_07, T08_02, T08_05, T08_01, T08_09, T08_08, T03_09, T06_1X, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01
#> Following 132 item(s) missed in data frame will be removed from Q matrix: 
#>     T07_08, T02_02, T07_01, T02_06, T07_05, T07_06, T02_07, T07_02, T07_04, T02_03, T07_09, T02_01, T02_04, T07_03, T02_05, T07_07, T07_10, T09_11, T01_05, T03_07, T03_01, T03_08, T03_05, T03_06, T03_03, T10_08, T11_06, T11_03, T11_01, T04_02, T04_01, T04_06, T04_04, T11_04, T10_01, T11_08, T10_02, T10_09, T10_07, T11_02, T04_05, T04_03, T11_05, T10_06, T11_10, T11_09, T04_07, T10_10, T04_09, T04_08, T11_0X, T06_05, T05_03, T05_06, T05_05, T06_02, T06_04, T06_03, T05_04, T05_02, T05_01, T06_01, T08_04, T08_06, T08_03, T08_07, T08_02, T08_05, T08_01, T08_09, T08_08, T03_09, T06_1X, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01
#> 1 subject(s) solved each item: P04491 (17 correct), P04491 (17 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).

# run 2 models
runs  <- runModel(mods)

# get the results
ress  <- getResults(runs)
#> |*****|
#> |-----|
#> |*****|
#> |-----|

# norm parameters are the item parameters from example 1
aPar  <- itemFromRes(res1)[,c("item", "est")]
#> Error: object 'res1' not found

# anchoring without exclusion of linking DIF items (DIF items will only be identified)
# the linking error is well above the critical threshold of 0.1
anch  <- equat1pl ( results = ress, prmNorm = aPar, excludeLinkingDif = FALSE,
         difBound = 0.64)
#> Found 2 model(s).
#>    Equating is executed for each dimension in each model separately.
#> Error: object 'aPar' not found

# anchoring with exclusion of linking DIF items: especially for listening, the
# linking constant changes substantially when linking DIF items are excluded
anch2 <- equat1pl ( results = ress, prmNorm = aPar, excludeLinkingDif = TRUE,
         difBound = 0.64, iterativ = FALSE)
#> Found 2 model(s).
#>    Equating is executed for each dimension in each model separately.
#> Error: object 'aPar' not found

# anchoring with iterative exclusion of linking DIF items
anch3 <- equat1pl ( results = ress, prmNorm = aPar, excludeLinkingDif = TRUE,
         difBound = 0.64, iterativ = TRUE)
#> Found 2 model(s).
#>    Equating is executed for each dimension in each model separately.
#> Error: object 'aPar' not found

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
#> Error: object 'anch3' not found
head(dfr$itempars)
#> Error: object 'dfr' not found
head(dfr$personpars)
#> Error: object 'dfr' not found


################################################################################
###      Example 5a: Linking for multiple models, including global domain    ###
################################################################################

# Example 5a: define und run multiple models according to different domains (item groups)
# and further linking/equating. Same as example 5, but extended for the 'global'
# domain.

# add the 'global' domain (reading and listening together) in the Q matrix
qMat2 <- qMat |> dplyr::mutate(global = 1)

# specify two models according to the two domains 'reading' and 'listening'
l3    <- splitModels ( qMatrix = qMat2, nCores = 1)
#> --------------------------------
#> splitModels: generating 3 models
#> ...
#> see <returned>$models
#> number of cores: 1
#> --------------------------------

# define 3 models
mods5 <- defineModel(dat = vera, id = "idstud", check.for.linking = TRUE,
         splittedModels = l3, software = "tam")
#> 
#> Specification of 'qMatrix' and 'person.groups' results in 3 model(s).
#> 
#> 
#> =========================================
#> Model No. 1
#>     Model name:           domainlistening
#>     Number of items:      16
#>     Number of persons:    133
#>     Number of dimensions: 1
#> =========================================
#> 
#> Warning! Model No. 1, model name: 'domainlistening': 118 from 134 items listed the Q matrix not found in data:
#>     T12_02, T13_15, T15_03, T13_04, T13_07, T13_16, T13_01, T15_02, T15_01, T15_15, T15_16, T15_09, T13_12, T13_13, T15_14, T13_08, T13_06, T15_12, T13_10, T13_17, T15_11, T13_05, T13_09, T15_10, T15_13, T15_17, T15_18, T13_18, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> Following 118 item(s) missed in data frame will be removed from Q matrix: 
#>     T12_02, T13_15, T15_03, T13_04, T13_07, T13_16, T13_01, T15_02, T15_01, T15_15, T15_16, T15_09, T13_12, T13_13, T15_14, T13_08, T13_06, T15_12, T13_10, T13_17, T15_11, T13_05, T13_09, T15_10, T15_13, T15_17, T15_18, T13_18, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> 
#> 
#> =======================================
#> Model No. 2
#>     Model name:           domainreading
#>     Number of items:      17
#>     Number of persons:    133
#>     Number of dimensions: 1
#> =======================================
#> 
#> Warning! Model No. 2, model name: 'domainreading': 132 from 149 items listed the Q matrix not found in data:
#>     T07_08, T02_02, T07_01, T02_06, T07_05, T07_06, T02_07, T07_02, T07_04, T02_03, T07_09, T02_01, T02_04, T07_03, T02_05, T07_07, T07_10, T09_11, T01_05, T03_07, T03_01, T03_08, T03_05, T03_06, T03_03, T10_08, T11_06, T11_03, T11_01, T04_02, T04_01, T04_06, T04_04, T11_04, T10_01, T11_08, T10_02, T10_09, T10_07, T11_02, T04_05, T04_03, T11_05, T10_06, T11_10, T11_09, T04_07, T10_10, T04_09, T04_08, T11_0X, T06_05, T05_03, T05_06, T05_05, T06_02, T06_04, T06_03, T05_04, T05_02, T05_01, T06_01, T08_04, T08_06, T08_03, T08_07, T08_02, T08_05, T08_01, T08_09, T08_08, T03_09, T06_1X, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01
#> Following 132 item(s) missed in data frame will be removed from Q matrix: 
#>     T07_08, T02_02, T07_01, T02_06, T07_05, T07_06, T02_07, T07_02, T07_04, T02_03, T07_09, T02_01, T02_04, T07_03, T02_05, T07_07, T07_10, T09_11, T01_05, T03_07, T03_01, T03_08, T03_05, T03_06, T03_03, T10_08, T11_06, T11_03, T11_01, T04_02, T04_01, T04_06, T04_04, T11_04, T10_01, T11_08, T10_02, T10_09, T10_07, T11_02, T04_05, T04_03, T11_05, T10_06, T11_10, T11_09, T04_07, T10_10, T04_09, T04_08, T11_0X, T06_05, T05_03, T05_06, T05_05, T06_02, T06_04, T06_03, T05_04, T05_02, T05_01, T06_01, T08_04, T08_06, T08_03, T08_07, T08_02, T08_05, T08_01, T08_09, T08_08, T03_09, T06_1X, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01
#> 1 subject(s) solved each item: P04491 (17 correct), P04491 (17 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> 
#> 
#> ================================
#> Model No. 3
#>     Model name:           global
#>     Number of items:      33
#>     Number of persons:    133
#>     Number of dimensions: 1
#> ================================
#> 
#> Warning! Model No. 3, model name: 'global': 250 from 283 items listed the Q matrix not found in data:
#>     T12_02, T07_08, T02_02, T07_01, T02_06, T07_05, T07_06, T02_07, T07_02, T07_04, T02_03, T07_09, T02_01, T02_04, T07_03, T02_05, T07_07, T07_10, T09_11, T01_05, T13_15, T15_03, T13_04, T13_07, T13_16, T13_01, T15_02, T15_01, T15_15, T15_16, T15_09, T13_12, T13_13, T15_14, T13_08, T13_06, T15_12, T13_10, T13_17, T15_11, T13_05, T13_09, T15_10, T15_13, T03_07, T03_01, T03_08, T03_05, T03_06, T03_03, T10_08, T11_06, T11_03, T11_01, T04_02, T04_01, T04_06, T04_04, T11_04, T10_01, T11_08, T10_02, T10_09, T10_07, T11_02, T04_05, T04_03, T11_05, T10_06, T11_10, T11_09, T04_07, T15_17, T15_18, T13_18, T10_10, T04_09, T04_08, T11_0X, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07, T06_05, T05_03, T05_06, T05_05, T06_02, T06_04, T06_03, T05_04, T05_02, T05_01, T06_01, T08_04, T08_06, T08_03, T08_07, T08_02, T08_05, T08_01, T08_09, T08_08, T03_09, T06_1X, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> Following 250 item(s) missed in data frame will be removed from Q matrix: 
#>     T12_02, T07_08, T02_02, T07_01, T02_06, T07_05, T07_06, T02_07, T07_02, T07_04, T02_03, T07_09, T02_01, T02_04, T07_03, T02_05, T07_07, T07_10, T09_11, T01_05, T13_15, T15_03, T13_04, T13_07, T13_16, T13_01, T15_02, T15_01, T15_15, T15_16, T15_09, T13_12, T13_13, T15_14, T13_08, T13_06, T15_12, T13_10, T13_17, T15_11, T13_05, T13_09, T15_10, T15_13, T03_07, T03_01, T03_08, T03_05, T03_06, T03_03, T10_08, T11_06, T11_03, T11_01, T04_02, T04_01, T04_06, T04_04, T11_04, T10_01, T11_08, T10_02, T10_09, T10_07, T11_02, T04_05, T04_03, T11_05, T10_06, T11_10, T11_09, T04_07, T15_17, T15_18, T13_18, T10_10, T04_09, T04_08, T11_0X, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07, T06_05, T05_03, T05_06, T05_05, T06_02, T06_04, T06_03, T05_04, T05_02, T05_01, T06_01, T08_04, T08_06, T08_03, T08_07, T08_02, T08_05, T08_01, T08_09, T08_08, T03_09, T06_1X, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).

# run 3 models
runs5 <- runModel(mods5)

# get the results
res5  <- getResults(runs5)
#> |*****|
#> |-----|
#> |*****|
#> |-----|
#> |*****|
#> |-----|

# only for illustration, we create arbitrary 'normed' parameters for anchoring.
# Each item now has two item parameter: one is domain-specific, one is the global
# item parameter. Hence, each item occurs two times in 'prmNrm'
prmNrm<- do.call("rbind", by(data = itemFromRes(ress),
         INDICES = itemFromRes(ress)[,"dimension"], FUN = function (d) {
             d    <- d[sample(1:nrow(d), size = round(nrow(d)/3), replace=FALSE), c("dimension", "item", "est")]
             const<- ifelse(d[1,"dimension"] == "domainlistening", 0.6, 0.4)
             err  <- ifelse(d[1,"dimension"] == "domainlistening", 0.75, 0.6)
             d[,"est"] <- d[,"est"] - const + rnorm ( nrow(d), 0, err)
             return(d) }))
prmGlo<- prmNrm |> dplyr::mutate(dimension = "global")
prmNrm<- rbind ( prmNrm, prmGlo)

# if the item identifier in 'prmNrm' is not unique, 'equat1pl' has to know which
# parameter in 'prmNrm' belongs to which dimension/domain. Hence, the dimension
# is added to 'prmNrm'.

# anchoring: if 'prmNrm' has more than 2 columns, the columns of 'prmNrm' must be
# specified in 'equat1pl'
anch5 <- equat1pl ( results = res5, prmNorm = prmNrm, item = "item", value = "est",
         domain = "dimension", excludeLinkingDif = FALSE, difBound = 0.64)
#> Found 3 model(s).
#>    Equating is executed for each dimension in each model separately.
#> 
#> ====================================================================================================
#>  
#> Model No. 1
#>     Model name:                domainlistening
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainlistening
#>     Number of linking items:   5
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainlistening': 1 of 5 items with linking |DIF| > 0.64 identified.
#> 
#>     item   dif linking.constant linkerror
#> 1 T14_04 1.464            -0.37     0.377
#> 
#> 
#> ====================================================================================================
#>  
#> Model No. 2
#>     Model name:                domainreading
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainreading
#>     Number of linking items:   6
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainreading': 3 of 6 items with linking |DIF| > 0.64 identified.
#> 
#>     item    dif linking.constant linkerror
#> 1 T01_01  0.664            0.067     0.232
#> 2 T01_02 -0.669            0.067     0.232
#> 3 T09_05  0.751            0.067     0.232
#> 
#> 
#> ====================================================================================================
#>  
#> Model No. 3
#>     Model name:                global
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   global
#>     Number of linking items:   11
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'global': 4 of 11 items with linking |DIF| > 0.64 identified.
#> 
#>     item    dif linking.constant linkerror
#> 1 T01_02 -0.782           -0.157     0.209
#> 2 T09_05  0.648           -0.157     0.209
#> 3 T12_09  0.791           -0.157     0.209
#> 4 T14_04  1.585           -0.157     0.209
#> 

# for illustration: consequences of linking the global domain if linking constants
# differ between domains
anch5a<- equat1pl ( results = res5, prmNorm = prmNrm, item = "item", value = "est",
         domain = "dimension", excludeLinkingDif = TRUE, iterativ = TRUE, difBound = 0.64)
#> Found 3 model(s).
#>    Equating is executed for each dimension in each model separately.
#> 
#> ====================================================================================================
#>  
#> Model No. 1
#>     Model name:                domainlistening
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainlistening
#>     Number of linking items:   5
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainlistening': 1 of 5 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T14_04'.
#>    Iteration 2: Exclude item 'T12_09'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                     -0.370     0.377
#> 2 iterativ    1       T14_04        1.464           -0.050     0.259
#> 3 iterativ    2       T12_09        0.714            0.175     0.181
#> 
#> 
#> ====================================================================================================
#>  
#> Model No. 2
#>     Model name:                domainreading
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainreading
#>     Number of linking items:   6
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainreading': 3 of 6 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T09_05'.
#>    Iteration 2: Exclude item 'T01_01'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.067     0.232
#> 2 iterativ    1       T09_05        0.751            0.202     0.231
#> 3 iterativ    2       T01_01        0.813            0.384     0.185
#> 
#> 
#> ====================================================================================================
#>  
#> Model No. 3
#>     Model name:                global
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   global
#>     Number of linking items:   11
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'global': 4 of 11 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T14_04'.
#>    Iteration 2: Exclude item 'T12_09'.
#>    Iteration 3: Exclude item 'T09_05'.
#>    Iteration 4: Exclude item 'T01_01'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                     -0.157     0.209
#> 2 iterativ    1       T14_04        1.585           -0.014     0.169
#> 3 iterativ    2       T12_09        0.847            0.072     0.162
#> 4 iterativ    3       T09_05        0.789            0.162     0.153
#> 5 iterativ    4       T01_01        0.821            0.267     0.128
#> 

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
dfr   <- transformToBista(equatingList = anch5, refPop = refPop, cuts=cuts )
#> The 'refPop' data.frame does not include information about reference population mean/SD on Bista metric. Values will be defaulted to 500/100.
#> Found duplicated entries in 'item-ID' column. This should only occur for subject 'math' in grade 3.
#> Cannot find 'global' entry in the 'domain' column. Cancel reshaping.
head(dfr$itempars)
#>             model   item       dimension Nvalid     itemP itemDiscrim       est
#> 1 domainlistening T12_01 domainlistening    133 0.9398496  0.02099527 -3.015202
#> 2 domainlistening T12_03 domainlistening    133 0.8872180  0.37141881 -2.290056
#> 3 domainlistening T12_04 domainlistening    133 0.8796992  0.23727591 -2.211710
#> 4 domainlistening T12_05 domainlistening    133 0.9849624  0.10028410 -4.485206
#> 5 domainlistening T12_07 domainlistening    133 0.8646617  0.33255757 -2.066257
#> 6 domainlistening T12_09 domainlistening    133 0.9548872  0.12657818 -3.330597
#>          se     infit    outfit estTransf estTransf625 estTransfBista
#> 1 0.3725945 1.0539760 1.2093755 -3.385366    -2.874540      175.86079
#> 2 0.2836654 0.9163307 0.7550022 -2.660219    -2.149393      256.62189
#> 3 0.2762515 0.9862945 0.9141946 -2.581874    -2.071048      265.34737
#> 4 0.7170350 1.0056418 0.7545260 -4.855369    -4.344543       12.14348
#> 5 0.2634642 0.9313443 0.8899342 -2.436420    -1.925595      281.54680
#> 6 0.4249207 1.0059954 1.0341609 -3.700760    -3.189935      140.73466
#>   traitLevel linkingConstant linkingMethod nLinkitems linkingError
#> 1          1      -0.3701634     Mean.Mean          5    0.3774444
#> 2          1      -0.3701634     Mean.Mean          5    0.3774444
#> 3          1      -0.3701634     Mean.Mean          5    0.3774444
#> 4          1      -0.3701634     Mean.Mean          5    0.3774444
#> 5          1      -0.3701634     Mean.Mean          5    0.3774444
#> 6          1      -0.3701634     Mean.Mean          5    0.3774444
#>   linkingErrorTransfBista linkingErrorTraitLevel    refMean     refSD
#> 1                42.03676             0.08832631 0.03587727 0.8978912
#> 2                42.03676             0.08832631 0.03587727 0.8978912
#> 3                42.03676             0.08832631 0.03587727 0.8978912
#> 4                42.03676             0.08832631 0.03587727 0.8978912
#> 5                42.03676             0.08832631 0.03587727 0.8978912
#> 6                42.03676             0.08832631 0.03587727 0.8978912
#>   refTransfMean refTransfSD
#> 1           500         100
#> 2           500         100
#> 3           500         100
#> 4           500         100
#> 5           500         100
#> 6           500         100
head(dfr$personpars)
#>             model idstud           group imp      value valueTransfBista
#> 1 domainlistening P04477 domainlistening pv1 -0.7975125         365.9578
#> 2 domainlistening P04478 domainlistening pv1 -0.2713932         424.5528
#> 3 domainlistening P06981 domainlistening pv5 -0.7468397         371.6013
#> 4 domainlistening P04480 domainlistening pv1 -0.7048242         376.2807
#> 5 domainlistening P07232 domainlistening pv3 -0.4266830         407.2578
#> 6 domainlistening P07233 domainlistening pv3 -0.2941525         422.0180
#>   traitLevel       dimension linkingError linkingErrorTransfBista
#> 1          2 domainlistening    0.3774444                42.03676
#> 2          2 domainlistening    0.3774444                42.03676
#> 3          2 domainlistening    0.3774444                42.03676
#> 4          2 domainlistening    0.3774444                42.03676
#> 5          2 domainlistening    0.3774444                42.03676
#> 6          2 domainlistening    0.3774444                42.03676
#>   linkingErrorTraitLevel
#> 1             0.07658122
#> 2             0.07658122
#> 3             0.07658122
#> 4             0.07658122
#> 5             0.07658122
#> 6             0.07658122


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

# Preparation: assume time of measurement 't1' corresponds to the year 2010.
# This is the year of the reference population. We consider both domains,
# reading and listening in a single call. Hence, the data.frame 'datT1' contains
# items of both domains. We add weights to the data as the definition of mean
# and sd of the reference population (population of 't1') necessitates weights.
datT1<- reshape2::dcast(subset ( trends, year == 2010),
        idstud+country+sex+ses+language+wgt~item, value.var="value")

# generate Q matrix
qMat <- unique(trends[ ,c("item","domain")])
qMat <- data.frame ( qMat[,"item", drop=FALSE], model.matrix(~domain-1, data = qMat))

# First step: item calibration in separate unidimensional models for each domain
# (the models are short and simple, so we don't need multicore)
modsT1<- splitModels ( qMatrix = qMat, nCores = 1)
#> --------------------------------
#> splitModels: generating 2 models
#> ..
#> see <returned>$models
#> number of cores: 1
#> --------------------------------

# define 2 models. Note: not all items of the Q matrix are present in the data.
# Items which occur only in the Q matrix will be ignored. 
defT1 <- defineModel(dat = datT1, id = "idstud", check.for.linking = TRUE,
         splittedModels = modsT1, software = "tam")
#> 
#> Specification of 'qMatrix' and 'person.groups' results in 2 model(s).
#> 
#> 
#> =========================================
#> Model No. 1
#>     Model name:           domainlistening
#>     Number of items:      51
#>     Number of persons:    4476
#>     Number of dimensions: 1
#> =========================================
#> 
#> Warning! Model No. 1, model name: 'domainlistening': 83 from 134 items listed the Q matrix not found in data:
#>     T12_11, T15_17, T15_18, T13_18, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> Following 83 item(s) missed in data frame will be removed from Q matrix: 
#>     T12_11, T15_17, T15_18, T13_18, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> 76 subject(s) do not solve any item:
#>    P00009 (6 false), P00864 (6 false), P02759 (21 false) ... 
#> 55 subject(s) solved each item: P00099 (10 correct), P00404 (10 correct), P02659 (11 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> 
#> 
#> =======================================
#> Model No. 2
#>     Model name:           domainreading
#>     Number of items:      80
#>     Number of persons:    4476
#>     Number of dimensions: 1
#> =======================================
#> 
#> Warning! Model No. 2, model name: 'domainreading': 69 from 149 items listed the Q matrix not found in data:
#>     T09_12, T01_08, T10_10, T04_09, T04_08, T11_0X, T08_09, T08_08, T03_09, T06_1X, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01
#> Following 69 item(s) missed in data frame will be removed from Q matrix: 
#>     T09_12, T01_08, T10_10, T04_09, T04_08, T11_0X, T08_09, T08_08, T03_09, T06_1X, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01
#> 17 subject(s) do not solve any item:
#>    P00173 (6 false), P01137 (7 false), P02863 (23 false) ... 
#> 97 subject(s) solved each item: P00078 (6 correct), P00389 (10 correct), P03338 (27 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).

# run (calibrate) the 2 models subsequently
runT1 <- runModel(defT1)

# get the results of the two unidimensional models
resT1 <- getResults(runT1, omitWle = TRUE, Q3 = FALSE)
#> |*****|
#> |-----|
#> |*****|
#> |-----|

# extract item parameters from the 'results' object
# t1 is the reference measurement occasion, i.e. no linking/equating is necessary
itemT1<- itemFromRes(resT1)

# The measurement point 't1' defines the reference population. We need to specify
# mean and SD, using weights.
eqRef <- equat1pl(resT1)
#> Found 2 model(s).
#>    Equating is executed for each dimension in each model separately.
#> No norm parameter defined ('prmNorm' is missing). Treat current sample as drawn from the reference population.

# transformation to the 'bista' metric
# Note: if the sample was drawn from the reference population in a weighted sample,
# mean and SD are not yet known. So we ignore the 'refPop' argument in 'transformToBista'
# and simply define the cut scores. The function then assumes that the parameter
# stem from the reference population and estimates its mean and sd.
cuts  <- list ( domainreading = list ( values = 390+0:3*75),
         domainlistening = list ( values = 360+0:3*85))
tfRef <- transformToBista ( equatingList = eqRef, cuts=cuts, vera=FALSE,
         weights = datT1[,c("idstud", "wgt")] )
#> 'refPop' was not defined. Treat current sample as drawn from the reference population.
#> The 'refPop' data.frame does not include information about reference population mean/SD on Bista metric. Values will be defaulted to 500/100.

# Second step: drawing plausible values separately for each country.
# A two-dimensional (reading/listening) model is specified separately for each
# person group (= each country) with item parameters fixed at their calibration
# values. Moreover, a latent regression model is used (in the actual 'Laendervergleich',
# regressors are principal components of previously imputed background variables). We use 'sex',
# 'ses' and 'language' as regressors. For convenience, 'ses' is scaled (mean = 0, sd = 1)
datT1[,"ses_scaled"] <- scale(datT1[,"ses"])[,1]

# have a look at 'sex' and 'language at home' in each country:
table(datT1[,c("country", "sex")])
#>           sex
#> country    female male
#>   countryA    787  811
#>   countryB    663  646
#>   countryC    795  774
table(datT1[,c("country", "language")])
#>           language
#> country    native nativeAndOther other
#>   countryA   1143            437    18
#>   countryB   1218             81    10
#>   countryC   1248            308    13

# Running second step: split models according to person groups
# ('all.persons' must be FALSE, otherwise the whole group would be treated as
# a separate distinct group.)
modT1P<- splitModels ( person.groups = datT1[,c("idstud", "country")],
         all.persons = FALSE, nCores = 1)
#> --------------------------------
#> splitModels: generating 3 models
#> ...
#> see <returned>$models
#> number of cores: 1
#> --------------------------------

# define the 3 country-specific 2-dimensional models, specifying latent regression
# model and fixed item parameters.
defT1P<- defineModel(dat = datT1, items = itemT1[,"item"], id = "idstud",
         check.for.linking = TRUE, splittedModels = modT1P, qMatrix = qMat,
         anchor = itemT1[,c("item", "est")],
         HG.var = c("sex", "ses_scaled", "language"),  software = "tam")
#> 
#> Specification of 'qMatrix' and 'person.groups' results in 3 model(s).
#> Error in eval(cl1[[u]]): object 'itemT1' not found

# run the 3 models (estimation needs approx. 20 seconds)
runT1P<- runModel(defT1P)
#> Error: object 'defT1P' not found

# get the results (to save time, item fit estimation is skipped)
resT1P<- getResults(runT1P, omitWle = TRUE, Q3 = FALSE)
#> Error: object 'runT1P' not found

# latent regression coefficients for the three countries and two dimensions
regcoefFromRes(resT1P, digits = 3)
#> Error: object 'resT1P' not found

# equating is not necessary, as the models run with fixed item parameters
# However, to prepare for the transformation on the 'bista' metric, run
# 'equat1pl' with empty arguments
ankT1P<- equat1pl ( results = resT1P)
#> Error: object 'resT1P' not found

# transformation to the 'bista' metric
# Note: if the sample was drawn from the reference population, mean and SD
# were just computed and captured in 'tfRef'.
dfrT1P<- transformToBista ( equatingList = ankT1P, refPop = tfRef[["refPop"]][,-2], cuts=cuts, vera=FALSE )
#> Error: object 'ankT1P' not found


################################################################################
###     Example 6a: Extend example 6 with trend estimation (now for t2)      ###
################################################################################

# Example 6a needs the objects (Q matrix, item parameters, ...) created in example 6
# Preparation: assume time of measurement 't2'.
datT2<- reshape2::dcast(subset ( trends, year == 2015),
        idstud+country+sex+ses+language~item, value.var="value")

# First step: item calibration in separate unidimensional models for each domain
modsT2<- splitModels ( qMatrix = qMat, nCores = 1)
#> --------------------------------
#> splitModels: generating 2 models
#> ..
#> see <returned>$models
#> number of cores: 1
#> --------------------------------

# define 2 models. Items which occur only in the Q matrix but not in the data
# will be ignored.
defT2 <- defineModel(dat = datT2, id = "idstud", check.for.linking = TRUE,
         splittedModels = modsT2, software = "tam")
#> 
#> Specification of 'qMatrix' and 'person.groups' results in 2 model(s).
#> 
#> 
#> =========================================
#> Model No. 1
#>     Model name:           domainlistening
#>     Number of items:      96
#>     Number of persons:    4516
#>     Number of dimensions: 1
#> =========================================
#> 
#> Warning! Model No. 1, model name: 'domainlistening': 38 from 134 items listed the Q matrix not found in data:
#>     T12_02, T13_16, T15_12, T15_11, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01
#> Following 38 item(s) missed in data frame will be removed from Q matrix: 
#>     T12_02, T13_16, T15_12, T15_11, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01
#> Found 105 cases with missings on all items.
#> Cases with missings on all items will be deleted.
#> 149 subject(s) do not solve any item:
#>    P04521 (6 false), P07332 (6 false), P08038 (19 false) ... 
#> 77 subject(s) solved each item: P04599 (10 correct), P06913 (10 correct), P06282 (13 correct) ... 
#> W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).
#>                   501 cases unconnected. Following items are unconnected: 
#>                   T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> 
#> 
#> =======================================
#> Model No. 2
#>     Model name:           domainreading
#>     Number of items:      119
#>     Number of persons:    4516
#>     Number of dimensions: 1
#> =======================================
#> 
#> Warning! Model No. 2, model name: 'domainreading': 30 from 149 items listed the Q matrix not found in data:
#>     T09_11, T01_05, T03_07, T10_08, T11_03, T04_01, T11_04, T04_07, T11_0X, T06_04, T06_03, T08_04, T08_05, T06_1X, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03
#> Following 30 item(s) missed in data frame will be removed from Q matrix: 
#>     T09_11, T01_05, T03_07, T10_08, T11_03, T04_01, T11_04, T04_07, T11_0X, T06_04, T06_03, T08_04, T08_05, T06_1X, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03
#> 56 subject(s) do not solve any item:
#>    P04511 (6 false), P05616 (10 false), P08143 (14 false) ... 
#> 114 subject(s) solved each item: P04509 (6 correct), P08241 (8 correct), P05805 (20 correct) ... 
#> W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).
#>                   4035 cases unconnected. Following items are unconnected: 
#>                   T02_02, T02_06, T02_07, T02_03, T02_01, T02_04, T02_05, T09_03, T09_05, T01_03, T09_02, T01_06, T01_01, T09_10, T01_02, T09_06, T01_04, T09_09, T09_08, T09_07, T09_04, T01_07, T09_12, T01_08, T03_01, T03_08, T03_05, T03_06, T03_03, T11_06, T11_01, T04_02, T04_06, T04_04, T10_01, T11_08, T10_02, T10_09, T10_07, T11_02, T04_05, T04_03, T11_05, T10_06, T11_10, T11_09, T10_10, T04_09, T04_08, T06_05, T05_03, T05_06, T05_05, T06_02, T05_04, T05_02, T05_01, T06_01, T08_06, T08_03, T08_07, T08_02, T08_01, T08_09, T08_08, T03_09, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).

# run 2 models for calibration
runT2 <- runModel(defT2)

# get the results
resT2 <- getResults(runT2)
#> |*****|
#> |-----|
#> |*****|
#> |-----|

# collect item parameters
itemT2<- itemFromRes(resT2)

# Second step: compute linking constant between 't1' and 't2' with the iterative
# exclusion of linking DIF items and computation of linking error. We use the
# 'itemT1' object created in example 6 for reference item parameters. The linking
# procedure is executed consecutively for listening and reading.
L.t1t2<- equat1pl ( results = resT2, prmNorm = itemT1[,c("item", "est")],
         excludeLinkingDif = TRUE, difBound = 0.64, iterativ = TRUE)
#> Found 2 model(s).
#>    Equating is executed for each dimension in each model separately.
#> 
#> ====================================================================================================
#>  
#> Model No. 1
#>     Model name:                domainlistening
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainlistening
#>     Number of linking items:   36
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainlistening': 2 of 36 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T15_02'.
#>    Iteration 2: Exclude item 'T14_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                     -0.243     0.082
#> 2 iterativ    1       T15_02       -1.873           -0.298     0.063
#> 3 iterativ    2       T14_04        1.675           -0.250     0.042
#> 
#> 
#> ====================================================================================================
#>  
#> Model No. 2
#>     Model name:                domainreading
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainreading
#>     Number of linking items:   68
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainreading': 4 of 68 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T02_07'.
#>    Iteration 2: Exclude item 'T01_07'.
#>    Iteration 3: Exclude item 'T02_04'.
#>    Iteration 4: Exclude item 'T05_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                     -0.232     0.061
#> 2 iterativ    1       T02_07        2.505           -0.195     0.050
#> 3 iterativ    2       T01_07        2.267           -0.160     0.036
#> 4 iterativ    3       T02_04        1.176           -0.142     0.032
#> 5 iterativ    4       T05_04        0.867           -0.129     0.029
#> 

# linking constant is negative: students performance at T2 is worse than T1
# see that DIF is non symmetrical for reading: the four items with the highest
# amount of DIF have all positive DIF values. DIF exclusion hence shrinks linking
# constant towards zero.
# Third step: transform item parameters of 't2' to the metric of 't1'
# We now need to specify the 'refPop' argument. We use the values from 't1' which
# serves as the reference. To capture linking errors in a separate data.frame
# within the returned list, we define the years of assessment
ref   <- tfRef[["refPop"]][,-2]
T.t1t2<- transformToBista ( equatingList = L.t1t2, refPop=ref, cuts = cuts,
         vera=FALSE, years = c(2010,2015))
#> Error in transformToBista(equatingList = L.t1t2, refPop = ref, cuts = cuts,     vera = FALSE, years = c(2010, 2015)): Invalid 'refPop'.

# The object 'T.t1t2' now contains transformed person and item parameters with
# original and transformed linking errors. See for example person parameter:
head(T.t1t2$personpars)
#> Error: object 'T.t1t2' not found

# Fourth step: drawing plausible values for 't2'. We use the transformed item
# parameters (captured in 'T.t1t2') for anchoring

# Running second step: split models according to person groups (countries)
# ('all.persons' must be FALSE, otherwise the whole group would be treated as
# a separate distinct group.)
modT2P<- splitModels ( person.groups = datT2[,c("idstud", "country")] ,
         all.persons = FALSE, nCores = 1)
#> --------------------------------
#> splitModels: generating 3 models
#> ...
#> see <returned>$models
#> number of cores: 1
#> --------------------------------

# define the 3 country-specific 2-dimensional models, specifying latent regression
# model and fixed item parameters. We used the transformed item parameters (captured
# in 'T.t1t2[["itempars"]]' --- using the 'estTransf' column) for anchoring.
# Again, 'ses' is scaled (mean = 0, sd = 1)
datT2[,"ses_scaled"] <- scale(datT2[,"ses"])[,1]
defT2P<- defineModel(dat = datT2, items = itemT2[,"item"], id = "idstud",
         check.for.linking = TRUE, splittedModels = modT2P, qMatrix = qMat,
         anchor = T.t1t2[["itempars"]][,c("item", "estTransf")],
         HG.var = c("sex", "ses_scaled", "language"), software = "tam")
#> 
#> Specification of 'qMatrix' and 'person.groups' results in 3 model(s).
#> Error in eval(cl1[[u]]): object 'itemT2' not found

# run the 3 models (estimation takes approx. 29 seconds)
runT2P<- runModel(defT2P)
#> Error: object 'defT2P' not found

# get the results
resT2P<- getResults(runT2P)
#> Error: object 'runT2P' not found

# equating is not necessary, as the models run with fixed item parameters
# However, to prepare for the transformation on the 'bista' metric, run
# 'equat1pl' with empty arguments
ankT2P<- equat1pl ( results = resT2P)
#> Error: object 'resT2P' not found

# transformation to the 'bista' metric, using the previously defined cut scores
# and the reference population mean and sd from 't1'
dfrT2P<- transformToBista ( equatingList = ankT2P, refPop=ref, cuts=cuts, vera=FALSE)
#> Error: object 'ankT2P' not found


################################################################################
###     Example 6b: Extend example 6a with trend estimation (now for t3)     ###
################################################################################

# Example 6b needs the objects (Q matrix, item parameters, ...) created in example 6
# and 6a. Preparation: assume time of measurement 't3'.
datT3<- reshape2::dcast(subset ( trends, year == 2020),
        idstud+country+sex+ses+language~item, value.var="value")

# First step: item calibration in separate unidimensional models for each domain
modsT3<- splitModels ( qMatrix = qMat, nCores = 1)
#> --------------------------------
#> splitModels: generating 2 models
#> ..
#> see <returned>$models
#> number of cores: 1
#> --------------------------------

# define 2 models. Items which occur only in the Q matrix but not in the data
# will be ignored.
defT3 <- defineModel(dat = datT3, id = "idstud", check.for.linking = TRUE,
         splittedModels = modsT3, software = "tam")
#> 
#> Specification of 'qMatrix' and 'person.groups' results in 2 model(s).
#> 
#> 
#> =========================================
#> Model No. 1
#>     Model name:           domainlistening
#>     Number of items:      119
#>     Number of persons:    4532
#>     Number of dimensions: 1
#> =========================================
#> 
#> Warning! Model No. 1, model name: 'domainlistening': 15 from 134 items listed the Q matrix not found in data:
#>     T12_02, T13_16, T15_12, T15_11, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07
#> Following 15 item(s) missed in data frame will be removed from Q matrix: 
#>     T12_02, T13_16, T15_12, T15_11, T16_08, T16_09, T16_11, T16_13, T16_01, T16_04, T16_10, T16_02, T16_12, T16_06, T16_07
#> Found 212 cases with missings on all items.
#> Cases with missings on all items will be deleted.
#> 185 subject(s) do not solve any item:
#>    P09749 (6 false), P12158 (6 false), P11510 (24 false) ... 
#> 87 subject(s) solved each item: P10065 (7 correct), P09684 (10 correct), P11607 (21 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> 
#> 
#> =======================================
#> Model No. 2
#>     Model name:           domainreading
#>     Number of items:      137
#>     Number of persons:    4532
#>     Number of dimensions: 1
#> =======================================
#> 
#> Warning! Model No. 2, model name: 'domainreading': 12 from 149 items listed the Q matrix not found in data:
#>     T09_11, T01_05, T03_07, T10_08, T11_03, T04_01, T11_04, T04_07, T06_04, T06_03, T08_04, T08_05
#> Following 12 item(s) missed in data frame will be removed from Q matrix: 
#>     T09_11, T01_05, T03_07, T10_08, T11_03, T04_01, T11_04, T04_07, T06_04, T06_03, T08_04, T08_05
#> Found 59 cases with missings on all items.
#> Cases with missings on all items will be deleted.
#> 83 subject(s) do not solve any item:
#>    P09337 (6 false), P09202 (10 false), P10569 (17 false) ... 
#> 137 subject(s) solved each item: P09247 (6 correct), P09425 (7 correct), P12928 (24 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).

# run 2 models
runT3 <- runModel(defT3)

# get the results
resT3 <- getResults(runT3)
#> |*****|
#> |-----|
#> |*****|
#> |-----|

# collect item parameters
itemT3<- itemFromRes(resT3)

# Second step: compute linking constant. We link the items of 't3' to the items
# of 't2' which are already transformed to the metric of 't1' (we use the item
# parameters which were used for plausible values imputation at 't2' as norm
# parameters)
L.t2t3<- equat1pl ( results = resT3, prmNorm = T.t1t2[["itempars"]][,c("item", "estTransf")],
         excludeLinkingDif = TRUE, difBound = 0.64, iterativ = TRUE)
#> Found 2 model(s).
#>    Equating is executed for each dimension in each model separately.
#> Error: object 'T.t1t2' not found

# linking constant is negative: students performance at T3 is worse than T1
# Third step: transform item parameters of 't3' to the common metric of 't1' and 't2'
# We already know the 'refPop' values.
ref   <- tfRef[["refPop"]][,-2]
T.t2t3<- transformToBista ( equatingList = L.t2t3, refPop=ref, cuts = cuts,
         vera=FALSE, years = c(2015,2020))
#> Error: object 'L.t2t3' not found

# Fourth step: drawing plausible values for 't3'. We use the transformed item
# parameters (captured in 'T.t2t3') for anchoring

# Running second step: split models according to person groups (countries)
# ('all.persons' must be FALSE, otherwise the whole group would be treated as
# a separate distinct group.)
modT3P<- splitModels ( person.groups = datT3[,c("idstud", "country")] ,
         all.persons = FALSE, nCores = 1)
#> --------------------------------
#> splitModels: generating 3 models
#> ...
#> see <returned>$models
#> number of cores: 1
#> --------------------------------

# define the 3 country-specific 2-dimensional models, specifying latent regression
# model and fixed item parameters. We used the transformed item parameters (captured
# in 'T.t2t3[["itempars"]]' --- using the 'estTransf' column) for anchoring.
# Again, 'ses' is scaled (mean = 0, sd = 1)
datT3[,"ses_scaled"] <- scale(datT3[,"ses"])[,1]
defT3P<- defineModel(dat = datT3, items = itemT3[,"item"], id = "idstud",
         check.for.linking = TRUE, splittedModels = modT3P, qMatrix = qMat,
         anchor = T.t2t3[["itempars"]][,c("item", "estTransf")],
         HG.var = c("sex", "ses_scaled", "language"), software = "tam")
#> 
#> Specification of 'qMatrix' and 'person.groups' results in 3 model(s).
#> Error in eval(cl1[[u]]): object 'itemT3' not found

# run the 3 models (estimation takes approx. 20 seconds)
runT3P<- runModel(defT3P)
#> Error: object 'defT3P' not found

# get the results
resT3P<- getResults(runT3P)
#> Error: object 'runT3P' not found

# equating is not necessary, as the models run with fixed item parameters
# However, to prepare for the transformation on the 'bista' metric, run
# 'equat1pl' with empty arguments
ankT3P<- equat1pl ( results = resT3P)
#> Error: object 'resT3P' not found

# transformation to the 'bista' metric, using the previously defined cut scores
# and the reference population mean and sd from 't1'
dfrT3P<- transformToBista ( equatingList = ankT3P, refPop=ref, cuts=cuts, vera=FALSE)
#> Error: object 'ankT3P' not found


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
#> Error: object 'dfrT1P' not found

# plausibles values of measurement occasion 2 ('t2'): add the year to the data.frame
persT2<- data.frame ( year = 2015,
         dfrT2P[["personpars"]][,c("idstud", "dimension", "imp", "value", "valueTransfBista", "traitLevel")],
         stringsAsFactors = FALSE)
#> Error: object 'dfrT2P' not found

# plausibles values of measurement occasion 3 ('t3'): add the year to the data.frame
persT3<- data.frame ( year = 2020,
         dfrT3P[["personpars"]][,c("idstud", "dimension", "imp", "value", "valueTransfBista", "traitLevel")],
         stringsAsFactors = FALSE)
#> Error: object 'dfrT3P' not found

# bind together in a common data.frame
pers  <- rbind(persT1, persT2, persT3)
#> Error: object 'persT1' not found

# merge background variables to plausible values data
# first we have to create the 'domain' column in plausible values data
pers[,"domain"] <- car::recode(pers[,"dimension"], "'domainlistening'='listening'; 'domainreading'='reading'")
#> Error in `[.data.frame`(pers, , "dimension"): undefined columns selected
pers[,"dimension"] <- NULL
pers  <- eatTools::mergeAttr(unique(trends[,c("year", "idclass", "idstud", "wgt", "jkzone", "jkrep", "domain", "country", "language", "ses", "sex")]),
         pers, by = c("year", "idstud", "domain"), all = FALSE, setAttr=FALSE)
#> Error in `[.data.frame`(y, , by.y, drop = FALSE): undefined columns selected

# collect original linking errors
# t1 vs. t2: linking errors were computed in example 6a.
let1t2<- T.t1t2[["linkingErrors"]]
#> Error: object 'T.t1t2' not found

# t2 vs. t3: linking errors were computed in example 6b.
let2t3<- T.t2t3[["linkingErrors"]]
#> Error: object 'T.t2t3' not found

# t1 vs. t3: linking errors were not yet computed: link t3 to t1 to create linking error template
L.t1t3<- equat1pl ( results = resT3, prmNorm = itemT1[,c("item", "est")],
         excludeLinkingDif = TRUE, difBound = 0.64, iterativ = TRUE)
#> Found 2 model(s).
#>    Equating is executed for each dimension in each model separately.
#> 
#> ====================================================================================================
#>  
#> Model No. 1
#>     Model name:                domainlistening
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainlistening
#>     Number of linking items:   36
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainlistening': 4 of 36 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T15_02'.
#>    Iteration 2: Exclude item 'T12_05'.
#>    Iteration 3: Exclude item 'T14_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                     -0.355     0.078
#> 2 iterativ    1       T15_02       -1.746           -0.406     0.061
#> 3 iterativ    2       T12_05        1.242           -0.370     0.051
#> 4 iterativ    3       T14_04       -0.886           -0.396     0.045
#> 
#> 
#> ====================================================================================================
#>  
#> Model No. 2
#>     Model name:                domainreading
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainreading
#>     Number of linking items:   68
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainreading': 4 of 68 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T01_07'.
#>    Iteration 2: Exclude item 'T02_07'.
#>    Iteration 3: Exclude item 'T09_07'.
#>    Iteration 4: Exclude item 'T04_04'.
#>    Iteration 5: Exclude item 'T07_10'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                     -0.305     0.048
#> 2 iterativ    1       T01_07        1.656           -0.280     0.041
#> 3 iterativ    2       T02_07       -0.849           -0.293     0.040
#> 4 iterativ    3       T09_07       -0.735           -0.304     0.039
#> 5 iterativ    4       T04_04        -0.67           -0.314     0.038
#> 6 iterativ    5       T07_10       -0.644           -0.325     0.037
#> 
         
# indirect linking ('chained' linking)
chain <- multiEquatError (eq.1_2=L.t1t2, eq.2_3=L.t2t3, eq.1_3=L.t1t3)
#> Error in multiEquatError(eq.1_2 = L.t1t2, eq.2_3 = L.t2t3, eq.1_3 = L.t1t3): unused arguments (eq.1_2 = L.t1t2, eq.2_3 = L.t2t3, eq.1_3 = L.t1t3)

# replace direct linking errors with indirect linking errors
L.t1t3<- replaceLinkingError (equatingList =L.t1t3, multiEquatError_output=chain)
#> Error: object 'chain' not found

# transform linking errors
ref   <- tfRef[["refPop"]][,-2]
tle   <- transformToBista ( equatingList = L.t1t3, refPop=ref, cuts = cuts,
         vera=FALSE, years = c(2010,2020))
#> Error in transformToBista(equatingList = L.t1t3, refPop = ref, cuts = cuts,     vera = FALSE, years = c(2010, 2020)): Invalid 'refPop'.
let1t3<- tle[["linkingErrors"]]
#> Error: object 'tle' not found

# bind all linking errors in a common data.frame
lErr  <- rbind(let1t2, let2t3, let1t3)
#> Error: object 'let1t2' not found
lErr[,"domain"] <- car::recode(lErr[,"domain"], "'domainlistening'='listening'; 'domainreading'='reading'")
#> Error: object 'lErr' not found


################################################################################
###  Example 6d: eatRep estimation with plausible values and linking errors  ###
################################################################################

# Example 6d needs the objects ('pers', 'lErr') created in example 6c
# load the 'eatRep' package ... note: needs eatRep version 0.14.0 or higher
library(eatRep)
#> Warning: package 'eatRep' was built under R version 4.4.3
#> Loading required package: survey
#> Warning: package 'survey' was built under R version 4.4.3
#> Loading required package: grid
#> Loading required package: Matrix
#> Loading required package: survival
#> 
#> Attaching package: 'survey'
#> The following object is masked from 'package:graphics':
#> 
#>     dotchart
#> Loading required package: BIFIEsurvey
#> |-----------------------------------------------------------------
#> | BIFIEsurvey 3.6-6 (2024-04-25 13:32:27)
#> | http://www.bifie.at                                             
#> |-----------------------------------------------------------------
#> Loading required package: progress
#> Loading required package: lavaan
#> Warning: package 'lavaan' was built under R version 4.4.3
#> This is lavaan 0.6-20
#> lavaan is FREE software! Please report any bugs.

# compute means for both countries with trend, for both domains separately,
# using replications methods (jackknife-1)
means <- by(data = pers, INDICES = pers[,"domain"], FUN = function ( dim ) {
         m <- repMean(datL = dim, ID="idstud", PSU = "idclass", type = "jk1",
              imp = "imp", groups = "country", dependent = "valueTransfBista",
              trend = "year", linkErr = lErr[which(lErr[,"domain"] == dim[1,"domain"]),],
              engine="BIFIEsurvey")
         r <- report(m, add = list(domain = dim[1,"domain"]))
         return(r)})
#> Error in `[.data.frame`(pers, , "domain"): undefined columns selected
means <- do.call("rbind", means)
#> Error: object 'means' not found
         
# additionally: differ the sex-specific means in each country from the sex-specific means
# in the whole population? Are the differences (male vs. female) in each country different
# from the difference (male vs. female) in the whole population?
means2<- by(data = pers, INDICES = pers[,"domain"], FUN = function ( dim ) {
         m <- repMean(datL = dim, ID="idstud", PSU = "idclass", type = "jk1",
              imp = "imp", groups = c("country","sex"), group.differences.by = "sex",
              group.splits = 0:1, cross.differences = TRUE,crossDiffSE.engine= "lm",
              dependent = "valueTransfBista", trend = "year",
              linkErr = lErr[which(lErr[,"domain"] == dim[1,"domain"]),],
              engine="BIFIEsurvey")
         r <- report(m, add = list(domain = dim[1,"domain"]))
         return(r)})
#> Error in `[.data.frame`(pers, , "domain"): undefined columns selected
means2<- do.call("rbind", means2)
#> Error: object 'means2' not found


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
#> Error in `[.data.frame`(pers, , "domain"): undefined columns selected
freqs <- do.call("rbind", freqs)
#> Error: object 'freqs' not found

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
#> Error in `[.data.frame`(pers, , "domain"): undefined columns selected
freqs1<- do.call("rbind", freqs1)
#> Error: object 'freqs1' not found

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
#> Error in `[.data.frame`(pers, , "domain"): undefined columns selected
freqs2<- do.call("rbind", freqs2)
#> Error: object 'freqs2' not found


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
#> --------------------------------
#> splitModels: generating 2 models
#> ..
#> see <returned>$models
#> number of cores: 1
#> --------------------------------

# lets specify a 2pl model with constraints: a common discrimination for all
# items belonging to the same task
slopes<- data.frame ( variable = qMat[,"item"],
         slope = as.numeric(as.factor(substr(qMat[,"item"],1,3))))

# prepare 2pl model
defT1 <- defineModel(dat = datT1, id = "idstud", check.for.linking = TRUE,
         splittedModels = modsT1, irtmodel = "2PL.groups", est.slopegroups = slopes,
         software = "tam")
#> 
#> Specification of 'qMatrix' and 'person.groups' results in 2 model(s).
#> Error in eval(cl1[[u]]): object 'slopes' not found

# run 2 models
runT1 <- runModel(defT1)

# get the results
resT1 <- getResults(runT1)
#> |*****|
#> |-----|
#> |*****|
#> |-----|

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
#> --------------------------------
#> splitModels: generating 3 models
#> ...
#> see <returned>$models
#> number of cores: 1
#> --------------------------------

# define the 3 country-specific 2-dimensional models, specifying latent regression
# model and fixed item and fixed slope parameters.
defT1P<- defineModel(dat = datT1, items = itemT1[,"item"], id = "idstud", irtmodel = "2PL",
         check.for.linking = TRUE, splittedModels = modT1P, qMatrix = qMat,
         anchor = itemT1[,c("item", "est")], fixSlopeMat = itemT1[,c("item", "estSlope")],
         HG.var = c("ses", "sex", "language"), software = "tam")
#> Error in `[.data.frame`(itemT1, , c("item", "estSlope")): undefined columns selected

# run the 3 models
runT1P<- runModel(defT1P)
#> Error: object 'defT1P' not found

# get the results
resT1P<- getResults(runT1P, Q3 = FALSE)
#> Error: object 'runT1P' not found


################################################################################
###   Example 8: anchored partial credit model excluding linking DIF (TAM)   ###
################################################################################

# load partial credit long format data
data(reading)
#> Warning: data set 'reading' not found

# transform into wide format
datW <- reshape2::dcast(reading[which(reading[,"type"] != "iglu"),],idstud+sex+language+country~item, value.var = "valueSum")
#> Error: object 'reading' not found

# only some items are partial credit, which ones?
pc.it<- reading[which(reading[,"valueSum"] > 1),"item"] |> unique()
#> Error: object 'reading' not found

# give values of partial credit items
# few observation of 1-category for item D223143 ... may cause convergence trouble
lapply(datW[,pc.it], FUN = table)
#> Error: object 'pc.it' not found

# combine the two lowest categories to avoid categories with too few observations
datW[,"D205143"] <- car::recode(datW[,"D205143"], "1=0; 2=1; 3=2; 4=3; 5=4")
#> Error in `[.data.frame`(datW, , "D205143"): undefined columns selected

# partial credit model ... we consider the female subgroup to be the reference population
# (norm population) that defines the scale and reference item parameters
dFema<- datW[which(datW[,"sex"] == "female"),]                                  ### data for females
lapply(dFema[,pc.it], FUN = table)                                              ### very few observations for some categories ... may cause convergence trouble
#> Error: object 'pc.it' not found
def1 <- defineModel(dat=dFema, items = -c(1:4), id=1, irtmodel = "PCM", software="tam", nodes = 21)
#> Warning: Found 89 non-numeric values in 1 of 35 items:
#> ℹ Items: 'country'
#> ℹ Non-numeric values: 'countryC', 'countryA', 'countryB'
#> Warning: Found unexpected class type(s) in item response columns:
#> ℹ 1 item columns of class 'character': 'country'
#> ℹ All item columns will be transformed to be 'numeric'. Recommend to edit your
#>   data manually prior to analysis
#> Warning: There was 1 warning in `mutate()`.
#> ℹ In argument: `country = (function (x, maintain.factor.scores = TRUE,
#>   force.string = TRUE, ...`.
#> Caused by warning:
#> ! Variable has been coerced to numeric, NAs have been induced.
#> Warning: 1 testitem without any values:
#> ℹ Remove these items from the data set: 'country'
#> Warning: 1 testitem with less than 50 valid responses.
#> ℹ These items are nevertheless kept in the data set: 'country'
#> Warning: 1 testitem is constants. Remove these items from the data set:
#> ℹ Item 'sexNum', only value '1' occurs: 89 valid responses.
#> Remove 2 test item(s) overall.
#> Dataset is completely linked.
#> Q matrix specifies 1 dimension(s).
run1 <- runModel(def1)

# plot polytomous item "D205143": category 3 isn't necessary
ind1 <- grep("D205143", run1$item$item)                                         ### where is item "D205143"
foo1 <- capture.output(plot(run1, items = ind1, type="items", export=FALSE, low=-6, high=6))
#> Error in plot.tam.mml(run1, items = ind1, type = "items", export = FALSE,     low = -6, high = 6): object 'theta0' not found
res1 <- getResults(run1)
#> |*****|
#> |-----|
it1  <- itemFromRes(res1)

# males are focus group: initial free estimation of item parameters
def2 <- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1, irtmodel = "PCM",software="tam")
#> Warning: Found 83 non-numeric values in 1 of 35 items:
#> ℹ Items: 'country'
#> ℹ Non-numeric values: 'countryC', 'countryA', 'countryB'
#> Warning: Found unexpected class type(s) in item response columns:
#> ℹ 1 item columns of class 'character': 'country'
#> ℹ All item columns will be transformed to be 'numeric'. Recommend to edit your
#>   data manually prior to analysis
#> Warning: There was 1 warning in `mutate()`.
#> ℹ In argument: `country = (function (x, maintain.factor.scores = TRUE,
#>   force.string = TRUE, ...`.
#> Caused by warning:
#> ! Variable has been coerced to numeric, NAs have been induced.
#> Warning: 1 testitem without any values:
#> ℹ Remove these items from the data set: 'country'
#> Warning: 1 testitem with less than 50 valid responses.
#> ℹ These items are nevertheless kept in the data set: 'country'
#> Warning: 1 testitem is constants. Remove these items from the data set:
#> ℹ Item 'sexNum', only value '0' occurs: 83 valid responses.
#> Remove 2 test item(s) overall.
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
run2 <- runModel(def2)
res2 <- getResults(run2)
#> |*****|
#> |-----|
it2  <- itemFromRes(res2)

# link males to females ... males perform worse
# 10 items with linking dif identified 
eq   <- equat1pl(results = res2, prmNorm = it1, item = "item", cat="category", value = "est", difBound=.64, iterativ = TRUE)
#> Error in equat1pl(results = res2, prmNorm = it1, item = "item", cat = "category",     value = "est", difBound = 0.64, iterativ = TRUE): unused argument (cat = "category")

# re-calibrate males with anchoring: use the DIF-cleaned set of anchor parameters for calibrating males on the reference scale
# variant 1: use the original item parameters for females and exclude linking dif items (it1_cleaned)
weg  <- eq[["items"]][["not_specified"]][["Dim1"]][["info"]][-1,"itemExcluded"]
#> Error: object 'eq' not found
weg  <- eatTools::whereAre(weg, paste(it1[,"item"], it1[,"category"], sep="_"))
#> Error: object 'weg' not found
it1C <- it1[-weg,]
#> Error: object 'weg' not found
def3A<- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1, irtmodel = "PCM",
        anchor = it1C, itemCol = "item", valueCol = "est", catCol = "category", software="tam")
#> Error in defineModel(dat = datW[which(datW[, "sex"] == "male"), ], items = -c(1:4),     id = 1, irtmodel = "PCM", anchor = it1C, itemCol = "item",     valueCol = "est", catCol = "category", software = "tam"): unused argument (catCol = "category")
run3A<- runModel(def3A)
#> Error: object 'def3A' not found
res3A<- getResults(run3A)
#> Error: object 'run3A' not found
it3A <- itemFromRes(res3A)                                                      ### all items except the ones with linking dif with equal item parameters? check 
#> Error: object 'res3A' not found
comp <- merge(it1[,c("item", "category", "est")], it3A[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"))
#> Error in `[.data.frame`(it1, , c("item", "category", "est")): undefined columns selected
equal<- na.omit(comp[,c("est_ref", "offset")])
#> Error: object 'comp' not found
stopifnot(all(equal[,1] == equal[,2]))                                          ### all item parameters without linking dif should be equal
#> Error: object 'equal' not found
lDif <- subset(comp, !is.na(est_foc))                                           ### all items with specific focus parameter must be included in linking DIF exclusion list 
#> Error: object 'comp' not found
stopifnot(all(paste(lDif[,"item"], lDif[,"category"], sep="_") %in% eq$items[["not_specified"]][["Dim1"]][["info"]][,"itemExcluded"]))
#> Error: object 'lDif' not found

# variant 2: use the item parameters for males (linking dif items excluded), transformed to the metric of females
def3B<- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1, irtmodel = "PCM",
        anchor = eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]],
        itemCol = "item", valueCol = "est", catCol = "category", software="tam")
#> Error in defineModel(dat = datW[which(datW[, "sex"] == "male"), ], items = -c(1:4),     id = 1, irtmodel = "PCM", anchor = eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]],     itemCol = "item", valueCol = "est", catCol = "category",     software = "tam"): unused argument (catCol = "category")
run3B<- runModel(def3B)
#> Error: object 'def3B' not found
res3B<- getResults(run3B)
#> Error: object 'run3B' not found
it3B <- itemFromRes(res3B)                                                      ### all items except the ones with linking dif with equal item parameters? check 
#> Error: object 'res3B' not found
link <- eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]][,c("item", "category", "est")]
#> Error: object 'eq' not found
comp <- merge(link, it3B[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"))
#> Error: object 'link' not found
equal<- na.omit(comp[,c("est_ref", "offset")])
#> Error: object 'comp' not found
stopifnot(all(equal[,1] == equal[,2]))                                          ### all item parameters without linking dif should be equal
#> Error: object 'equal' not found
lDif <- subset(comp, !is.na(est_foc))                                           ### all items with specific focus paraeter must be included in linking DIF exclusion list 
#> Error: object 'comp' not found
stopifnot(all(paste(lDif[,"item"], lDif[,"category"], sep="_") %in% eq$items[["not_specified"]][["Dim1"]][["info"]][,"itemExcluded"]))
#> Error: object 'lDif' not found

# transform to Bista metric
eq4  <- equat1pl(results = res3B)                                               ### reference population mean and SD
#> Error: object 'res3B' not found
refP <- data.frame(domain = "Dim1", m = 0.0389, sd = 1.07108, stringsAsFactors = FALSE)
cuts <- list ( Dim1 = list(values = 390+0:3*75))
tf4  <- transformToBista(equatingList=eq4, refPop=refP, cuts=cuts, vera = FALSE)
#> Error: object 'eq4' not found


################################################################################
### Example 8.1: same like example 8 (anchored 1pl PCM), now using Conquest  ###
################################################################################

# load partial credit long format data
data(reading)
#> Warning: data set 'reading' not found

# transform into wide format
datW <- reshape2::dcast(reading[which(reading[,"type"] != "iglu"),],idstud+sex+language+country~item, value.var = "valueSum")
#> Error: object 'reading' not found

# only some items are partial credit, which ones?
pc.it<- reading[which(reading[,"valueSum"] > 1),"item"] |> unique()
#> Error: object 'reading' not found

# give values of partial credit items
# few observation of 1-category for item D223143 ... may cause convergence trouble
lapply(datW[,pc.it], FUN = table)
#> Error: object 'pc.it' not found

# combine the two lowest categories to avoid categories with too few observations
datW[,"D205143"] <- car::recode(datW[,"D205143"], "1=0; 2=1; 3=2; 4=3; 5=4")
#> Error in `[.data.frame`(datW, , "D205143"): undefined columns selected

# partial credit model ... we consider the female subgroup to be the reference population
# (norm population) that defines the scale and reference item parameters
dFema<- datW[which(datW[,"sex"] == "female"),]                                  ### data for females
lapply(dFema[,pc.it], FUN = table)                                              ### very few observations for some categories ... may cause convergence trouble
#> Error: object 'pc.it' not found
def1 <- defineModel(dat=dFema, items = -c(1:4), id=1, model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest_females", dir=tempdir())
#> Error in defineModel(dat = dFema, items = -c(1:4), id = 1, model.statement = "item+item*step",     nodes = 21, analysis.name = "pcm_conquest_females", dir = tempdir()): Assertion on 'conquest.folder' failed: Directory expected, but file in place: 'C:/Users/grewered/AppData/Local/Programs/R/R-4.4.2/library/eatModel/exec/console_Feb2007.exe'.
run1 <- runModel(def1)
res1 <- getResults(run1)
#> |*****|
#> |-----|
it1  <- itemFromRes(res1)

# males are focus group: initial free estimation of item parameters
def2 <- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1, model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest_males", dir=tempdir())
#> Error in defineModel(dat = datW[which(datW[, "sex"] == "male"), ], items = -c(1:4),     id = 1, model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest_males",     dir = tempdir()): Assertion on 'conquest.folder' failed: Directory expected, but file in place: 'C:/Users/grewered/AppData/Local/Programs/R/R-4.4.2/library/eatModel/exec/console_Feb2007.exe'.
run2 <- runModel(def2)
res2 <- getResults(run2)
#> |*****|
#> |-----|
it2  <- itemFromRes(res2)

# link males to females ... males perform worse
# 11 items with linking dif identified
eq   <- equat1pl(results = res2, prmNorm = it1, item = "item", cat="category", value = "est", difBound=.64, iterativ = TRUE)
#> Error in equat1pl(results = res2, prmNorm = it1, item = "item", cat = "category",     value = "est", difBound = 0.64, iterativ = TRUE): unused argument (cat = "category")

# re-calibrate males with anchoring: use the DIF-cleaned set of anchor parameters for calibrating males on the reference scale
# variant 1: use the original item parameters for females and exclude linking dif items (it1_cleaned)
# please note that to date in Conquest only dichotomous items can be used for anchoring
weg  <- eq[["items"]][["pcm_conquest_males"]][["Dim1"]][["info"]][-1,"itemExcluded"]
#> Error: object 'eq' not found
weg  <- eatTools::whereAre(weg, paste(it1[,"item"], it1[,"category"], sep="_"))
#> Error: object 'weg' not found
it1C <- it1[-weg,]
#> Error: object 'weg' not found
pcmIt<- unique(subset(it1C, category == "Cat2")[,"item"])                       ### additionally exclude pcm items
#> Error: object 'it1C' not found
it1C <- it1C[-eatTools::whereAre(pcmIt, it1C[,"item"]),]
#> Error: object 'it1C' not found
def3A<- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1,
        model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest_males_anchored1", dir=tempdir(),
        anchor = it1C[,c("item", "est")])
#> Error in defineModel(dat = datW[which(datW[, "sex"] == "male"), ], items = -c(1:4),     id = 1, model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest_males_anchored1",     dir = tempdir(), anchor = it1C[, c("item", "est")]): Assertion on 'conquest.folder' failed: Directory expected, but file in place: 'C:/Users/grewered/AppData/Local/Programs/R/R-4.4.2/library/eatModel/exec/console_Feb2007.exe'.
run3A<- runModel(def3A)
#> Error: object 'def3A' not found
res3A<- getResults(run3A)
#> Error: object 'run3A' not found
it3A <- itemFromRes(res3A)                                                      ### all dichotomous items except the ones with linking dif with equal item parameters? check
#> Error: object 'res3A' not found
comp <- merge(it1[,c("item", "category", "est")], it3A[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"))
#> Error in `[.data.frame`(it1, , c("item", "category", "est")): undefined columns selected
equal<- na.omit(comp[,c("est_ref", "offset")])
#> Error: object 'comp' not found
stopifnot(all(equal[,1] == equal[,2]))                                          ### all item parameters without linking dif should be equal
#> Error: object 'equal' not found
### TO DO: all items with specific focus parameter must be included in linking DIF exclusion list or must be pcm items

# variant 2: use the item parameters for males (linking dif items excluded), transformed to the metric of females
ank3B<- eq[["items"]][["pcm_conquest_males"]][["Dim1"]][["cleanedLinkItemPars"]]
#> Error: object 'eq' not found
pcmIt<- unique(subset(ank3B, category == "Cat2")[,"item"])                      ### additionally exclude pcm items
#> Error: object 'ank3B' not found
ank3B<- ank3B[-eatTools::whereAre(pcmIt, ank3B[,"item"]),]
#> Error: object 'ank3B' not found
def3B<- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1,
        model.statement = "item+item*step", analysis.name = "pcm_conquest_males_anchored2",
        dir=tempdir(), anchor = ank3B[,c("item", "est")])
#> Error in defineModel(dat = datW[which(datW[, "sex"] == "male"), ], items = -c(1:4),     id = 1, model.statement = "item+item*step", analysis.name = "pcm_conquest_males_anchored2",     dir = tempdir(), anchor = ank3B[, c("item", "est")]): Assertion on 'conquest.folder' failed: Directory expected, but file in place: 'C:/Users/grewered/AppData/Local/Programs/R/R-4.4.2/library/eatModel/exec/console_Feb2007.exe'.
run3B<- runModel(def3B)
#> Error: object 'def3B' not found
res3B<- getResults(run3B)
#> Error: object 'run3B' not found
it3B <- itemFromRes(res3B)                                                      ### all items except the ones with linking dif with equal item parameters? check
#> Error: object 'res3B' not found
link <- eq[["items"]][["pcm_conquest_males"]][["Dim1"]][["cleanedLinkItemPars"]][,c("item", "category", "est")]
#> Error: object 'eq' not found
comp <- merge(link, it3B[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"))
#> Error: object 'link' not found
equal<- na.omit(comp[,c("est_ref", "offset")])
#> Error: object 'comp' not found
stopifnot(all(equal[,1] == equal[,2]))                                          ### all item parameters without linking dif should be equal
#> Error: object 'equal' not found
### TO DO: all items with specific focus parameter must be included in linking DIF exclusion list or must be pcm items

# transform to Bista metric
eq4  <- equat1pl(results = res3B)                                               ### reference population mean and SD
#> Error: object 'res3B' not found
refP <- data.frame(domain = "Dim1", m = 0.0389, sd = 1.07108, stringsAsFactors = FALSE)
cuts <- list ( Dim1 = list(values = 390+0:3*75))
tf4  <- transformToBista(equatingList=eq4, refPop=refP, cuts=cuts, vera = FALSE)
#> Error: object 'eq4' not found


################################################################################
### Example 9: (restricted) generalized partial credit model with anchoring  ###
################################################################################

# load partial credit long format data
data(reading)
#> Warning: data set 'reading' not found

# transform into wide format
datW <- reshape2::dcast(reading[which(reading[,"type"] != "iglu"),],idstud+sex+language+country~item, value.var = "valueSum")
#> Error: object 'reading' not found

# combine the two lowest categories to avoid categories with too few observations
datW[,"D205143"] <- car::recode(datW[,"D205143"], "1=0; 2=1; 3=2; 4=3; 5=4")
#> Error in `[.data.frame`(datW, , "D205143"): undefined columns selected

# restricted generalized partial credit model ... we consider the female subgroup to be the reference population
# we assume that task D223 and D224 are created for students with special educational need and therefore have
# different slope. Non-SPF items should have slope = 1
dFema<- datW[which(datW[,"sex"] == "female"),]                                  ### data for females
items<- colnames(dFema)[-c(1:4)]
items<- data.frame(item = items, task = substr(items,1,4), slopeGrp = as.numeric(substr(items,1,4) %in% c("D223", "D224")), stringsAsFactors = FALSE) |> dplyr::mutate(slope = car::recode(slopeGrp, "0=1; else=NA"))
items[,"dim"] <- car::recode(substr(items[,"item"], 1,2),"'D0'='norm'; 'D2'='pilot'")
items<- data.frame ( items, model.matrix(~dim-1, data = items), stringsAsFactors = FALSE)
def1T<- defineModel(dat=dFema, items = items[,"item"], id=1, irtmodel = "GPCM.groups", software="tam", fixSlopeMat = na.omit(items[,c("item", "slope")]), est.slopegroups = items[,c("item", "slopeGrp")], nodes = 21, fac.oldxsi =0.6, increment.factor=1.05)
#> Error in match.arg(arg = irtmodel, choices = c("1PL", "2PL", "PCM", "PCM2",     "RSM", "GPCM", "2PL.groups", "GPCM.design", "3PL")): 'arg' should be one of "1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", "2PL.groups", "GPCM.design", "3PL"
run1T<- runModel(def1T)
#> Error: object 'def1T' not found
res1T<- getResults(run1T)
#> Error: object 'run1T' not found
it1T <- itemFromRes(res1T)
#> Error: object 'res1T' not found

# males are focus group: initial free estimation of item parameters
def2T<- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1, irtmodel = "GPCM.groups",software="tam", fixSlopeMat = na.omit(items[,c("item", "slope")]), est.slopegroups = items[,c("item", "slopeGrp")], nodes = 21)
#> Error in match.arg(arg = irtmodel, choices = c("1PL", "2PL", "PCM", "PCM2",     "RSM", "GPCM", "2PL.groups", "GPCM.design", "3PL")): 'arg' should be one of "1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", "2PL.groups", "GPCM.design", "3PL"
run2T<- runModel(def2T)
#> Error: object 'def2T' not found
res2T<- getResults(run2T)
#> Error: object 'run2T' not found
it2T <- itemFromRes(res2T)
#> Error: object 'res2T' not found
unique(round(it1T[,"estSlope"],3)); unique(round(it2T[,"estSlope"],3))          ### average discrimination (reg vs. spf) differs for males, but not for females
#> Error: object 'it1T' not found

# use only slope=1 items for 1pl linking
eq   <- equat1pl(results = res2T, prmNorm = it1T[which(it1T[,"estSlope"] ==1),], item = "item", cat="category", value = "est", difBound=.64, iterativ = TRUE)
#> Error in equat1pl(results = res2T, prmNorm = it1T[which(it1T[, "estSlope"] ==     1), ], item = "item", cat = "category", value = "est", difBound = 0.64,     iterativ = TRUE): unused argument (cat = "category")

# variant 1
# use the DIF-cleaned set of original anchor parameters for calibrating males on the reference scale
weg  <- eq[["items"]][["not_specified"]][["Dim1"]][["info"]][-1,"itemExcluded"]
#> Error: object 'eq' not found
weg  <- eatTools::whereAre(weg, paste(it1T[,"item"], it1T[,"category"], sep="_"))
#> Error: object 'weg' not found
it1C <- subset(it1T[-weg,], estSlope == 1)
#> Error: object 'it1T' not found
def3T<- defineModel(eatTools::na_omit_selection(dat=datW[which(datW[,"sex"] == "male"),],"language"), 
        items = -c(1:4), id=1, irtmodel = "GPCM.groups", anchor = it1C, HG.var = c("language", "country"), 
        itemCol = "item", valueCol = "est", catCol = "category", fixSlopeMat = na.omit(items[,c("item", "slope")]),
        qMatrix = items[,c("item", "dimnorm", "dimpilot")], est.slopegroups = items[,c("item", "slopeGrp")], nodes = 21, software="tam")
#> Error in defineModel(eatTools::na_omit_selection(dat = datW[which(datW[,     "sex"] == "male"), ], "language"), items = -c(1:4), id = 1,     irtmodel = "GPCM.groups", anchor = it1C, HG.var = c("language",         "country"), itemCol = "item", valueCol = "est", catCol = "category",     fixSlopeMat = na.omit(items[, c("item", "slope")]), qMatrix = items[,         c("item", "dimnorm", "dimpilot")], est.slopegroups = items[,         c("item", "slopeGrp")], nodes = 21, software = "tam"): unused argument (catCol = "category")
run3T<- runModel(def3T)
#> Error: object 'def3T' not found
res3T<- getResults(run3T)
#> Error: object 'run3T' not found
it3T <- itemFromRes(res3T)                                                      ### all items except the ones with linking dif with equal item parameters? check 
#> Error: object 'res3T' not found
comp <- merge(subset(it1T[,c("item", "category", "est", "estSlope")],estSlope == 1), it3T[,c("item", "category", "est", "offset")], by=c("item", "category"), all=TRUE, suffixes = c("_ref", "_foc"))
#> Error: object 'it1T' not found
equal<- na.omit(comp[,c("est_ref", "offset")])
#> Error: object 'comp' not found
stopifnot(all(equal[,1] == equal[,2]))                                          ### all item parameters without linking dif should be equal
#> Error: object 'equal' not found
lDif <- subset(comp, !is.na(est_foc))
#> Error: object 'comp' not found
stopifnot(all(paste(subset(lDif, !is.na(est_ref))[,"item"], subset(lDif, !is.na(est_ref))[,"category"], sep="_") %in% eq$items[["not_specified"]][["Dim1"]][["info"]][,"itemExcluded"]))
#> Error: object 'lDif' not found

# variant 2: use the 1pl item parameters for males (linking dif items excluded), transformed to the metric of females
def4T<- defineModel(eatTools::na_omit_selection(dat=datW[which(datW[,"sex"] == "male"),],"language"), 
        items = -c(1:4), id=1, irtmodel = "GPCM.groups",  anchor = eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]],
        HG.var = c("language", "country"), itemCol = "item", valueCol = "est", catCol = "category", fixSlopeMat = na.omit(items[,c("item", "slope")]),
        qMatrix = items[,c("item", "dimnorm", "dimpilot")], est.slopegroups = items[,c("item", "slopeGrp")], nodes = 21, software="tam")
#> Error in defineModel(eatTools::na_omit_selection(dat = datW[which(datW[,     "sex"] == "male"), ], "language"), items = -c(1:4), id = 1,     irtmodel = "GPCM.groups", anchor = eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]],     HG.var = c("language", "country"), itemCol = "item", valueCol = "est",     catCol = "category", fixSlopeMat = na.omit(items[, c("item",         "slope")]), qMatrix = items[, c("item", "dimnorm", "dimpilot")],     est.slopegroups = items[, c("item", "slopeGrp")], nodes = 21,     software = "tam"): unused argument (catCol = "category")
run4T<- runModel(def4T)
#> Error: object 'def4T' not found
res4T<- getResults(run4T)
#> Error: object 'run4T' not found
it4T <- itemFromRes(res4T)                                                      ### all items except the ones with linking dif with equal item parameters? check 
#> Error: object 'res4T' not found
link <- eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]][,c("item", "category", "est")]
#> Error: object 'eq' not found
comp <- merge(link, it4T[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"), all=TRUE)
#> Error: object 'link' not found
equal<- na.omit(comp[,c("est_ref", "offset")])
#> Error: object 'comp' not found
stopifnot(all(equal[,1] == equal[,2]))                                          ### all item parameters without linking dif should be equal
#> Error: object 'equal' not found
lDif <- subset(comp, !is.na(est_foc))                                           ### all items with specific focus paraeter must be included in linking DIF exclusion list 
#> Error: object 'comp' not found
stopifnot(all(paste(subset(lDif, !is.na(est_ref))[,"item"], subset(lDif, !is.na(est_ref))[,"category"], sep="_") %in% eq$items[["not_specified"]][["Dim1"]][["info"]][,"itemExcluded"]))
#> Error: object 'lDif' not found

# transform to Bista metric
eq4  <- equat1pl(results = res4T)                                               ### reference population mean and SD
#> Error: object 'res4T' not found
refP <- data.frame(domain = c("dimnorm", "dimpilot"), m = 0.0389, sd = 1.07108, stringsAsFactors = FALSE)
cuts <- list ( dimnorm = list(values = 390+0:3*75), dimpilot = list(values = 390+0:3*75))
tf4  <- transformToBista(equatingList=eq4, refPop=refP, cuts=cuts, vera = FALSE)
#> Error: object 'eq4' not found


################################################################################
###       Example 10: GPCM, similar like example 9, but now using mirt       ###
################################################################################

# load package Hmisc if necesary
if(!"Hmisc" %in% .packages()) {library(Hmisc)}
#> Warning: package 'Hmisc' was built under R version 4.4.3
#> 
#> Attaching package: 'Hmisc'
#> The following object is masked from 'package:survey':
#> 
#>     deff
#> The following objects are masked from 'package:base':
#> 
#>     format.pval, units

# load partial credit long format data
data(reading)
#> Warning: data set 'reading' not found

# transform into wide format
datW <- reshape2::dcast(reading[which(reading[,"type"] != "iglu"),],idstud+sex+language+country~item, value.var = "valueSum")
#> Error: object 'reading' not found

# only some items are partial credit, which ones?
pc.it<- reading[which(reading[,"valueSum"] > 1),"item"] |> unique()
#> Error: object 'reading' not found

# combine the two lowest categories to avoid categories with too few observations
datW[,"D205143"] <- car::recode(datW[,"D205143"], "1=0; 2=1; 3=2; 4=3; 5=4")
#> Error in `[.data.frame`(datW, , "D205143"): undefined columns selected

# restricted generalized partial credit model ... we consider the female subgroup to be the reference population
# we assume that task D223 and D224 are created for students with special educational need and therefore have
# different slope (gpcm). Non-SPF items should have slope = 1
dFema<- datW[which(datW[,"sex"] == "female"),]                                  ### data for females
items<- colnames(dFema)[-c(1:4)]
items<- data.frame(item = items, task = substr(items,1,4), slopeGrp = as.numeric(substr(items,1,4) %in% c("D223", "D224")), stringsAsFactors = FALSE) |> dplyr::mutate(irtmod = NA)

# add model information to 'item'
items[intersect(intersect(which(items[,"slopeGrp"] == 0), which(items[,"task"] %nin% c("D223", "D224"))), which(items[,"item"] %nin% pc.it)),"irtmod"] <- "Rasch"
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': object 'pc.it' not found
items[intersect(intersect(which(items[,"slopeGrp"] == 1), which(items[,"task"] %in% c("D223", "D224"))), which(items[,"item"] %nin% pc.it)),"irtmod"] <- "2PL"
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': object 'pc.it' not found

# this is partial credit which is specified in 'mirt' with 'Rasch' statement as the slope equals 1 like in dichotomous Rasch models
items[intersect(intersect(which(items[,"slopeGrp"] == 0), which(items[,"task"] %nin% c("D223", "D224"))), which(items[,"item"] %in% pc.it)),"irtmod"] <- "Rasch"
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': object 'pc.it' not found

# generalized partial credit for unidimensional models
items[intersect(which(items[,"slopeGrp"] == 1), which(items[,"item"] %in% pc.it)),"irtmod"] <- "gpcm"
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': object 'pc.it' not found

# define arbitrary dimensions
items[,"dim"] <- car::recode(substr(items[,"item"], 1,2),"'D0'='norm'; 'D2'='pilot'")
items<- data.frame ( items, model.matrix(~dim-1, data = items), stringsAsFactors = FALSE)

# procedure similar to conquest/tam
def1 <- defineModel(dat=eatTools::na_omit_selection(dFema,varsToOmitIfNA="language"), items = items[,"item"], 
        id=1,irtmodel = items[,c("item", "irtmod")], software="mirt")
#> Error in match.arg(arg = irtmodel, choices = c("1PL", "2PL", "PCM", "PCM2",     "RSM", "GPCM", "2PL.groups", "GPCM.design", "3PL")): 'arg' must be NULL or a character vector
skel1<- runModel(def1, onlySkeleton = TRUE)                                     ### check the model "skeleton"
#> Error in runModel(def1, onlySkeleton = TRUE): unused argument (onlySkeleton = TRUE)
run1 <- runModel(def1)                                                          ### run the model finally 
res1 <- getResults(run1)
#> |*****|
#> |-----|
it1  <- itemFromRes(res1)

# males are focus group: initial free estimation of item parameters
def2 <- defineModel(dat=eatTools::na_omit_selection(datW[which(datW[,"sex"] == "male"),],varsToOmitIfNA="language"), 
        items = -c(1:4), id=1, irtmodel = items[,c("item", "irtmod")],software="mirt")
#> Error in match.arg(arg = irtmodel, choices = c("1PL", "2PL", "PCM", "PCM2",     "RSM", "GPCM", "2PL.groups", "GPCM.design", "3PL")): 'arg' must be NULL or a character vector
run2 <- runModel(def2)
res2 <- getResults(run2)
#> |*****|
#> |-----|
it2  <- itemFromRes(res2)

# link males to females, using only 1pl items ... males perform worse. 10 items with linking dif identified 
eq   <- equat1pl(results = res2, prmNorm = subset(it1,estSlope ==1), item = "item", cat="category", value = "est", difBound=.64, iterativ = TRUE)
#> Error in equat1pl(results = res2, prmNorm = subset(it1, estSlope == 1),     item = "item", cat = "category", value = "est", difBound = 0.64,     iterativ = TRUE): unused argument (cat = "category")

# re-calibrate males with anchoring: use the DIF-cleaned set of anchor parameters for calibrating males on the reference scale
# variant 1: use the original (1pl) item parameters for females and exclude linking dif items (it1_cleaned)
weg  <- eq[["items"]][["not_specified"]][["Dim1"]][["info"]][-1,"itemExcluded"]
#> Error: object 'eq' not found
weg  <- eatTools::whereAre(weg, paste(it1[,"item"], it1[,"category"], sep="_"))
#> Error: object 'weg' not found
it1C <- subset(it1[-weg,], estSlope == 1)
#> Error: object 'weg' not found
def3A<- defineModel(dat=eatTools::na_omit_selection(datW[which(datW[,"sex"] == "male"),],varsToOmitIfNA="language"), items = -c(1:4), id=1,  irtmodel = items[,c("item", "irtmod")],
        anchor = it1C, itemCol = "item", valueCol = "est", catCol = "category", software="mirt")
#> Error in defineModel(dat = eatTools::na_omit_selection(datW[which(datW[,     "sex"] == "male"), ], varsToOmitIfNA = "language"), items = -c(1:4),     id = 1, irtmodel = items[, c("item", "irtmod")], anchor = it1C,     itemCol = "item", valueCol = "est", catCol = "category",     software = "mirt"): unused argument (catCol = "category")
run3A<- runModel(def3A)
#> Error: object 'def3A' not found
res3A<- getResults(run3A)
#> Error: object 'run3A' not found
it3A <- itemFromRes(res3A)                                                      ### all items except the ones with linking dif with equal item parameters? check 
#> Error: object 'res3A' not found
comp <- merge(subset(it1, estSlope==1)[,c("item", "category", "est")], subset(it3A, !is.na(offset))[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"))
#> Error in eval(e, x, parent.frame()): object 'estSlope' not found
stopifnot(all(round(comp[,"est_ref"],5) == round(comp[,"offset"],5)))           ### all item parameters without linking dif should be equal
#> Error: object 'comp' not found
# all items with specific focus parameter must be included in linking DIF exclusion list 
stopifnot(all((paste(comp[,"item"], comp[,"category"],sep="_") %in% eq$items[["not_specified"]][["Dim1"]][["info"]][-1,"itemExcluded"]) == FALSE))
#> Error: object 'comp' not found

# variant 2: use the item parameters for males (linking dif items excluded), transformed to the metric of females
def3B<- defineModel(dat=eatTools::na_omit_selection(datW[which(datW[,"sex"] == "male"),],varsToOmitIfNA="language"), items = -c(1:4), id=1,  irtmodel = items[,c("item", "irtmod")],
        anchor = eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]],
        itemCol = "item", valueCol = "est", catCol = "category", software="mirt")
#> Error in defineModel(dat = eatTools::na_omit_selection(datW[which(datW[,     "sex"] == "male"), ], varsToOmitIfNA = "language"), items = -c(1:4),     id = 1, irtmodel = items[, c("item", "irtmod")], anchor = eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]],     itemCol = "item", valueCol = "est", catCol = "category",     software = "mirt"): unused argument (catCol = "category")
run3B<- runModel(def3B)
#> Error: object 'def3B' not found
res3B<- getResults(run3B)
#> Error: object 'run3B' not found
it3B <- itemFromRes(res3B)                                                      ### all items except the ones with linking dif with equal item parameters? check 
#> Error: object 'res3B' not found
link <- eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]][,c("item", "category", "est", "estSlope")]
#> Error: object 'eq' not found
stopifnot(all(link[,"estSlope"] == 1))
#> Error: object 'link' not found
comp <- merge(link, it3B[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"))
#> Error: object 'link' not found
stopifnot(all(round(comp[,"est_ref"],5) == round(comp[,"offset"],5)))           ### all item parameters without linking dif should be equal
#> Error: object 'comp' not found
# all items with specific focus parameter must be included in linking DIF exclusion list 
stopifnot(all((paste(comp[,"item"], comp[,"category"],sep="_") %in% eq$items[["not_specified"]][["Dim1"]][["info"]][-1,"itemExcluded"]) == FALSE))
#> Error: object 'comp' not found

# transform to Bista metric
eq4  <- equat1pl(results = res3B)                                               ### reference population mean and SD
#> Error: object 'res3B' not found
refP <- data.frame(domain = "Dim1", m = 0.0389, sd = 1.07108, stringsAsFactors = FALSE)
cuts <- list ( Dim1 = list(values = 390+0:3*75))
tf4  <- transformToBista(equatingList=eq4, refPop=refP, cuts=cuts, vera = FALSE)
#> Error: object 'eq4' not found


################################################################################
###    Example 11: differential item functioning in partial credit models    ###
################################################################################

# load partial credit long format data
data(reading)
#> Warning: data set 'reading' not found

# create a small wide format exemplary data set, only booklet 08
dw  <- reshape2::dcast(subset(reading, bookletID == "TH08"), idstud+sex~item, value.var = "valueSum") |>
       dplyr::mutate(sex = car::recode(sex, "'male'=0; 'female'=1", as.factor=FALSE)) |>
       eatTools::na_omit_selection("sex")
#> Error: object 'reading' not found

# define the model: the specification resembles exemple 9 in tam.mml help page 
defT<- defineModel(dat = dw, items = -c(1:2), DIF.var="sex", id = "idstud",  irtmodel = "PCM",software="tam")
#> Error: object 'dw' not found
runT<- runModel(defT)
#> Error: object 'defT' not found
resT<- getResults(runT)
#> Error: object 'runT' not found

# the same model in conquest (recommended)
defC<- defineModel(dat = dw, items = -c(1:2), model.statement = "item+item*step - sex + item*sex", DIF.var="sex", id = "idstud", analysis.name = "dif_pcm", dir=tempdir())
#> Error: object 'dw' not found
runC<- runModel(defC, wait=FALSE)
#> Error: object 'defC' not found
resC<- getResults(runC)
#> Error: object 'runC' not found
it  <- itemFromRes(resC)
#> Error: object 'resC' not found
shw <- get.shw(file.path(tempdir(), "dif_pcm.shw"), dif.term = "item*sex")
#> Error in get.shw(file.path(tempdir(), "dif_pcm.shw"), dif.term = "item*sex"): Assertion on 'file' failed: File does not exist: 'C:\Users\grewered\AppData\Local\Temp\RtmpoLrBif/dif_pcm.shw'.


################################################################################
###    Example 12: Norm study within the realm of partial credit models      ###
################################################################################

# The following example mimics the steps the required to norm the educational standards
# when the underlying model is a partial credit model.

# load partial credit long format data. This large scale assessment data surveyed in
# an incomplete block design is representative of a domain covered by the new standards.
# We assume thta the sample was drawn from the norm population (reference population)
data(reading)
#> Warning: data set 'reading' not found

# transform the data to the wide format
datW <- reshape2::dcast(reading[which(reading[,"type"] != "iglu"),],idstud+sex+language+country~item, value.var = "valueSum")
#> Error: object 'reading' not found

# combine the two lowest categories to avoid categories with too few observations
datW[,"D205143"] <- car::recode(datW[,"D205143"], "1=0; 2=1; 3=2; 4=3; 5=4")
#> Error in `[.data.frame`(datW, , "D205143"): undefined columns selected

# partial credit model. We draw 15 plausible values to achieve greater precision
# in estimating the distribution of individuals.
def1 <- defineModel(dat=datW, items = -c(1:4), id=1, irtmodel = "PCM", software="tam", nodes = 21, n.plausible=15)
#> Warning: Found 172 non-numeric values in 1 of 35 items:
#> ℹ Items: 'country'
#> ℹ Non-numeric values: 'countryC', 'countryA', 'countryB'
#> Warning: Found unexpected class type(s) in item response columns:
#> ℹ 1 item columns of class 'character': 'country'
#> ℹ All item columns will be transformed to be 'numeric'. Recommend to edit your
#>   data manually prior to analysis
#> Warning: There was 1 warning in `mutate()`.
#> ℹ In argument: `country = (function (x, maintain.factor.scores = TRUE,
#>   force.string = TRUE, ...`.
#> Caused by warning:
#> ! Variable has been coerced to numeric, NAs have been induced.
#> Warning: 1 testitem without any values:
#> ℹ Remove these items from the data set: 'country'
#> Warning: 1 testitem with less than 50 valid responses.
#> ℹ These items are nevertheless kept in the data set: 'country'
#> Remove 1 test item(s) overall.
#> Dataset is completely linked.
#> Q matrix specifies 1 dimension(s).
run1 <- runModel(def1)

# use ntheta = 40000 for higher precision
res1 <- getResults(run1, ntheta = 40000, theta.model = FALSE)
#> |***************|
#> |---------------|

# estimate mean and standard deviation of the reference population from the plausible values
# use the long format by toWideFormat = FALSE. Load the eatRep package for mean/SD computation
library(eatRep)
pv1  <- pvFromRes(res1, toWideFormat = FALSE)
mSD  <- repMean(datL = pv1, ID="idstud", imp="imp", dependent="value")
#> 1 analyse(s) overall according to: 'group.splits = 0'.
#> Assume unnested structure with 15 imputations.
#> 
mSDR <- report2(mSD, round = FALSE)[["plain"]]; dim(mSDR)
#> [1]  3 10

# use mean/SD to transform the item and person parameters. Since the cur scores have not yet
# been established - as they are still being defined in the standard-setting process - we are
# assigning arbitrary values for now. Since the sample can be considered as being drawn from
# the reference population, no equating procedure is necessary. The function is therefore
# simply passed through. The mean and standard deviation of the educational standards metric
# are set to 500 and 100, respectively.
eq   <- equat1pl(res1)
#> Found 1 model(s).
#>    Equating is executed for each dimension in each model separately.
#> No norm parameter defined ('prmNorm' is missing). Treat current sample as drawn from the reference population.
trnsf<- transformToBista(equatingList=eq, refPop= data.frame(domain="Dim1",
        m = mSDR[which(mSDR[,"parameter"] == "mean"),"est"],
        sd = mSDR[which(mSDR[,"parameter"] == "sd"),"est"]),
        cuts = list(Dim1 = list(values = c(400,600))), vera=FALSE)
#> The 'refPop' data.frame does not include information about reference population mean/SD on Bista metric. Values will be defaulted to 500/100.
#> mutmasslicher fehler.
```
