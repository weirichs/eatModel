# 1pl linking errors for three measurement occasions

Function combines linking errors based on the 1pl linking procedure
according to
[`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md)
for three measurement occasions.

## Usage

``` r
multiEquatError (eq.1_2, eq.2_3, eq.1_3, dependentDIF =FALSE, verbose=TRUE)
```

## Arguments

- eq.1_2:

  The object returned by
  [`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md)
  for linking the first and second measurement occasion.

- eq.2_3:

  The object returned by
  [`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md)
  for linking the second and third measurement occasion.

- eq.1_3:

  The object returned by
  [`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md)
  for linking the first and third measurement occasion.

- dependentDIF:

  Logical. If DIF is assumed to be correlated between measurement
  occasions one and two versus two and three, this covariance will be
  taken into account. If the assumption of linear dependency does not
  seem reasonable, this option possibly results in overfitting and
  underestimation of the true linking error.

- verbose:

  Logical: print linking informations to console?

## Value

A list of data.frames with linking errors, according to the number of
domains. The output is intended to be used as argument of the
replaceLinkingError function.

## Author

Sebastian Weirich

## References

Monseur, C., & Berezner, A. (2007). The Computation of Equating Errors
in International Surveys in Education. *Journal of Applied Measurement
8*(3): 323-335.

## Examples

``` r
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
#> 17 subject(s) do not solve any item:
#>    P00173 (6 false), P01137 (7 false), P02863 (23 false) ... 
#> 97 subject(s) solved each item: P00078 (6 correct), P00389 (10 correct), P03338 (27 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Error in prepareDatasets(namen.all.hg = cbc[["namen.all.hg"]], dat = cpsc[["dat"]],     software = software, allNam = cbc[["allNam"]], use.letters = use.letters): object 'var.char' not found
items   <- lapply(results, itemFromRes)
#> Error: object 'results' not found
eq.1_2  <- equat1pl(items[[1]][,c("item", "est")], items[[2]][,c("item", "est")],difBound = 0.64, iterativ = TRUE)
#> Error: object 'items' not found
eq.2_3  <- equat1pl(items[[2]][,c("item", "est")], items[[3]][,c("item", "est")],difBound = 0.64, iterativ = TRUE)
#> Error: object 'items' not found
eq.1_3  <- equat1pl(items[[1]][,c("item", "est")], items[[3]][,c("item", "est")],difBound = 0.64, iterativ = TRUE)
#> Error: object 'items' not found
lErrors <- multiEquatError(eq.1_2, eq.2_3, eq.1_3)
#> Error: object 'eq.1_2' not found

# dependent DIF
lErrDif <- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)
#> Error: object 'eq.1_2' not found

# Following code demonstrates the replacement of old linking error in the 'equatingList' object
# for further processing
# direct linking 1 to 3 (1 is reference)
eq1.vs.3<- equat1pl(results[[3]], prmNorm = items[[1]][,c("item", "est")], difBound = 0.64, iterativ = TRUE)
#> Error: object 'results' not found

# replace 'direct' linking error with 'indirect' linking error from 'multiEquatError()'
eq1.vs.3<- replaceLinkingError (equatingList=eq1.vs.3, multiEquatError_output=lErrors)
#> Error: object 'lErrors' not found


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
#> --------------------------------
#> splitModels: generating 2 models
#> ..
#> see <returned>$models
#> number of cores: 1
#> --------------------------------
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
#> 76 subject(s) do not solve any item:
#>    P00009 (6 false), P00864 (6 false), P02759 (21 false) ... 
#> 55 subject(s) solved each item: P00099 (10 correct), P00404 (10 correct), P02659 (11 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Error in prepareDatasets(namen.all.hg = cbc[["namen.all.hg"]], dat = cpsc[["dat"]],     software = software, allNam = cbc[["allNam"]], use.letters = use.letters): object 'var.char' not found
eq.1_2  <- equat1pl(results2[[1]], itemFromRes(results2[[2]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE)
#> Error: object 'results2' not found
eq.2_3  <- equat1pl(results2[[2]], itemFromRes(results2[[3]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE)
#> Error: object 'results2' not found
eq.1_3  <- equat1pl(results2[[1]], itemFromRes(results2[[3]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE)
#> Error: object 'results2' not found
lErrors2<- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)
#> Error: object 'eq.1_2' not found

# direct linking 1 to 3 (1 is reference)
it1     <- itemFromRes(results2[[1]])
#> Error: object 'results2' not found
eq1.vs.3<- equat1pl(results2[[3]], prmNorm = it1[,c("item", "est")], difBound = 0.64, iterativ = TRUE)
#> Error: object 'results2' not found

# replace 'direct' linking error with 'indirect' linking error from 'multiEquatError()'
eq1.vs.3<- replaceLinkingError (equatingList=eq1.vs.3, multiEquatError_output=lErrors2)
#> Error: object 'lErrors2' not found


################################################################################
###    Example 3: Jackknife linking for three measurement and two domains    ###
################################################################################

# we reuse the previously created object 'results2' here
eq.1_2  <- equat1pl(results2[[1]], itemFromRes(results2[[2]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE)
#> Error: object 'results2' not found
eq.2_3  <- equat1pl(results2[[2]], itemFromRes(results2[[3]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE)
#> Error: object 'results2' not found
eq.1_3  <- equat1pl(results2[[1]], itemFromRes(results2[[3]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE)
#> Error: object 'results2' not found
lErrors3<- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)
#> Error: object 'eq.1_2' not found


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
#> 
#> ====================================================================================================
#>  
#> Model No. NA
#>     Model name:                Dim1
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   Dim1
#>     Number of linking items:   3
#>     Linking method:            Mean.Mean
#> 
#>   linking.constant linkerror
#> 1           -0.144     0.138
#> 
eq.2_3 <- equat1pl(e2, e3,difBound = 0.64, iterativ = TRUE)
#> 
#> ====================================================================================================
#>  
#> Model No. NA
#>     Model name:                Dim1
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   Dim1
#>     Number of linking items:   5
#>     Linking method:            Mean.Mean
#> 
#>   linking.constant linkerror
#> 1           -0.072     0.117
#> 
eq.1_3 <- equat1pl(e1, e3,difBound = 0.64, iterativ = TRUE)
#> 
#> ====================================================================================================
#>  
#> Model No. NA
#>     Model name:                Dim1
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   Dim1
#>     Number of linking items:   4
#>     Linking method:            Mean.Mean
#> 
#>   linking.constant linkerror
#> 1           -0.149     0.079
#> 
multiEquatError(eq.1_2, eq.2_3, eq.1_3)
#> 
#> Dimension 'Dim1': Direct linking error of the three combinations of measurements: 
#>         mzp1 mzp2 N.Items    SD   Var linkerror chained mzp_1.vs.3 mzp_1.vs.2.vs.3
#>       1    1    2       3 0.238 0.057     0.138      NA         NA              NA
#>       2    1    3       4 0.158 0.025     0.079    0.18     -0.149          -0.216
#>       3    2    3       5 0.261 0.068     0.117      NA         NA              NA
#> $Dim1
#>      trend31  trend3221       le31    le3221
#> 1 -0.1488454 -0.2160025 0.07918052 0.1804314
#> 
```
