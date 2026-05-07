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
#> Q matrix specifies 1 dimension(s).
#> Getting standard errors with the tam.se function: 0.4 secs
#> Getting WLEs calling tam.wle from getTamWles: 0.7 secs
#> |*****|
#> |-----|
#> Getting PVs calling tam.pv from getTamPVs: 1.7 secs
#> Getting Q3 statistic calling tam.modelfit from getTamQ3: 0.5 secs
#> 56 subject(s) do not solve any item:
#>    P04511 (6 false), P05616 (10 false), P08143 (14 false) ... 
#> 114 subject(s) solved each item: P04509 (6 correct), P08241 (8 correct), P05805 (20 correct) ... 
#> W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).
#>                   481 cases unconnected. Following items are unconnected: 
#>                   T07_01, T07_02, T07_03, T07_04, T07_05, T07_06, T07_07, T07_08, T07_09, T07_10
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> Getting standard errors with the tam.se function: 0.6 secs
#> Getting infit parameters calling tam.fit from getTamInfit: 0.3 secs
#> Getting WLEs calling tam.wle from getTamWles: 0.9 secs
#> |*****|
#> |-----|
#> Getting PVs calling tam.pv from getTamPVs: 1.7 secs
#> Getting Q3 statistic calling tam.modelfit from getTamQ3: 1.2 secs
#> 83 subject(s) do not solve any item:
#>    P09337 (6 false), P09202 (10 false), P10569 (17 false) ... 
#> 137 subject(s) solved each item: P09247 (6 correct), P09425 (7 correct), P12928 (24 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> Getting standard errors with the tam.se function: 0.8 secs
#> Getting infit parameters calling tam.fit from getTamInfit: 0.3 secs
#> Getting WLEs calling tam.wle from getTamWles: 1.3 secs
#> |*****|
#> |-----|
#> Getting PVs calling tam.pv from getTamPVs: 1.7 secs
#> Getting Q3 statistic calling tam.modelfit from getTamQ3: 1.4 secs
items   <- lapply(results, itemFromRes)
eq.1_2  <- equat1pl(items[[1]][,c("item", "est")], items[[2]][,c("item", "est")],difBound = 0.64, iterativ = TRUE)
#> 
#> ====================================================================================================
#>  
#> Model No. NA
#>     Model name:                Dim1
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   Dim1
#>     Number of linking items:   68
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'Dim1': 4 of 68 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T02_07'.
#>    Iteration 2: Exclude item 'T01_07'.
#>    Iteration 3: Exclude item 'T02_04'.
#>    Iteration 4: Exclude item 'T05_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.232     0.061
#> 2 iterativ    1       T02_07       -2.505            0.195     0.050
#> 3 iterativ    2       T01_07       -2.267            0.160     0.036
#> 4 iterativ    3       T02_04       -1.176            0.142     0.032
#> 5 iterativ    4       T05_04       -0.867            0.129     0.029
#> 
eq.2_3  <- equat1pl(items[[2]][,c("item", "est")], items[[3]][,c("item", "est")],difBound = 0.64, iterativ = TRUE)
#> 
#> ====================================================================================================
#>  
#> Model No. NA
#>     Model name:                Dim1
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   Dim1
#>     Number of linking items:   119
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'Dim1': 7 of 119 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T02_07'.
#>    Iteration 2: Exclude item 'T05_01'.
#>    Iteration 3: Exclude item 'T23_14'.
#>    Iteration 4: Exclude item 'T04_08'.
#>    Iteration 5: Exclude item 'T03_03'.
#>    Iteration 6: Exclude item 'T10_10'.
#>    Iteration 7: Exclude item 'T02_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.098     0.040
#> 2 iterativ    1       T02_07        3.381            0.126     0.029
#> 3 iterativ    2       T05_01       -1.045            0.117     0.028
#> 4 iterativ    3       T23_14       -0.992            0.109     0.026
#> 5 iterativ    4       T04_08        0.915            0.117     0.025
#> 6 iterativ    5       T03_03       -0.687            0.111     0.025
#> 7 iterativ    6       T10_10         0.68            0.117     0.024
#> 8 iterativ    7       T02_04        0.654            0.123     0.024
#> 
eq.1_3  <- equat1pl(items[[1]][,c("item", "est")], items[[3]][,c("item", "est")],difBound = 0.64, iterativ = TRUE)
#> 
#> ====================================================================================================
#>  
#> Model No. NA
#>     Model name:                Dim1
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   Dim1
#>     Number of linking items:   68
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'Dim1': 4 of 68 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T01_07'.
#>    Iteration 2: Exclude item 'T02_07'.
#>    Iteration 3: Exclude item 'T09_07'.
#>    Iteration 4: Exclude item 'T04_04'.
#>    Iteration 5: Exclude item 'T07_10'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.305     0.048
#> 2 iterativ    1       T01_07       -1.656            0.280     0.041
#> 3 iterativ    2       T02_07        0.849            0.293     0.040
#> 4 iterativ    3       T09_07        0.734            0.304     0.039
#> 5 iterativ    4       T04_04         0.67            0.314     0.038
#> 6 iterativ    5       T07_10        0.644            0.325     0.037
#> 
lErrors <- multiEquatError(eq.1_2, eq.2_3, eq.1_3)
#> 
#> Dimension 'Dim1': Direct linking error of the three combinations of measurements: 
#>         mzp1 mzp2 N.Items    SD   Var linkerror chained mzp_1.vs.3 mzp_1.vs.2.vs.3
#>       1    1    2      64 0.234 0.055     0.029      NA         NA              NA
#>       2    1    3      63 0.295 0.087     0.037   0.038      0.325           0.251
#>       3    2    3     112 0.253 0.064     0.024      NA         NA              NA

# dependent DIF
lErrDif <- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)
#> Dimension 'Dim1': Negative covariance of -0.013 between DIF for measurement 1 vs. 2, and DIF for measurement 2 vs. 3. Covariance is set to 0 to avoid underestimation of the true linking error.
#> 
#> Dimension 'Dim1': Direct linking error of the three combinations of measurements: 
#>         mzp1 mzp2 N.Items    SD   Var linkerror chained mzp_1.vs.3 mzp_1.vs.2.vs.3
#>       1    1    2      64 0.234 0.055     0.029      NA         NA              NA
#>       2    1    3      63 0.295 0.087     0.037   0.038      0.325           0.251
#>       3    2    3     112 0.253 0.064     0.024      NA         NA              NA

# Following code demonstrates the replacement of old linking error in the 'equatingList' object
# for further processing
# direct linking 1 to 3 (1 is reference)
eq1.vs.3<- equat1pl(results[[3]], prmNorm = items[[1]][,c("item", "est")], difBound = 0.64, iterativ = TRUE)
#> Found 1 model(s).
#>    Equating is executed for each dimension in each model separately.
#> 
#> ====================================================================================================
#>  
#> Model No. 1
#>     Model name:                not_specified
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   Dim1
#>     Number of linking items:   68
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'Dim1': 4 of 68 items with linking |DIF| > 0.64 identified.
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

# replace 'direct' linking error with 'indirect' linking error from 'multiEquatError()'
eq1.vs.3<- replaceLinkingError (equatingList=eq1.vs.3, multiEquatError_output=lErrors)
#> Dimension 'Dim1': Replace old linking error 0.0372 with 0.0378


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
#> 17 subject(s) do not solve any item:
#>    P00173 (6 false), P01137 (7 false), P02863 (23 false) ... 
#> 97 subject(s) solved each item: P00078 (6 correct), P00389 (10 correct), P03338 (27 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> Getting standard errors with the tam.se function: 0.3 secs
#> Getting WLEs calling tam.wle from getTamWles: 0.3 secs
#> |*****|
#> |-----|
#> Getting PVs calling tam.pv from getTamPVs: 1.6 secs
#> Getting Q3 statistic calling tam.modelfit from getTamQ3: 0.2 secs
#> Getting standard errors with the tam.se function: 0.4 secs
#> Getting WLEs calling tam.wle from getTamWles: 0.6 secs
#> |*****|
#> |-----|
#> Getting PVs calling tam.pv from getTamPVs: 1.7 secs
#> Getting Q3 statistic calling tam.modelfit from getTamQ3: 0.5 secs
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
#>     Number of items:      96
#>     Number of persons:    4516
#>     Number of dimensions: 1
#> =========================================
#> 
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
#> 56 subject(s) do not solve any item:
#>    P04511 (6 false), P05616 (10 false), P08143 (14 false) ... 
#> 114 subject(s) solved each item: P04509 (6 correct), P08241 (8 correct), P05805 (20 correct) ... 
#> W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).
#>                   481 cases unconnected. Following items are unconnected: 
#>                   T07_01, T07_09, T07_07, T07_02, T07_06, T07_08, T07_05, T07_03, T07_10, T07_04
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> Getting standard errors with the tam.se function: 0.5 secs
#> Getting infit parameters calling tam.fit from getTamInfit: 0.3 secs
#> Getting WLEs calling tam.wle from getTamWles: 0.7 secs
#> |*****|
#> |-----|
#> Getting PVs calling tam.pv from getTamPVs: 1.7 secs
#> Getting Q3 statistic calling tam.modelfit from getTamQ3: 0.7 secs
#> Getting standard errors with the tam.se function: 0.7 secs
#> Getting infit parameters calling tam.fit from getTamInfit: 0.3 secs
#> Getting WLEs calling tam.wle from getTamWles: 1 secs
#> |*****|
#> |-----|
#> Getting PVs calling tam.pv from getTamPVs: 1.7 secs
#> Getting Q3 statistic calling tam.modelfit from getTamQ3: 1.1 secs
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
#>     Number of items:      119
#>     Number of persons:    4532
#>     Number of dimensions: 1
#> =========================================
#> 
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
#> Found 59 cases with missings on all items.
#> Cases with missings on all items will be deleted.
#> 83 subject(s) do not solve any item:
#>    P09337 (6 false), P09202 (10 false), P10569 (17 false) ... 
#> 137 subject(s) solved each item: P09247 (6 correct), P09425 (7 correct), P12928 (24 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 1 dimension(s).
#> Getting standard errors with the tam.se function: 0.7 secs
#> Getting infit parameters calling tam.fit from getTamInfit: 0.4 secs
#> Getting WLEs calling tam.wle from getTamWles: 0.7 secs
#> |*****|
#> |-----|
#> Getting PVs calling tam.pv from getTamPVs: 1.9 secs
#> Getting Q3 statistic calling tam.modelfit from getTamQ3: 1.1 secs
#> Getting standard errors with the tam.se function: 0.9 secs
#> Getting infit parameters calling tam.fit from getTamInfit: 0.5 secs
#> Getting WLEs calling tam.wle from getTamWles: 1.4 secs
#> |*****|
#> |-----|
#> Getting PVs calling tam.pv from getTamPVs: 1.9 secs
#> Getting Q3 statistic calling tam.modelfit from getTamQ3: 1.5 secs
eq.1_2  <- equat1pl(results2[[1]], itemFromRes(results2[[2]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE)
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
#> 1 iterativ    0                                      0.243     0.082
#> 2 iterativ    1       T15_02        1.872            0.298     0.063
#> 3 iterativ    2       T14_04       -1.674            0.250     0.042
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
#> 1 iterativ    0                                      0.232     0.061
#> 2 iterativ    1       T02_07       -2.505            0.195     0.050
#> 3 iterativ    2       T01_07       -2.267            0.160     0.036
#> 4 iterativ    3       T02_04       -1.176            0.142     0.032
#> 5 iterativ    4       T05_04       -0.867            0.129     0.029
#> 
eq.2_3  <- equat1pl(results2[[2]], itemFromRes(results2[[3]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE)
#> Found 2 model(s).
#>    Equating is executed for each dimension in each model separately.
#> 
#> ====================================================================================================
#>  
#> Model No. 1
#>     Model name:                domainlistening
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainlistening
#>     Number of linking items:   96
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainlistening': 6 of 96 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T19_11'.
#>    Iteration 2: Exclude item 'T14_04'.
#>    Iteration 3: Exclude item 'T20_30'.
#>    Iteration 4: Exclude item 'T15_18'.
#>    Iteration 5: Exclude item 'T20_03'.
#>    Iteration 6: Exclude item 'T12_05'.
#>    Iteration 7: Exclude item 'T20_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.052     0.052
#> 2 iterativ    1       T19_11        2.534            0.078     0.046
#> 3 iterativ    2       T14_04         2.53            0.105     0.037
#> 4 iterativ    3       T20_30        1.749            0.124     0.033
#> 5 iterativ    4       T15_18        1.435            0.139     0.029
#> 6 iterativ    5       T20_03        0.704            0.147     0.028
#> 7 iterativ    6       T12_05       -0.678            0.140     0.028
#> 8 iterativ    7       T20_04        0.654            0.147     0.027
#> 
#> 
#> ====================================================================================================
#>  
#> Model No. 2
#>     Model name:                domainreading
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainreading
#>     Number of linking items:   119
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainreading': 7 of 119 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T02_07'.
#>    Iteration 2: Exclude item 'T05_01'.
#>    Iteration 3: Exclude item 'T23_14'.
#>    Iteration 4: Exclude item 'T04_08'.
#>    Iteration 5: Exclude item 'T03_03'.
#>    Iteration 6: Exclude item 'T10_10'.
#>    Iteration 7: Exclude item 'T02_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.098     0.040
#> 2 iterativ    1       T02_07        3.381            0.126     0.029
#> 3 iterativ    2       T05_01       -1.045            0.117     0.028
#> 4 iterativ    3       T23_14       -0.992            0.109     0.026
#> 5 iterativ    4       T04_08        0.915            0.117     0.025
#> 6 iterativ    5       T03_03       -0.687            0.111     0.025
#> 7 iterativ    6       T10_10         0.68            0.117     0.024
#> 8 iterativ    7       T02_04        0.654            0.123     0.024
#> 
eq.1_3  <- equat1pl(results2[[1]], itemFromRes(results2[[3]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE)
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
#> 1 iterativ    0                                      0.355     0.078
#> 2 iterativ    1       T15_02        1.744            0.406     0.061
#> 3 iterativ    2       T12_05       -1.244            0.370     0.051
#> 4 iterativ    3       T14_04        0.885            0.396     0.045
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
#> 1 iterativ    0                                      0.305     0.048
#> 2 iterativ    1       T01_07       -1.656            0.280     0.041
#> 3 iterativ    2       T02_07        0.849            0.293     0.040
#> 4 iterativ    3       T09_07        0.734            0.304     0.039
#> 5 iterativ    4       T04_04         0.67            0.314     0.038
#> 6 iterativ    5       T07_10        0.644            0.325     0.037
#> 
lErrors2<- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)
#> Dimension 'domainlistening': Negative covariance of -0.021 between DIF for measurement 1 vs. 2, and DIF for measurement 2 vs. 3. Covariance is set to 0 to avoid underestimation of the true linking error.
#> 
#> Dimension 'domainlistening': Direct linking error of the three combinations of measurements: 
#>         mzp1 mzp2 N.Items    SD   Var linkerror chained mzp_1.vs.3 mzp_1.vs.2.vs.3
#>       1    1    2      34 0.246 0.060     0.042      NA         NA              NA
#>       2    1    3      33 0.256 0.065     0.045    0.05      0.396           0.397
#>       3    2    3      89 0.253 0.064     0.027      NA         NA              NA
#> Dimension 'domainreading': Negative covariance of -0.013 between DIF for measurement 1 vs. 2, and DIF for measurement 2 vs. 3. Covariance is set to 0 to avoid underestimation of the true linking error.
#> 
#> Dimension 'domainreading': Direct linking error of the three combinations of measurements: 
#>         mzp1 mzp2 N.Items    SD   Var linkerror chained mzp_1.vs.3 mzp_1.vs.2.vs.3
#>       1    1    2      64 0.234 0.055     0.029      NA         NA              NA
#>       2    1    3      63 0.295 0.087     0.037   0.038      0.325           0.251
#>       3    2    3     112 0.253 0.064     0.024      NA         NA              NA

# direct linking 1 to 3 (1 is reference)
it1     <- itemFromRes(results2[[1]])
eq1.vs.3<- equat1pl(results2[[3]], prmNorm = it1[,c("item", "est")], difBound = 0.64, iterativ = TRUE)
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

# replace 'direct' linking error with 'indirect' linking error from 'multiEquatError()'
eq1.vs.3<- replaceLinkingError (equatingList=eq1.vs.3, multiEquatError_output=lErrors2)
#> Dimension 'domainlistening': Replace old linking error 0.0445 with 0.0499
#> Dimension 'domainreading': Replace old linking error 0.0372 with 0.0378


################################################################################
###    Example 3: Jackknife linking for three measurement and two domains    ###
################################################################################

# we reuse the previously created object 'results2' here
eq.1_2  <- equat1pl(results2[[1]], itemFromRes(results2[[2]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE)
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
#>     Number of testlets:        4
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainlistening': 2 of 36 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T15_02'.
#>    Iteration 2: Exclude item 'T14_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.243     0.113
#> 2 iterativ    1       T15_02        1.872            0.298     0.054
#> 3 iterativ    2       T14_04       -1.674            0.250     0.052
#> 
#> 
#> ====================================================================================================
#>  
#> Model No. 2
#>     Model name:                domainreading
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainreading
#>     Number of linking items:   68
#>     Number of testlets:        11
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainreading': 4 of 68 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T02_07'.
#>    Iteration 2: Exclude item 'T01_07'.
#>    Iteration 3: Exclude item 'T02_04'.
#>    Iteration 4: Exclude item 'T05_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.232     0.045
#> 2 iterativ    1       T02_07       -2.505            0.195     0.042
#> 3 iterativ    2       T01_07       -2.267            0.160     0.033
#> 4 iterativ    3       T02_04       -1.176            0.142     0.033
#> 5 iterativ    4       T05_04       -0.867            0.129     0.035
#> 
eq.2_3  <- equat1pl(results2[[2]], itemFromRes(results2[[3]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE)
#> Found 2 model(s).
#>    Equating is executed for each dimension in each model separately.
#> 
#> ====================================================================================================
#>  
#> Model No. 1
#>     Model name:                domainlistening
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainlistening
#>     Number of linking items:   96
#>     Number of testlets:        9
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainlistening': 6 of 96 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T19_11'.
#>    Iteration 2: Exclude item 'T14_04'.
#>    Iteration 3: Exclude item 'T20_30'.
#>    Iteration 4: Exclude item 'T15_18'.
#>    Iteration 5: Exclude item 'T20_03'.
#>    Iteration 6: Exclude item 'T12_05'.
#>    Iteration 7: Exclude item 'T20_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.052     0.055
#> 2 iterativ    1       T19_11        2.534            0.078     0.050
#> 3 iterativ    2       T14_04         2.53            0.105     0.047
#> 4 iterativ    3       T20_30        1.749            0.124     0.041
#> 5 iterativ    4       T15_18        1.435            0.139     0.047
#> 6 iterativ    5       T20_03        0.704            0.147     0.044
#> 7 iterativ    6       T12_05       -0.678            0.140     0.042
#> 8 iterativ    7       T20_04        0.654            0.147     0.041
#> 
#> 
#> ====================================================================================================
#>  
#> Model No. 2
#>     Model name:                domainreading
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainreading
#>     Number of linking items:   119
#>     Number of testlets:        15
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainreading': 7 of 119 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T02_07'.
#>    Iteration 2: Exclude item 'T05_01'.
#>    Iteration 3: Exclude item 'T23_14'.
#>    Iteration 4: Exclude item 'T04_08'.
#>    Iteration 5: Exclude item 'T03_03'.
#>    Iteration 6: Exclude item 'T10_10'.
#>    Iteration 7: Exclude item 'T02_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.098     0.035
#> 2 iterativ    1       T02_07        3.381            0.126     0.033
#> 3 iterativ    2       T05_01       -1.045            0.117     0.029
#> 4 iterativ    3       T23_14       -0.992            0.109     0.029
#> 5 iterativ    4       T04_08        0.915            0.117     0.027
#> 6 iterativ    5       T03_03       -0.687            0.111     0.026
#> 7 iterativ    6       T10_10         0.68            0.117     0.023
#> 8 iterativ    7       T02_04        0.654            0.123     0.024
#> 
eq.1_3  <- equat1pl(results2[[1]], itemFromRes(results2[[3]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE)
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
#>     Number of testlets:        4
#>     Linking method:            Mean.Mean
#> 
#> Dimension 'domainlistening': 4 of 36 items with linking |DIF| > 0.64 identified.
#>    Iteration 1: Exclude item 'T15_02'.
#>    Iteration 2: Exclude item 'T12_05'.
#>    Iteration 3: Exclude item 'T14_04'.
#> 
#>     method iter itemExcluded DIF.excluded linking.constant linkerror
#> 1 iterativ    0                                      0.355     0.079
#> 2 iterativ    1       T15_02        1.744            0.406     0.070
#> 3 iterativ    2       T12_05       -1.244            0.370     0.059
#> 4 iterativ    3       T14_04        0.885            0.396     0.053
#> 
#> 
#> ====================================================================================================
#>  
#> Model No. 2
#>     Model name:                domainreading
#>     Number of dimension(s):    1
#>     Name(s) of dimension(s):   domainreading
#>     Number of linking items:   68
#>     Number of testlets:        11
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
#> 1 iterativ    0                                      0.305     0.059
#> 2 iterativ    1       T01_07       -1.656            0.280     0.057
#> 3 iterativ    2       T02_07        0.849            0.293     0.059
#> 4 iterativ    3       T09_07        0.734            0.304     0.056
#> 5 iterativ    4       T04_04         0.67            0.314     0.052
#> 6 iterativ    5       T07_10        0.644            0.325     0.052
#> 
lErrors3<- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)
#> Dimension 'domainlistening': Negative covariance of -0.021 between DIF for measurement 1 vs. 2, and DIF for measurement 2 vs. 3. Covariance is set to 0 to avoid underestimation of the true linking error.
#> 
#> Dimension 'domainlistening': Direct linking error of the three combinations of measurements: 
#>         mzp1 mzp2 N.Items    SD   Var linkerror chained mzp_1.vs.3 mzp_1.vs.2.vs.3
#>       1    1    2      34 0.246 0.060     0.052      NA         NA              NA
#>       2    1    3      33 0.256 0.065     0.053   0.066      0.396           0.397
#>       3    2    3      89 0.253 0.064     0.041      NA         NA              NA
#> Dimension 'domainreading': Negative covariance of -0.013 between DIF for measurement 1 vs. 2, and DIF for measurement 2 vs. 3. Covariance is set to 0 to avoid underestimation of the true linking error.
#> 
#> Dimension 'domainreading': Direct linking error of the three combinations of measurements: 
#>         mzp1 mzp2 N.Items    SD   Var linkerror chained mzp_1.vs.3 mzp_1.vs.2.vs.3
#>       1    1    2      64 0.234 0.055     0.035      NA         NA              NA
#>       2    1    3      63 0.295 0.087     0.052   0.043      0.325           0.251
#>       3    2    3     112 0.253 0.064     0.024      NA         NA              NA


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
