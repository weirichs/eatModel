# Compares two [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md) objects

Function facilitates quality checks by comparing two objects returned by
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
whether there are functionally equivalent.

## Usage

``` r
compareDefineModelObjects (ref, tar, round = TRUE, digits = 3)
```

## Arguments

- ref:

  The reference object returned by
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md).

- tar:

  The target object returned by
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md).

- round:

  Logical. Round decimal values before the comparison?

- digits:

  Only necessary if `round = TRUE`. Number of digits for rounding

## Value

No value is returned. Function prints messages to the console if
inequalities occur.

## Examples

``` r
data(trends)

# first reshape the data for the first time of measurement set into wide format
datW <- reshape2::dcast(trends[which(trends[,"year"] == 2010),],
                        idstud+sex+ses+language~item, value.var="value")
                        
# second, create the q matrix from the long format data frame
qMat <- unique(trends[ ,c("item","domain")])
qMat <- data.frame ( qMat[,"item", drop=FALSE], model.matrix(~domain-1, data = qMat))

# model 1
def1 <- defineModel(datW, items = -c(1:4), id="idstud", qMatrix = qMat, software="tam",
        HG.var = c("sex", "ses"), method="gauss")
#> Following 152 item(s) missed in data frame will be removed from Q matrix: 
#>     T12_11, T09_12, T01_08, T15_17, T15_18, T13_18, T10_10, T04_09, T04_08, T11_0X, T08_09, T08_08, T03_09, T06_1X, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> Background variable(s) 'sex' of class 
#>     'factor' will be converted to indicator variables.
#>     Variable 'sex' was converted to 1 indicator(s) with name(s) 'sexmale'.
#> 5 subject(s) solved each item: P00106 (17 correct), P00939 (17 correct), P00393 (20 correct) ... 
#> Dataset is completely linked.
#> 'gauss' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.
#> Q matrix specifies 2 dimension(s).

# model 2 (misspecified)
def2 <- defineModel(datW, items = -c(1:6), id="idstud", qMatrix = qMat[-c(20,35,40),],
        software="tam", HG.var = c("sex", "ses", "language"), method="quadrature", nodes = 30)
#> Following 154 item(s) missed in data frame will be removed from Q matrix: 
#>     T01_01, T01_02, T12_11, T09_12, T01_08, T15_17, T15_18, T13_18, T10_10, T04_09, T04_08, T11_0X, T08_09, T08_08, T03_09, T06_1X, T19_06, T19_08, T19_05, T19_11, T19_01, T19_02, T19_30, T19_03, T19_07, T19_04, T19_09, T26_03, T26_05, T26_07, T26_09, T26_06, T26_10, T26_08, T26_02, T26_04, T21_11, T18_01, T21_01, T18_05, T21_12, T18_13, T21_10, T18_11, T21_05, T21_13, T18_08, T18_07, T21_09, T18_04, T18_02, T21_02, T21_06, T18_09, T21_04, T18_12, T21_14, T21_03, T18_10, T27_01, T27_05, T27_07, T27_06, T27_04, T27_02, T27_03, T30_06, T29_04, T29_02, T30_02, T30_03, T30_01, T29_05, T29_03, T24_05, T24_08, T24_02, T24_10, T24_06, T24_24, T24_25, T24_09, T24_21, T24_23, T28_08, T28_09, T28_02, T28_03, T28_04, T28_01, T28_10, T28_06, T23_10, T23_13, T23_09, T23_15, T23_03, T23_04, T23_14, T23_07, T23_08, T23_02, T23_06, T23_05, T25_05, T25_22, T22_03, T25_06, T22_12, T25_09, T25_10, T25_08, T22_11, T22_02, T25_11, T25_07, T25_14, T25_12, T25_04, T22_08, T25_23, T22_10, T22_06, T25_13, T22_01, T20_05, T20_09, T20_30, T20_08, T20_12, T20_06, T20_07, T20_02, T20_13, T20_03, T20_01, T20_04, T20_11, T20_10, T32_04, T31_03, T31_01, T31_04, T32_03, T32_02, T32_01, T17_02, T17_09, T17_10, T17_04, T17_03, T17_08, T17_06, T17_01
#> Following 3 item(s) missed in Q matrix will be removed from data: 
#>     T02_06, T09_05, T09_10
#> Background variable(s) 'sex', 'language' of class 
#>     'factor', 'factor' will be converted to indicator variables.
#>     Variable 'sex' was converted to 1 indicator(s) with name(s) 'sexmale'.
#>     Variable 'language' was converted to 2 indicator(s) with name(s) 'languagenativeAndOther', 'languageother'.
#> 3 subject(s) do not solve any item:
#>    P04281 (11 false), P04290 (11 false), P00782 (12 false) ... 
#> 5 subject(s) solved each item: P00106 (15 correct), P00939 (15 correct), P00393 (18 correct) ... 
#> Dataset is completely linked.
#> Q matrix specifies 2 dimension(s).

# compare both outputs
compareDefineModelObjects (def1, def2)
#> Additional common variables (beyond the 'by'-variables) found: 'value'. Add suffixes '.x', '.y' to these variables in the result data.frame.
#> 2 of 130 unit(s) of merging variable 'item' from data set 'target' not included in data set 'reference'.
#> 5 of 133 unit(s) of merging variable 'item' from data set 'reference' not included in data set 'target'.
#> 8952 of 165773 unit(s) of merging variable combination 'ID'+'item' from data set 'target' not included in data set 'reference'.
#> 4356 of 161177 unit(s) of merging variable combination 'ID'+'item' from data set 'reference' not included in data set 'target'.
#> Mismatch for entry 'variablen' of 'all.Names' object: 'Lengths (131, 126) differ (string compare on first 126)'. '126 string mismatches'
#> Mismatch for entry 'HG.var' of 'all.Names' object: 'Lengths (2, 4) differ (string compare on first 2)'. '2 string mismatches'
#> Additional common variables (beyond the 'by'-variables) found: 'value'. Add suffixes '.x', '.y' to these variables in the result data.frame.
#> 5 of 131 unit(s) of merging variable 'item' from data set 'reference Q matrix' not included in data set 'target Q matrix'.
#> 10 of 262 unit(s) of merging variable combination 'item'+'dim' from data set 'reference Q matrix' not included in data set 'target Q matrix'.
#> Q matrix mismatch: ''is.NA' value mismatch: 10 in current 0 in target
#> Control mismatch for object 'nodes': 'Numeric: lengths (20, 30) differ'
```
