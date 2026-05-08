# Read ConQuest showfiles

Function reads Conquest showfiles and transforms them into an R list of
data frames.

## Usage

``` r
get.shw(file, dif.term, split.dif = TRUE, abs.dif.bound = 0.6,
    sig.dif.bound = 0.3, p.value = 0.9)
```

## Arguments

- file:

  Character string with the name of the ConQuest show file.

- dif.term:

  Optional: Character string. Name of the term considered to be
  DIF-term. Must match corresponding term in showfile.

- split.dif:

  Logical: When TRUE, DIF-Parameter are only given for Reference group.

- abs.dif.bound:

  When DIF-Parameter are evaluated, this specifies the critical value
  for absolute DIF.

- sig.dif.bound:

  When DIF-Parameter are evaluated, this specifies the critical value
  for confidence interval DIF.

- p.value:

  When DIF-Parameter are evaluated, this specifies the critical p-value
  for confidence interval DIF.

## Details

Function searches for ‘TERM’-statements in Conquest showfile and reads
the tables associated with. If one statement, for example
`"item*gender"` is specified to contain DIF analyses, the Conquest
output concerning this term is read in and some *additional* analyses
are conducted. To compute the absolute DIF, item parameter estimates of
the `"item*gender"` table are doubled. Confidence intervals for 90, 95
and 99 percent are computed via the standard error of item parameter
estimates in the `"item*gender"` table. If both criteria - absolute DIF
exceeds `abs.dif.bound` and the confidence interval does not include
`sig.dif.bound`, item is considered to have DIF.

## Value

A list of data frames, named by the ‘TERM’-statements in Conquest
showfile, plus additional data frame named `"regression"` with
regression coefficients when latent linear regression model was
specified in Conquest analysis, plus an additional data frame named
`"cov.structure"` with covariance and correlation matrix of latent
dimensions (if uni-dimensional model is specified, the variance of the
latent dimension is given instead), plus an additional data frame named
`"final.deviance"` with information about the final deviance and the
number of estimated parameters intended for model comparison. If one
term was specified as DIF-statement, the corresponding data frame is
augmented with additional columns for confidence intervals and
indicators specifying significant DIF. Each data frame corresponding to
a ‘TERM’ statement contains following columns:

- No.:

  Item number

- item:

  Name of item

- ESTIMATE:

  Estimated item difficulty

- ERROR:

  Standard error of the estimated item difficulty

- MNSQ:

  Item's outfit

- CI:

  The lower bound of the confidence interval if the item has an ideal
  outfit of 1. This is to check whether the empirical outfit lies within
  or without the confidence interval of an item with hypothetically
  ideal fit.

- CI.1:

  The upper bound of the corresponding confidence interval

- T:

  The *T* value concerning the outfit

- MNSQ.1:

  Item's infit

- CI.2:

  The lower bound of the confidence interval if the item has an ideal
  infit of 1. This is to check whether the empirical infit lies within
  or without the confidence interval of an item with hypothetically
  ideal fit.

- CI.3:

  The upper bound of the corresponding confidence interval

- T.1:

  The *T* value concerning the infit

- filename:

  Name of the file which was read in.

When latent regression was specified, the last element of the returned
list is a data frame with regression coefficients, corresponding to the
number of dimensions and the number of regressors. Regressor names,
regression coefficients and its standard errors are given for each
dimension. Rows represent the regressors, columns represent the latent
dimension to which the regression is fitted.

## Examples

``` r
file <- system.file("extdata", "twodim.shw", package = "eatModel")
shw  <- get.shw(file)
#> Found TERM 1: 'item' 
#> Found 2 dimension(s): Dimension 1, Dimension 2
#> Found 4 regressor(s).
#> Warning: NAs introduced by coercion
#> Warning: NAs introduced by coercion
```
