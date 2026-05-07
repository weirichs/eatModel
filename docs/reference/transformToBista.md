# Transformation of item and person parameters to the Bista metric.

Function uses output of
[`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md)
to provide two data.frames, one for the item parameters on the bista
metric and one for the person parameter (PVs) on the bista metric.

## Usage

``` r
transformToBista ( equatingList, refPop, cuts, weights = NULL, defaultM = 500, 
                   defaultSD = 100, q3bound= .2, roman = FALSE, vera = TRUE,
                   idVarName = NULL, years=NULL)
```

## Arguments

- equatingList:

  The object returned by
  [`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md).

- refPop:

  Optional: Data frame with at least three columns. First column
  indicates the domain name. Note that this name must match the domain
  names in the output of
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).
  Second column contains the mean of the reference population. Third
  column contains the standard deviation of the reference population.
  Fourth column optionally contains the transformed mean on the Bista
  metric of the reference population. Fifth column optionally contains
  the transformed standard deviation on the Bista metric of the
  reference population. If the fourth and fifth columns are missing,
  values will be defaulted to 500/100. If `refPop` is not specified,
  mean and SD will be computed from the data, optionally using weights
  (if the `weights` argument is specified).

- cuts:

  A named list with cut scores. Names of the list must match the domain
  names in the output of
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).
  Each element of the list is a named list with one or two elements—the
  cut scores (in ascending order) and (optionally) the labels of the
  stages. The first element (cut scores) must be named `values`. The
  second element (labels of the stages) must be named `labels`. See the
  examples of
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
  for further details.

- weights:

  Optional: a data.frame with two columns, first column is person
  identifier, second columns is individual caseweight. Necessary for the
  transformation of linking error for (ordered) factors and/or if
  descriptives of the reference population should be computed directly
  from the data. See the examples of
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
  for further details.

- defaultM:

  Mean of the reference population on the “bista” metric.

- defaultSD:

  Standard deviation of the reference population in the “bista” metric.

- q3bound:

  Define the absolute boundary of Q3 values which should be captured in
  the item parameter list according to the guidelines of the
  “Vergleichsarbeiten”.

- roman:

  Logical: Use roman numbers for competence level column in the
  shortened item parameter table dedicated for the “Vergleichsarbeiten”?

- vera:

  Logical: Prepare item parameter list according to the guidelines of
  the “Vergleichsarbeiten”?

- idVarName:

  Optional: character vector of individual student id. This is only to
  provide compatibility with older package versions. Specification of
  this argument is only necessary if the function gives an error.

- years:

  Optional: numeric vector with two elements, indicating the both years
  of assessment. Only necessary if additional linking error object
  should be created.

## Value

A list with three data frames: the first one contains original and
transformed item parameters and competence levels. The second one
contains original and transformed person parameters and competence
levels. The third one contains transformation information.

## Author

Sebastian Weirich

## Examples

``` r
# see example 5, 6, and 6a in the help file of defineModel()
```
