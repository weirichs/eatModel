# Gives regression coefficients (betas) from the object returned by getResults()

Gives coefficients from latent regression model (conditioning model) in
a data frame.

## Usage

``` r
regcoefFromRes(resultsObj, digits = NULL)
```

## Arguments

- resultsObj:

  The object returned by
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).

- digits:

  Optional: integer value if coefficients should be rounded.

## Value

A list of data frame(s) with three columns, containing model name,
dimension name, and the corresponding regression coefficients.

## Examples

``` r
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
regs <- regcoefFromRes(res)
#> No regression coefficients found in results object.
```
