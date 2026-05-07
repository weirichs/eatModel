# Gives correlation estimates from the object returned by getResults()

Gives correlation estimates for multidimensional models.

## Usage

``` r
correlationFromRes(resultsObj, digits = NULL)
```

## Arguments

- resultsObj:

  The object returned by
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).

- digits:

  Optional: integer value if coefficients should be rounded.

## Value

A list of data frame(s) with contain(s) correlation matrices.
