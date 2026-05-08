# Extracts the Q3 matrices from the object returned by getResults()

Collects Q3 matrices in a list of data frames.

## Usage

``` r
q3FromRes(resultsObj, out = c("wide", "long" ), triangular = FALSE)
```

## Arguments

- resultsObj:

  The object returned by
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).

- out:

  Specifies the output format. `"wide"` gives a triangular matrix,
  `"long"` gives a long format table with the residual correlation in a
  separate column.

- triangular:

  Logical: should the wide-format matrix be arranged in triangular
  shape?

## Value

A list of data frames, one data.frame per model.

## Examples

``` r
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
q3   <- q3FromRes(res)
```
