# Collects all plausible values from the object returned by getResults()

Collects plausible values in a wide or long format data frame.

## Usage

``` r
pvFromRes(resultsObj, toWideFormat = TRUE, idVarName = NULL, verbose=TRUE)
```

## Arguments

- resultsObj:

  The object returned by
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).

- toWideFormat:

  Logical: Should the plausible values be arranged in the wide format?
  Alternatively, they will appear in the long format.

- idVarName:

  Optional: character vector of individual student id. This is only to
  provide compatibility with older package versions. Specification of
  this argument is only necessary if the function gives an error.

- verbose:

  Print messages to console?

## Value

A data frame in the wide or long format with several columns, containing
model name, model dimension(s), student ID, and several columns with the
plausible values.

## Examples

``` r
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
pv   <- pvFromRes(res)
```
