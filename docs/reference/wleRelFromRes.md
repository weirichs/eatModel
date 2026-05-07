# Gives WLE reliability from the object returned by getResults()

WLE reliability in a data frame.

## Usage

``` r
wleRelFromRes(resultsObj)
```

## Arguments

- resultsObj:

  The object returned by
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).

## Value

A data frame with three columns, containing model name, dimension name,
and the corresponding WLE reliability.

## Examples

``` r
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
wleRel <- wleRelFromRes(res)
```
