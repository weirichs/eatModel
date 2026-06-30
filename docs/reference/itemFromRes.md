# Collect all item results from the object returned by getResults()

Collect items results in a wide data frame.

## Usage

``` r
itemFromRes(resultsObj)
```

## Arguments

- resultsObj:

  The object returned by
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).

## Value

A data frame in the wide format with several columns, containing model
name, item name, dimension, number of valid responses per item,
percentage of correct responses, item discrimination, item difficulty
parameter, optional 2pl discrimination parameters, standard errors,
infit and outfit.

## Examples

``` r
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
items<- itemFromRes(res)
```
