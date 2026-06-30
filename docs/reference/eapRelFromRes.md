# Gives EAP reliability from the object returned by getResults()

EAP reliability in a data frame.

## Usage

``` r
eapRelFromRes(resultsObj)
```

## Arguments

- resultsObj:

  The object returned by
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).

## Value

A data frame with three columns, model name, domain name, and the EAP
reliability.

## Examples

``` r
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
eapRel <- eapRelFromRes(res)
```
