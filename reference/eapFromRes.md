# Collects all EAP results from the object returned by getResults()

Collects EAP results in a wide data frame.

## Usage

``` r
eapFromRes(resultsObj, idVarName = NULL, verbose = TRUE)
```

## Arguments

- resultsObj:

  The object returned by
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).

- idVarName:

  Optional: character vector of individual student id. This is only to
  provide compatibility with older package versions. Specification of
  this argument is only necessary if the function gives an error.

- verbose:

  Print messages to console?

## Value

A data frame with several columns, model name, student ID, dimension
name, group name, EAP and the corresponding standard error.

## Examples

``` r
### read exemplary results object
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)
eap  <- eapFromRes(res)
```
