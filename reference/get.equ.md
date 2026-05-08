# Reads equivalence table created in Conquest analysis.

Reads Conquest files comprising equivalence tables for MLE or WLE
parameters.

## Usage

``` r
get.equ(file)
```

## Arguments

- file:

  Character string with the name of the Conquest equ file (This file is
  requested using the `"equivalence"` statement in the conquest syntax
  file).

## Value

A list of \\n+1\\ elements, with \\n\\ the number of dimensions in the
analysis. Each element is a data frame, whose name corresponds to the
name of the dimension the values belongs to. All data frames except the
last one give the transformation of each possible raw score to the WLE
or MLE score including its standard error. First column in each data
frame contains the raw score, second column the transformed WLE or MLE
score, third columns its standard error. The last element of the list
give some sparse information about the model specifications.

## References

See pp. 162 of Wu, M.L., Adams, R.J., Wilson, M.R., & Haldane, S.A.
(2007). *ACER ConQuest Version 2.0. Generalised Item Response Modeling
Software.* Camberwell, Victoria: ACER Press.

## Examples

``` r
file <- system.file("extdata", "twodim.equ", package = "eatModel")
equ  <- get.equ(file)
#> Found 2 dimension(s).
```
