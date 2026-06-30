# check for consistency of context variables and items

Function checks whether some types of context vars, i.e. group, DIF and
weighting variables, are consistent with item variables. The function is
mainly used for internal consistency checks.

## Usage

``` r
checkContextVars (x, varname, type = c("weight", "DIF", "group", "HG"), itemdata,
                  suppressAbort = FALSE, internal = FALSE)
```

## Arguments

- x:

  A vector with values of the context variable (e.g., DIF variable)

- varname:

  Optional: name of the context variable

- type:

  Type of the context variable with following entries allowed: DIF,
  group, HG, or weight.

- itemdata:

  data.frame with item responses

- suppressAbort:

  Logical: should the function suppress abort if inconsistencies occur?

- internal:

  Logical: is only used for internal use. Recommend to set to FALSE.

## Value

A list

## Examples

``` r
data(trends)
# first reshape the data for the first time of measurement set into wide format
datW <- reshape2::dcast(trends[which(trends[,"year"] == 2010),],
                        idstud+sex+ses+language~item, value.var="value")
chk1 <- checkContextVars(datW[,"language"], "language", type="DIF",
                         itemdata = datW[,-c(1:4)])
#> Warning: Following 14 items are constants in DIF variable 'language', group other:
#>    T01_01, T05_04, T07_04, T07_07, T07_08, T07_10, T09_04, T09_05, T09_06, T10_08, T12_05, T13_06, T15_10, T16_04
#> Warning: row names were found from a short variable and have been discarded
chk1$info
#>     varname varlevel nCases     type   vars value nValue
#> 1  language    other     41 constant T01_01     1      9
#> 2  language    other     41 constant T05_04     1      8
#> 3  language    other     41 constant T07_04     0      5
#> 4  language    other     41 constant T07_07     0      5
#> 5  language    other     41 constant T07_08     1      5
#> 6  language    other     41 constant T07_10     0      5
#> 7  language    other     41 constant T09_04     0      5
#> 8  language    other     41 constant T09_05     1      5
#> 9  language    other     41 constant T09_06     1      5
#> 10 language    other     41 constant T10_08     0      9
#> 11 language    other     41 constant T12_05     1     14
#> 12 language    other     41 constant T13_06     0     13
#> 13 language    other     41 constant T15_10     0     17
#> 14 language    other     41 constant T16_04     1      9
```
