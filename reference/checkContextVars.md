# check for consistency of context variables and items

Function checks whether some types of context vars, i.e. group, DIF and
weighting variables, are consistent with item variables. The function is
mainly used for internal consistency checks.

## Usage

``` r
checkContextVars (x, varname, type = c("weight", "DIF", "group", "HG"), itemdata,
                  suppressAbort = FALSE, internal = FALSE, renam)
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

- renam:

  Optional: A data.frame with two columns containing the renamed item
  names. This is necessary because, when using Conquest, item names
  cannot exceed ten characters. If there are more than ten characters,
  the items are temporarily renamed internally, and after the Conquest
  analysis is complete, the renaming is reversed. To ensure that the
  console output contains the original item names, a renaming back to
  the original is also performed for each console output. The
  specification of the `renam` argument is only necessary when the
  function is called internally; otherwise, the value `NULL` should
  always be used here.

## Value

A list

## Examples

``` r
data(trends)
# first reshape the data for the first time of measurement set into wide format
datW <- reshape2::dcast(trends[which(trends[,"year"] == 2010),],
                        idstud+sex+ses+language~item, value.var="value")
chk1 <- checkContextVars(datW[,"language"], "language", type="DIF",
                         itemdata = datW[,-c(1:4)], renam=NULL)
#> Warning: Following 14 items are constants in DIF variable 'language', group other:
#>    T01_01, T05_04, T07_04, T07_07, T07_08, T07_10, T09_04, T09_05, T09_06, T10_08, T12_05, T13_06, T15_10, T16_04
#> Warning: For 58 items, some response categories in some DIF groups have less than 3 valid responses: 'T01_03', 'T01_04', 'T01_05', 'T01_07', 'T02_01', 'T02_02', 'T02_06', 'T02_07', 'T03_03', 'T04_02', 'T04_06', 'T04_07', 'T05_02', 'T06_03', 'T06_05', 'T07_01', 'T07_02', 'T07_03', 'T07_05', 'T07_06', 'T07_09', 'T08_01', 'T08_02', 'T08_03', 'T08_05', 'T08_06', 'T09_02', 'T09_03', 'T09_07', 'T09_08', 'T09_09', 'T09_10', 'T09_11', 'T10_02', 'T10_06', 'T10_07', 'T11_01', 'T11_03', 'T11_04', 'T11_06', 'T11_09', 'T11_10', 'T12_01', 'T12_03', 'T12_04', 'T12_06', 'T12_09', 'T13_08', 'T14_01', 'T14_03', 'T14_04', 'T14_06', 'T15_02', 'T16_02', 'T16_07', 'T16_09', 'T16_10', 'T16_11'. 
#>    Remove these items because otherwise the IRT DIF model probably will crash.
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
