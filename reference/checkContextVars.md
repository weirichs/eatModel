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
#> Error in checkContextVars(datW[, "language"], "language", type = "DIF",     itemdata = datW[, -c(1:4)]): argument "renam" is missing, with no default
chk1$info
#> Error: object 'chk1' not found
```
