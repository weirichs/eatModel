# Read Conquest WLE or MLE Output Files

Read Conquest files comprising maximum likelihood estimates (MLE) or
weighted likelihood estimates (WLE).

## Usage

``` r
get.wle(file)
```

## Arguments

- file:

  Character string with the name of the ConQuest MLE or WLE file.

## Value

A data frame with one row per person containing the following columns:

- case:

  Case number

- ID:

  Identifier for this case

- n.solved:

  Number of items this person answered correctly

- n.total:

  Number of total items presented to this person

- wle:

  WLE or MLE estimate. The last number of the columns name indicates the
  dimension the WLE or MLE estimate belongs to.

- std.wle:

  Standard error of WLE or MLE estimate. The last number of the columns
  name indicates the dimension the WLE or MLE estimate belongs to.

## References

See pp. 230 of Wu, M.L., Adams, R.J., Wilson, M.R., & Haldane, S.A.
(2007). *ACER ConQuest Version 2.0. Generalised Item Response Modeling
Software.* Camberwell, Victoria: ACER Press.

## Examples

``` r
file <- system.file("extdata", "twodim.wle", package = "eatModel")
wle  <- get.wle(file)
#> Found valid WLEs of 2907 person(s) for 2 dimension(s).
```
