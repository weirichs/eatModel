# Read Conquest Plausible Values Output Files

This function reads Conquest plausible value files and automatically
identifies the number of cases, the number of plausible values and the
number of dimensions.

## Usage

``` r
get.plausible(file, quiet = FALSE, forConquestResults = FALSE)
```

## Arguments

- file:

  Character string with the name of the ConQuest plausible values file.

- quiet:

  Logical: Suppress printing messages on console?

- forConquestResults:

  This argument only applies if `get.plausible` is called by
  `getResults` to enhance performance.

## Value

A data frame with one row per person containing the following columns:

- case:

  Case number

- ID:

  Identifier for this case

- pv:

  Plausible values columns. Columns are named pv.\[number of
  PV\].\[number of dimension\]. For example, pv.5.2 refers to the 5th
  plausible value of the second dimension.

- eap:

  Expected a posteriori (EAP) ability estimate for each person. Columns
  are named eap_Dim.\[number of dimension\]. For example, eap_Dim.2
  refers to the eap estimate of the second dimension.

- se.eap:

  Standard error of the EAP estimate. Columns are named
  se.eap_Dim.\[number of dimension\]. For example, se.eap_Dim.2 refers
  to the standard error of the EAP estimate of the second dimension.

## Examples

``` r
file <- system.file("extdata", "twodim.pvl", package = "eatModel")
pv   <- get.plausible(file)
#> 2907 persons and 2 dimensions(s) found.
#> 5 plausible values were drawn for each person on each dimension.
```
