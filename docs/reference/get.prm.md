# Read Conquest prm (parameter) Output Files

Read Conquest files comprising item parameters (prm).

## Usage

``` r
get.prm(file)
```

## Arguments

- file:

  Character string with the name of the ConQuest prm file (This file is
  requested using the `"export par"` statement in the conquest syntax
  file).

## Value

A data frame with three columns:

- Case:

  Case number

- item:

  Identifier for this item

- parameter:

  parameter estimate for this item

## Examples

``` r
file <- system.file("extdata", "twodim.prm", package = "eatModel")
prm  <- get.prm(file)
```
