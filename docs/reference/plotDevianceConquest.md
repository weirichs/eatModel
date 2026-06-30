# Deviance Plot for Conquest log files

Plots the deviance change in every iteration.

## Usage

``` r
plotDevianceConquest(logFile, omitUntil=1, reverse=TRUE, change=TRUE)
```

## Arguments

- logFile:

  Character string of the Conquest log file

- omitUntil:

  An optional value indicating number of iterations to be omitted for
  plotting.

- reverse:

  A logical indicating whether the deviance change should be multiplied
  by minus 1. The default is `TRUE`.

- change:

  An optional logical indicating whether deviance change or the deviance
  should be plotted.

## Author

Martin Hecht, Sebastian Weirich
