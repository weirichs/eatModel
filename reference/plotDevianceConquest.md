# Deviance Plot for Conquest log files

Plots the deviance change in every iteration.

## Usage

``` r
plotDevianceConquest(logFile, omitUntil=1, reverse=TRUE, change=TRUE, adaptWindow = TRUE)
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

- adaptWindow:

  Logical: Should the size of the output window be adjusted to
  accommodate the additional information about the R session and
  available memory? This option should not be used if the plots are not
  to be displayed directly but written to a PDF file.

## Author

Martin Hecht, Sebastian Weirich
