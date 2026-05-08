# Deviance Plot for mirt objects

Plots the deviance change in every iteration for IRT models estimated
with [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.html)

## Usage

``` r
plotDevianceMirt(mirt.obj, omitUntil=1, adaptWindow = TRUE)
```

## Arguments

- mirt.obj:

  The object returned by
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.html)

- omitUntil:

  An optional value indicating number of iterations to be omitted for
  plotting.

- adaptWindow:

  Logical: Should the size of the output window be adjusted to
  accommodate the additional information about the R session and
  available memory? This option should not be used if the plots are not
  to be displayed directly but written to a PDF file.
