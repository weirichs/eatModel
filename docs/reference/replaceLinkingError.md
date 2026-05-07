# Output interface for [`multiEquatError`](https://weirichs.github.io/eatModel/reference/multiEquatError.md)

[`multiEquatError`](https://weirichs.github.io/eatModel/reference/multiEquatError.md)
computes linking errors for three measurement occasions.
`replaceLinkingError` serves to integrate this linking error into the
object created by
[`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md).
The linking error then can be transformed with subsequent function
[`transformToBista`](https://weirichs.github.io/eatModel/reference/transformToBista.md).
For an illustration of the workflow, see the examples included in the
help page of
[`multiEquatError`](https://weirichs.github.io/eatModel/reference/multiEquatError.md).

## Usage

``` r
replaceLinkingError (equatingList, multiEquatError_output, verbose = TRUE, digits = 4)
```

## Arguments

- equatingList:

  The object returned by
  [`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md),
  for measurement occasion 1 vs. 3.

- multiEquatError_output:

  The object returned by
  [`multiEquatError`](https://weirichs.github.io/eatModel/reference/multiEquatError.md)
  which contains the linking error for measurement occasion 1 vs. 3.

- verbose:

  Logical: Print information about old (original) and new (replaced)
  linking error to console?

- digits:

  Number of decimals for printing information about old (original) and
  new (replaced) linking error

## Value

The 'equatingList' object with replaced linking errors. This object can
be processed subsequently by
[`transformToBista`](https://weirichs.github.io/eatModel/reference/transformToBista.md)

## Examples

``` r
# see the examples in the help page of 'multiEquatError'
```
