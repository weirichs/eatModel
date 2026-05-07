# Run IRT model specified by 'defineModel' using Conquest, TAM, ort mirt

Start the estimation of an IRT model defined by
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md).
`runModel` should be called with the argument returned by
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md).

## Usage

``` r
runModel(defineModelObj, show.output.on.console = FALSE, 
    show.dos.console = TRUE, wait = TRUE, onlySkeleton = FALSE)
```

## Arguments

- defineModelObj:

  The object returned by
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md).

- show.output.on.console:

  Applies only if
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
  previously was called with `software = "conquest"`. Logical: Should
  the output of the conquest console be printed on the R console during
  estimation?

- show.dos.console:

  Applies only if
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
  previously was called with `software = "conquest"`. Logical: Should
  the output of the conquest console be printed on screen?

- wait:

  Applies only if
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
  previously was called with `software = "conquest"`. A logical (not NA)
  indicating whether the R interpreter should wait for the command to
  finish, or run it asynchronously.

- onlySkeleton:

  Applies only if
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
  previously was called with `software = "mirt"`. If `TRUE`, the IRT
  model is not estimated. Instead, only the model specification is
  returned as a data.frame, allowing the user to check whether the model
  has been specified as desired.

## Value

If
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
previously was called with `software = "tam"`, the returned value is
nearly identically to the corresponding TAM output. Accordingly, if
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
previously was called with `software = "mirt"`, the returned value is
nearly identically to the corresponding mirt output. If
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
previously was called with `software = "conquest"`, the returned value
contains only internally used information useful for
[`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).

## Author

Sebastian Weirich

## Examples

``` r
# see examples in the help file of defineModel()
```
