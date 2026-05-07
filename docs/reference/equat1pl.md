# 1pl equating with optional elimination of linking DIF items

Function does the 1pl linking according to
[`equating.rasch`](https://rdrr.io/pkg/sirt/man/equating.rasch.html)
from the `sirt` package. Moreover, optional elimination of items with
linking DIF is allowed and linking error may be estimated via the
jackknife method if testlets are specified.

## Usage

``` r
equat1pl(results, prmNorm, item = NULL, domain = NULL, testlet = NULL, cat=NULL,
         value = NULL, excludeLinkingDif = TRUE, difBound = 1, iterativ = FALSE,
         method = c("Mean.Mean", "Haebara", "Stocking.Lord", "robust", "Haberman"), itemF = NULL,
         domainF = NULL, testletF = NULL, valueF = NULL, estimation=c("OLS", "BSQ", "HUB", "MED", "LTS", "L1", "L0"),
         b_trim=Inf, lts_prop=.5)
```

## Arguments

- results:

  The object returned by
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).
  Alternatively, a data.frame with item parameters of the focus group.
  In this case, additional arguments (`itemF`, `domainF`, and `valueF`)
  have to be defined.

- prmNorm:

  Data frame with normed anchor item parameters. Data frame must have at
  least two columns: items and item difficulties. Use the further
  arguments `item`, `domain`, `testlet` and `value` to define the column
  in which the corresponding parameter can be found. If `item`,
  `domain`, `testlet` and `value` is NULL, `prmNorm` must have only two
  columns: First column items, second column item difficulties. Column
  names then are irrelevant.

- item:

  Optional: Give the number or name of the item identifier column in
  prmNorm.

- domain:

  Optional: Give the number or name of the domain name in prmNorm. Only
  necessary if item identifiers are not unique in `prmNorm` (for
  example, if one item occurs several times, with a global item
  parameter and a domain-specific item parameter. Domain names in
  prmNorm must match dimension names in the object returned by
  [`getResults`](https://weirichs.github.io/eatModel/reference/getResults.md).

- testlet:

  Optional: Give the number or name of the testlet name in prmNorm. Only
  necessary if linking errors should be estimated via the jackknife
  method.

- cat:

  Optional: Give the number or name of the category column in `prmNorm`.
  Only necessary for partial credit models. See example 8 in the help
  file of
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)

- value:

  Optional: Give the number or name of the parameter column in
  `prmNorm`.

- excludeLinkingDif:

  Logical. Should items with linking DIF excluded?

- difBound:

  Defines the boundary. Items with absolute linking DIF greater than the
  boundary will be removed from the linking procedure.

- iterativ:

  Logical. Should the exclusion of linking DIF items executed in an
  iterative loop? (i.e. start with all items and compute linking
  constant, than remove the item with the most pronounced DIF and
  compute linking constant, and so on, until no item is left with
  \|DIF\|\> `difBound`.

- method:

  Linking method. If `"Mean.Mean"`, `"Haebara"`, or `"Stocking.Lord"`,
  the function
  [`equating.rasch`](https://rdrr.io/pkg/sirt/man/equating.rasch.html)
  from the `sirt` package is called. If `"robust"`, the function
  [`linking.robust`](https://rdrr.io/pkg/sirt/man/linking.robust.html)
  from the `sirt` package is called. If `"Haberman"`, the function
  [`linking.haberman`](https://rdrr.io/pkg/sirt/man/linking.haberman.html)
  from the `sirt` package is called.

- itemF:

  Optional: Give the number or name of the item column in results.

- domainF:

  Optional: Give the number or name of the domain column in results.

- testletF:

  Optional: Give the number or name of the testlet column in teh results
  object. Only necessary if linking errors should be estimated via the
  jackknife method.

- valueF:

  Optional: Give the number or name of the parameter column in results.

- estimation:

  Applies only if `method` equals `"Haberman"`. Estimation method. See
  the help page of `linking.haberman` from the `sirt` package for
  further details.

- b_trim:

  Applies only if `method` equals `"Haberman"`. Trimming parameter for
  item slopes. See the help page of `linking.haberman` from the `sirt`
  package for further details.

- lts_prop:

  Applies only if `method` equals `"Haberman"`. Proportion of retained
  observations in `"LTS"` regression estimation. See the help page of
  `linking.haberman` from the `sirt` package for further details.

## Value

A list with equating information intended for further transformation by
the
[`transformToBista`](https://weirichs.github.io/eatModel/reference/transformToBista.md)
function.

## Examples

``` r
# see example 5, 6, 6a, and 8 in the help file of defineModel()
```
