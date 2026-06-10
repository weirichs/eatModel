# Contains the transformation rules from the IQB studies

This functions allows the transformation rules (anchor parameters,
population reference values, cut scores) for the large-scale studies
conducted at the IQB to be retrieved from a centrally stored database.
The function arguments are used to specify in detail for which study,
subject, competency areas, and mode the values are to be requested.

## Usage

``` r
getTrafo(dataBase = system.file("extdata", "trafo.rda", package = "eatModel"),
         mode=c("paper","pc"), grade=c("primary", "secondary"),
         subject = c("math", "deu", "eng", "frz", "bio", "che", "phy"),
         domain = c("all", "GL", "ZO", "RF", "MS", "GM", "DHW", "ZA", "ME",
                    "FZ", "DZ", "lesen", "hoeren", "ortho", "sg",
                    "CE", "CF", "PE", "PF", "BE", "BF"),
         study = c("bt", "vera"))
```

## Arguments

- dataBase:

  Either the path to the R object containing the database, or the
  database object itself, which has already been loaded.

- mode:

  Specify the desired mode, e.g., paper-based assessment or
  computer-based assessment. Note that the data for computer-based
  studies will be added to the database gradually as they become
  available.

- grade:

  Specify for which grade level parameters are to be requested.

- subject:

  The desired subject. Only for the science subjects (i.e., bio, phy,
  che), it is possible to use more than obe subject. See the examples
  for more details.

- domain:

  The desired domain. One or several domains can be chosen. Use `"all"`,
  if you want to list all domains for the selected subject.

- study:

  Since some large-scale studies use different transformation rules than
  VERA, the specific study to which the transformation rules refer must
  be specified here. For some subjects (such as the science subjects
  biology, chemistry, and physics), however, there are no transformation
  rules for VERA. For other areas (such as the global domain in
  mathematics for the primary level), the transformation rules for
  Bildungstrend and VERA are identical.

## Details

The function returns anchor parameters, the mean, and the standard
deviation of the reference population, as well as cut scores for various
subjects, domains, grade levels, and modes, in the format required by
[`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md)
or
[`transformToBista`](https://weirichs.github.io/eatModel/reference/transformToBista.md).
The values returned by the function can be used directly as arguments
for the
[`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md)
and
[`transformToBista`](https://weirichs.github.io/eatModel/reference/transformToBista.md)
functions. For example, the returned list object `anchor` can be used as
an argument for the
[`equat1pl`](https://weirichs.github.io/eatModel/reference/equat1pl.md)
function. The returned list objects `refPop` and `cuts` can both be used
as arguments for the
[`transformToBista`](https://weirichs.github.io/eatModel/reference/transformToBista.md)
function.

The information stored in the database is based on this file:
I:/Methoden/10_sonstige Materialien/Transformationsvorschriften.xlsx

The following is a brief description of the domains stored in the
database.

**Mode: paper**:

- **Deutsch primary:** Originally, the scores from the 2007 norm study
  and re-standardization were used for this purpose. The
  re-standardization became necessary because the mean scores differed
  implausibly widely between the pilot study and the norm study. In the
  2011 Laendervergleich, it became apparent that the link to the norm
  parameters was functioning poorly, particularly for the listening
  domain, so the Laendervergleich was used to establish a new norm
  reference. Orthografy was included in the Laendervergleich only in a
  self-weighted subsample excluding students with special educational
  needs (SEN). Since this student group is also part of the population
  but was not included in the sampling until 2016, the population mean
  from 2011 was shifted by the difference that would have resulted had
  SEN students been included in the sampling. As a result, the cut
  scores also had to be shifted. This procedure, known at the IQB as the
  "backward trend," results in non-integer cut scores. Due to
  psychometric problems and severe violation of local stochastic
  independence assumption, the "Sprachgebrauch" domain was not included
  in the Laendervergleich and Bildungstrend (BT) studies. Consequently,
  only the norm values used for VERA are available for this domain.

- **Math primary:** There was no "backward trend" in mathematics. The
  cut scores for all domains are also identical here. The same applies
  to the mean and standard deviation of the reference population.

- **Deutsch secondary:** The reference population for VERA is the 2008
  norming sample. For the Bildungstrend (BT), the reference is the
  BT2015, not—as one might think—the Laendervergleich 2009. The reason,
  once again, is that no SEN students were included in the 2009
  laendervergleich. This also leads to a "backward trend" and shifted
  cut scores, but only for the Bildungstrend, not for VERA.

## Value

A list with four objects.

- anchor:

  The data.frame with anchor parameters

- refPop:

  A data.frame with mean and standard deviation for the reference
  population

- cuts:

  A list with cut scores

- info:

  Character string with short information about refrence population

## References

## Examples

``` r
# transformation rules for all domains subject deutsch, primary level
primDeu <- getTrafo(mode="paper", grade="primary", subject = "deu",
           domain = "all", study = "vera")
#> Warning: cannot open compressed file 'I:/Methoden/10_sonstige Materialien/trafo.rda', probable reason 'No such file or directory'
#> Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection

# transformation rules for all science domains
secScien<- getTrafo(mode="paper", grade="secondary", subject = c("bio", "che", "phy"),
           domain = "all", study = "bt")
#> Warning: cannot open compressed file 'I:/Methoden/10_sonstige Materialien/trafo.rda', probable reason 'No such file or directory'
#> Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
```
