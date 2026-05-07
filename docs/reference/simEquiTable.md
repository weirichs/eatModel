# Computes equivalence table for the Rasch and partial credit model

Function provides the equivalence table for unidimensional 1pl models,
specifying the individual competence level for each possible total score
of the test.

## Usage

``` r
simEquiTable  ( anchor, item = NULL, cat = NULL, value = NULL, mRef,
       sdRef, addConst = 500, multConst = 100, cutScores)
```

## Arguments

- anchor:

  A data frame with anchor parameters on the logit scale, transformed to
  the metric of the reference population. Data frame must have at least
  two columns: The first column contains the names of all anchored
  items. The second column contains anchor parameters. If the
  equivalence table is to be generated for a partial credit model, an
  additional column for the respective item category must be specified
  in the data.frame. The data.frame will then have three columns. In
  this case, the additional arguments `item`, `cat`, and `value` must be
  specified. The various item parameters for a partial credit item must
  therefore be listed one below the other in the data.frame. See the
  second example for further details.

- item:

  Optional: Give the number or name of the item identifier column in
  `anchor`. Only necessary for partial credit models.

- cat:

  Optional: Give the number or name of the category column in `anchor`.
  Only necessary for partial credit models.

- value:

  Optional: Give the number or name of the parameter column in `anchor`.
  Only necessary for partial credit models.

- mRef:

  Scalar: mean of the reference population.

- sdRef:

  Scalar: Standard deviation of the reference population.

- addConst:

  Additive constant for parameter transformation.

- multConst:

  Multiplicative constant for parameter transformation.

- cutScores:

  Named list of one or two elements. "values" is a numeric vector of cut
  scores (increasing), "labels" is an optional character vector of cut
  score labels. Note that "labels" (if specified) has to be of one more
  length than "values".

## Value

A list of two data frames, including the complete table and the reduced
table with the following 5 columns.

- score:

  Students raw score

- est:

  Estimated individual WLE according to the raw score.

- bista:

  Transformed WLE

- ks:

  competence level

## Author

Johannes Schult, Sebastian Weirich

## Examples

``` r
### Example 1: equivalence table for Rasch models
# read anchor parameter
file <- system.file("extdata", "results.rda", package = "eatModel")
load(file)

# use domain 'reading'
prm  <- subset(itemFromRes(res), model == "komplesen")

# use bista cut scores
cuts   <- list ( values = 390+0:3*75, labels = c("I", "II", "III", "IV", "V") )

# create the equivalence table
ret <- simEquiTable( anchor = prm[,c("item", "est")], cutScores = cuts , mRef = 0.039, sdRef = 1.071)

### Example 2: equivalence table for partial credit model
# load partial credit data
data(reading)
#> Warning: data set 'reading' not found

# To speed up the process, the tests will be administered only for Test
# Booklet 8. This booklet contains questions with some partial credit items.
d   <- subset(reading, bookletID == "TH08")
#> Error: object 'reading' not found
dw  <- reshape2::dcast(d, idstud~item, value.var = "valueSum")
#> Error: object 'd' not found
defT<- defineModel(dat = dw, items = -1, id = "idstud",  irtmodel = "PCM",software="tam")
#> Error: object 'dw' not found
runT<- runModel(defT)
#> Error: object 'defT' not found
resT<- getResults(runT, omitPV=TRUE, Q3 = FALSE)
#> Error: object 'runT' not found
it  <- itemFromRes(resT)
#> Error: object 'resT' not found

# create equivalence table with arbitrary cuts and reference values
ret2<- simEquiTable( anchor = it, item = "item", cat="category", value = "est",
       cutScores = list ( values = c(400, 600)), mRef = 0.047, sdRef = 1.181)
#> Error in simEquiTable(anchor = it, item = "item", cat = "category", value = "est",     cutScores = list(values = c(400, 600)), mRef = 0.047, sdRef = 1.181): unused arguments (item = "item", cat = "category", value = "est")
```
