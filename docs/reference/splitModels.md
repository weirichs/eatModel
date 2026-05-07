# Create several models by splitting the qMatrix and/or person.groups

[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
is programmed to define a single model. With `splitModels` several
models can be set up. The output of `splitModels` can be directly passed
to the `splittedModels` argument of
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)

## Usage

``` r
splitModels ( qMatrix = NULL , person.groups = NULL , 
    split = c ( "qMatrix" , "person.groups" ) , add = NULL , cross = NULL , 
    all.persons = TRUE , all.persons.lab = "all" , 
    person.split.depth = 0:length(person.groups[,-1,drop=FALSE]) , 
    full.model.names = TRUE , model.name.elements = c ( "dim" , "group" , "cross" ) , 
    include.var.name = FALSE , env = FALSE , nCores=NULL , mcPackage = c("future", "parallel"),
    GBcore=NULL , verbose = TRUE )
```

## Arguments

- qMatrix:

  Same argument as in
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md):

  Optional: A named data frame indicating how items should be grouped to
  dimensions. The first column contains the names of all items and
  should be named item. The other columns contain dimension definitions
  and should be named with the respective dimension names. A positive
  value (e.g., 1 or 2 or 1.4) indicates the loading weight with which an
  item loads on the dimension, a value of 0 indicates that the
  respective item does not load on this dimension. If no q matrix is
  specified by the user, an unidimensional structure is assumed.

- person.groups:

  data.frame, first row must be person ID, further columns contain group
  categories, e.g. data.frame ( "id" = 1:10 , "sex" = sample ( c (
  "male" , "female" ) , 10 , replace = TRUE ) )

- split:

  character, possible values and their consequences:

  `NULL`: qMatrix and person.groups are not split, one model with
  original qMatrix and all persons is set up `"qMatrix"`: qMatrix is
  split into single dimensions, number of created models equals number
  of dimensions `"person.groups"`: person.groups is split into single
  groups, number of created models equals number of all combinations of
  groups (with at least one person) `c("qMatrix","person.groups")`:
  default, both qMatrix and person.groups is split and single dimensions
  and single groups are crossed, number of created models equals number
  of dimension multiplied with number of all combinations of groups

- add:

  list of elements with single values, names of elements should be
  arguments of
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md),
  elements are the value that is passed when running
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md);
  elements in `add` are used for all models; e.g. list ( "software" =
  "conquest" , "nodes" = 15 ), that means that all models are estimated
  with "conquest" and 15 nodes

- cross:

  list of elements with several values, names of elements should be
  arguments of
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md),
  elements are the value that is passed when running
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md);
  elements in `cross` are crossed into models; e.g. list ( "software" =
  c("conquest","tam") , "nodes" = c(15,30) ), now all models are set up
  to run once with "conquest" and once with "tam", and with 15 and 30
  nodes

- all.persons:

  logical (default: TRUE), for each group variable in `person.groups` an
  "all" category is included

- all.persons.lab:

  character, name of the "all" category

- person.split.depth:

  integer, depth of group splits, 0: global all persons are included, 1:
  groups of all variables are included, 2: groups of all pairs of
  variables are included, n: groups of n variables are included. Can be
  a vector with more than one argument, e.g. for 3 variables, the full
  number of splits (which is also the default) can be obtained by
  c(0,1,2,3); this creates a model with all persons (0), all groups of
  all variables (1), groups from pairs of variables (2), and groups from
  combining all 3 variables (3). Using `person.split.depth` usually
  makes most sense if `all.persons=TRUE`; if `all.persons=FALSE` the
  depth equals the number of variables (if another depth is set, no
  splits will be performed).

- full.model.names:

  logical (default: TRUE), model names are derived from
  `model.name.elements`; if FALSE models are numbered in ascending order

- model.name.elements:

  character, elements that model names are built of, possible values:
  "dim" , "group" , "add" , "cross"; default: c ( "dim" , "group" ,
  "cross" ) , that means that model names include the name of the
  dimension(s) , group(s) , and parameter values that are crossed in

- include.var.name:

  logical (default: FALSE), include the name of the variable when
  building model names; e.g. (FALSE) "science\_\_sex.female\_\_conquest"
  , (TRUE) "dim.science\_\_group.sex.female\_\_software.conquest"

- env:

  logical (default: FALSE) (FALSE) returns a list with two elements:
  data.frame with model information (model overview), list with model
  specifications (TRUE) returns a list with two elements: data.frame
  with model information (model overview), list of environments with
  model specifications set as objects for intended subsequent use with
  [`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
  use `env=FALSE`

- nCores:

  integer (default: NULL), number of cores to use for subsequent data
  preparation, model estimation and results compilation

- mcPackage:

  Which package should be used for local host definition in multicore
  processing? If R version \< 3.4, `"parallel"` is recommended. If R
  version \>= 3.4, `"future"` is recommended.

- GBcore:

  numeric (default: NULL), maximum RAM usage per core in giga bytes

- verbose:

  logical (default: TRUE), print progress

## Value

depending on `env` either: (env=FALSE) returns a list with two elements:
data.frame with model information (model overview), list with model
specifications (env=TRUE) returns a list with two elements: data.frame
with model information (model overview), list of environments with model
specifications set as objects for intended subsequent use with
[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)
use `env=FALSE`

## Author

Martin Hecht

## See also

[`defineModel`](https://weirichs.github.io/eatModel/reference/defineModel.md)

## Examples

``` r
# see also examples in 'defineModel'
# example qMatrix
qMatrix <- data.frame ( "item" = 1:4 , "science" = c(1,1,0,0) ,
           "math" = c(0,0,1,1) , stringsAsFactors = FALSE )

# example person.groups
person.groups <- data.frame ( "person" = 1:4 , "state" = rep(c("Berlin","Bavaria"),2) ,
                 "sex" = c(rep("female",2),rep("male",2)) , stringsAsFactors = FALSE )

# Example 1: one 2-dimensional model with all persons (no split)
m01 <- splitModels ( qMatrix=qMatrix, person.groups=person.groups, split=NULL )
#> --------------------------------
#> splitModels: generating 1 models
#> .
#> see <returned>$models
#> number of cores: 1
#> --------------------------------
m01$models
#>   model.no                      model.name                    model.subpath
#> 1        1 science_math__state.all_sex.all ./science_math/state.all_sex.all
#>            dim Ndim             group Ngroup
#> 1 science_math    2 state.all_sex.all      1

# Example 2: split qMatrix to create two unidimensional models, each with all persons
m02 <- splitModels ( qMatrix=qMatrix, person.groups=person.groups, split=c("qMatrix") )
#> --------------------------------
#> splitModels: generating 2 models
#> ..
#> see <returned>$models
#> number of cores: 2
#> --------------------------------
m02$models
#>   model.no                 model.name               model.subpath     dim Ndim
#> 1        1 science__state.all_sex.all ./science/state.all_sex.all science    1
#> 2        2    math__state.all_sex.all    ./math/state.all_sex.all    math    1
#>               group Ngroup
#> 1 state.all_sex.all      1
#> 2 state.all_sex.all      1

# Example 3: split person.groups to create 2-dimensional models, each with a
# subgroup of persons
m03 <- splitModels ( qMatrix=qMatrix, person.groups=person.groups,
       split=c("person.groups") )
#> --------------------------------
#> splitModels: generating 9 models
#> .........
#> see <returned>$models
#> number of cores: 9
#> --------------------------------
m03$models
#>   model.no                             model.name
#> 1        1  science_math__state.Berlin_sex.female
#> 2        2    science_math__state.Berlin_sex.male
#> 3        3     science_math__state.Berlin_sex.all
#> 4        4 science_math__state.Bavaria_sex.female
#> 5        5   science_math__state.Bavaria_sex.male
#> 6        6    science_math__state.Bavaria_sex.all
#> 7        7     science_math__state.all_sex.female
#> 8        8       science_math__state.all_sex.male
#> 9        9        science_math__state.all_sex.all
#>                             model.subpath          dim Ndim
#> 1  ./science_math/state.Berlin_sex.female science_math    2
#> 2    ./science_math/state.Berlin_sex.male science_math    2
#> 3     ./science_math/state.Berlin_sex.all science_math    2
#> 4 ./science_math/state.Bavaria_sex.female science_math    2
#> 5   ./science_math/state.Bavaria_sex.male science_math    2
#> 6    ./science_math/state.Bavaria_sex.all science_math    2
#> 7     ./science_math/state.all_sex.female science_math    2
#> 8       ./science_math/state.all_sex.male science_math    2
#> 9        ./science_math/state.all_sex.all science_math    2
#>                      group Ngroup
#> 1  state.Berlin_sex.female      1
#> 2    state.Berlin_sex.male      1
#> 3     state.Berlin_sex.all      1
#> 4 state.Bavaria_sex.female      1
#> 5   state.Bavaria_sex.male      1
#> 6    state.Bavaria_sex.all      1
#> 7     state.all_sex.female      1
#> 8       state.all_sex.male      1
#> 9        state.all_sex.all      1

# Example 4: split both qMatrix and person.groups to create unidimensional
# models for all subgroups
m04 <- splitModels ( qMatrix=qMatrix, person.groups=person.groups,
       split=c("qMatrix","person.groups") )
#> ---------------------------------
#> splitModels: generating 18 models
#> ..................
#> see <returned>$models
#> number of cores: 18
#> ---------------------------------
m04$models
#>    model.no                        model.name
#> 1         1  science__state.Berlin_sex.female
#> 2         2    science__state.Berlin_sex.male
#> 3         3     science__state.Berlin_sex.all
#> 4         4 science__state.Bavaria_sex.female
#> 5         5   science__state.Bavaria_sex.male
#> 6         6    science__state.Bavaria_sex.all
#> 7         7     science__state.all_sex.female
#> 8         8       science__state.all_sex.male
#> 9         9        science__state.all_sex.all
#> 10       10     math__state.Berlin_sex.female
#> 11       11       math__state.Berlin_sex.male
#> 12       12        math__state.Berlin_sex.all
#> 13       13    math__state.Bavaria_sex.female
#> 14       14      math__state.Bavaria_sex.male
#> 15       15       math__state.Bavaria_sex.all
#> 16       16        math__state.all_sex.female
#> 17       17          math__state.all_sex.male
#> 18       18           math__state.all_sex.all
#>                         model.subpath     dim Ndim                    group
#> 1   ./science/state.Berlin_sex.female science    1  state.Berlin_sex.female
#> 2     ./science/state.Berlin_sex.male science    1    state.Berlin_sex.male
#> 3      ./science/state.Berlin_sex.all science    1     state.Berlin_sex.all
#> 4  ./science/state.Bavaria_sex.female science    1 state.Bavaria_sex.female
#> 5    ./science/state.Bavaria_sex.male science    1   state.Bavaria_sex.male
#> 6     ./science/state.Bavaria_sex.all science    1    state.Bavaria_sex.all
#> 7      ./science/state.all_sex.female science    1     state.all_sex.female
#> 8        ./science/state.all_sex.male science    1       state.all_sex.male
#> 9         ./science/state.all_sex.all science    1        state.all_sex.all
#> 10     ./math/state.Berlin_sex.female    math    1  state.Berlin_sex.female
#> 11       ./math/state.Berlin_sex.male    math    1    state.Berlin_sex.male
#> 12        ./math/state.Berlin_sex.all    math    1     state.Berlin_sex.all
#> 13    ./math/state.Bavaria_sex.female    math    1 state.Bavaria_sex.female
#> 14      ./math/state.Bavaria_sex.male    math    1   state.Bavaria_sex.male
#> 15       ./math/state.Bavaria_sex.all    math    1    state.Bavaria_sex.all
#> 16        ./math/state.all_sex.female    math    1     state.all_sex.female
#> 17          ./math/state.all_sex.male    math    1       state.all_sex.male
#> 18           ./math/state.all_sex.all    math    1        state.all_sex.all
#>    Ngroup
#> 1       1
#> 2       1
#> 3       1
#> 4       1
#> 5       1
#> 6       1
#> 7       1
#> 8       1
#> 9       1
#> 10      1
#> 11      1
#> 12      1
#> 13      1
#> 14      1
#> 15      1
#> 16      1
#> 17      1
#> 18      1

# Example 5: set "software"="conquest" and "method"="montecarlo" for all models
m05 <- splitModels ( qMatrix=qMatrix, person.groups=person.groups ,
       add = list ( "software"="conquest" , "method"="montecarlo" ) )
#> ---------------------------------
#> splitModels: generating 18 models
#> ..................
#> see <returned>$models
#> number of cores: 18
#> ---------------------------------
m05$models
#>    model.no                        model.name
#> 1         1  science__state.Berlin_sex.female
#> 2         2    science__state.Berlin_sex.male
#> 3         3     science__state.Berlin_sex.all
#> 4         4 science__state.Bavaria_sex.female
#> 5         5   science__state.Bavaria_sex.male
#> 6         6    science__state.Bavaria_sex.all
#> 7         7     science__state.all_sex.female
#> 8         8       science__state.all_sex.male
#> 9         9        science__state.all_sex.all
#> 10       10     math__state.Berlin_sex.female
#> 11       11       math__state.Berlin_sex.male
#> 12       12        math__state.Berlin_sex.all
#> 13       13    math__state.Bavaria_sex.female
#> 14       14      math__state.Bavaria_sex.male
#> 15       15       math__state.Bavaria_sex.all
#> 16       16        math__state.all_sex.female
#> 17       17          math__state.all_sex.male
#> 18       18           math__state.all_sex.all
#>                         model.subpath     dim Ndim                    group
#> 1   ./science/state.Berlin_sex.female science    1  state.Berlin_sex.female
#> 2     ./science/state.Berlin_sex.male science    1    state.Berlin_sex.male
#> 3      ./science/state.Berlin_sex.all science    1     state.Berlin_sex.all
#> 4  ./science/state.Bavaria_sex.female science    1 state.Bavaria_sex.female
#> 5    ./science/state.Bavaria_sex.male science    1   state.Bavaria_sex.male
#> 6     ./science/state.Bavaria_sex.all science    1    state.Bavaria_sex.all
#> 7      ./science/state.all_sex.female science    1     state.all_sex.female
#> 8        ./science/state.all_sex.male science    1       state.all_sex.male
#> 9         ./science/state.all_sex.all science    1        state.all_sex.all
#> 10     ./math/state.Berlin_sex.female    math    1  state.Berlin_sex.female
#> 11       ./math/state.Berlin_sex.male    math    1    state.Berlin_sex.male
#> 12        ./math/state.Berlin_sex.all    math    1     state.Berlin_sex.all
#> 13    ./math/state.Bavaria_sex.female    math    1 state.Bavaria_sex.female
#> 14      ./math/state.Bavaria_sex.male    math    1   state.Bavaria_sex.male
#> 15       ./math/state.Bavaria_sex.all    math    1    state.Bavaria_sex.all
#> 16        ./math/state.all_sex.female    math    1     state.all_sex.female
#> 17          ./math/state.all_sex.male    math    1       state.all_sex.male
#> 18           ./math/state.all_sex.all    math    1        state.all_sex.all
#>    Ngroup software     method
#> 1       1 conquest montecarlo
#> 2       1 conquest montecarlo
#> 3       1 conquest montecarlo
#> 4       1 conquest montecarlo
#> 5       1 conquest montecarlo
#> 6       1 conquest montecarlo
#> 7       1 conquest montecarlo
#> 8       1 conquest montecarlo
#> 9       1 conquest montecarlo
#> 10      1 conquest montecarlo
#> 11      1 conquest montecarlo
#> 12      1 conquest montecarlo
#> 13      1 conquest montecarlo
#> 14      1 conquest montecarlo
#> 15      1 conquest montecarlo
#> 16      1 conquest montecarlo
#> 17      1 conquest montecarlo
#> 18      1 conquest montecarlo

# Example 6: cross "nodes"=c(1000,5000) and "seed"=c(1234,4321) into all models
m06 <- splitModels ( qMatrix=qMatrix, person.groups=person.groups ,
       add = list ( "software"="conquest" , "method"="montecarlo" ) ,
       cross = list ( "nodes"=c(1000,5000) , "seed"=c(1234,4321) ) )
#> ---------------------------------
#> splitModels: generating 72 models
#> ........................................................................
#> see <returned>$models
#> number of cores: 32
#> ---------------------------------
m06$models
#>    model.no                                    model.name
#> 1         1  science__state.Berlin_sex.female__1000__1234
#> 2         2  science__state.Berlin_sex.female__1000__4321
#> 3         3  science__state.Berlin_sex.female__5000__1234
#> 4         4  science__state.Berlin_sex.female__5000__4321
#> 5         5    science__state.Berlin_sex.male__1000__1234
#> 6         6    science__state.Berlin_sex.male__1000__4321
#> 7         7    science__state.Berlin_sex.male__5000__1234
#> 8         8    science__state.Berlin_sex.male__5000__4321
#> 9         9     science__state.Berlin_sex.all__1000__1234
#> 10       10     science__state.Berlin_sex.all__1000__4321
#> 11       11     science__state.Berlin_sex.all__5000__1234
#> 12       12     science__state.Berlin_sex.all__5000__4321
#> 13       13 science__state.Bavaria_sex.female__1000__1234
#> 14       14 science__state.Bavaria_sex.female__1000__4321
#> 15       15 science__state.Bavaria_sex.female__5000__1234
#> 16       16 science__state.Bavaria_sex.female__5000__4321
#> 17       17   science__state.Bavaria_sex.male__1000__1234
#> 18       18   science__state.Bavaria_sex.male__1000__4321
#> 19       19   science__state.Bavaria_sex.male__5000__1234
#> 20       20   science__state.Bavaria_sex.male__5000__4321
#> 21       21    science__state.Bavaria_sex.all__1000__1234
#> 22       22    science__state.Bavaria_sex.all__1000__4321
#> 23       23    science__state.Bavaria_sex.all__5000__1234
#> 24       24    science__state.Bavaria_sex.all__5000__4321
#> 25       25     science__state.all_sex.female__1000__1234
#> 26       26     science__state.all_sex.female__1000__4321
#> 27       27     science__state.all_sex.female__5000__1234
#> 28       28     science__state.all_sex.female__5000__4321
#> 29       29       science__state.all_sex.male__1000__1234
#> 30       30       science__state.all_sex.male__1000__4321
#> 31       31       science__state.all_sex.male__5000__1234
#> 32       32       science__state.all_sex.male__5000__4321
#> 33       33        science__state.all_sex.all__1000__1234
#> 34       34        science__state.all_sex.all__1000__4321
#> 35       35        science__state.all_sex.all__5000__1234
#> 36       36        science__state.all_sex.all__5000__4321
#> 37       37     math__state.Berlin_sex.female__1000__1234
#> 38       38     math__state.Berlin_sex.female__1000__4321
#> 39       39     math__state.Berlin_sex.female__5000__1234
#> 40       40     math__state.Berlin_sex.female__5000__4321
#> 41       41       math__state.Berlin_sex.male__1000__1234
#> 42       42       math__state.Berlin_sex.male__1000__4321
#> 43       43       math__state.Berlin_sex.male__5000__1234
#> 44       44       math__state.Berlin_sex.male__5000__4321
#> 45       45        math__state.Berlin_sex.all__1000__1234
#> 46       46        math__state.Berlin_sex.all__1000__4321
#> 47       47        math__state.Berlin_sex.all__5000__1234
#> 48       48        math__state.Berlin_sex.all__5000__4321
#> 49       49    math__state.Bavaria_sex.female__1000__1234
#> 50       50    math__state.Bavaria_sex.female__1000__4321
#> 51       51    math__state.Bavaria_sex.female__5000__1234
#> 52       52    math__state.Bavaria_sex.female__5000__4321
#> 53       53      math__state.Bavaria_sex.male__1000__1234
#> 54       54      math__state.Bavaria_sex.male__1000__4321
#> 55       55      math__state.Bavaria_sex.male__5000__1234
#> 56       56      math__state.Bavaria_sex.male__5000__4321
#> 57       57       math__state.Bavaria_sex.all__1000__1234
#> 58       58       math__state.Bavaria_sex.all__1000__4321
#> 59       59       math__state.Bavaria_sex.all__5000__1234
#> 60       60       math__state.Bavaria_sex.all__5000__4321
#> 61       61        math__state.all_sex.female__1000__1234
#> 62       62        math__state.all_sex.female__1000__4321
#> 63       63        math__state.all_sex.female__5000__1234
#> 64       64        math__state.all_sex.female__5000__4321
#> 65       65          math__state.all_sex.male__1000__1234
#> 66       66          math__state.all_sex.male__1000__4321
#> 67       67          math__state.all_sex.male__5000__1234
#> 68       68          math__state.all_sex.male__5000__4321
#> 69       69           math__state.all_sex.all__1000__1234
#> 70       70           math__state.all_sex.all__1000__4321
#> 71       71           math__state.all_sex.all__5000__1234
#> 72       72           math__state.all_sex.all__5000__4321
#>                                   model.subpath     dim Ndim
#> 1   ./science/state.Berlin_sex.female/1000/1234 science    1
#> 2   ./science/state.Berlin_sex.female/1000/4321 science    1
#> 3   ./science/state.Berlin_sex.female/5000/1234 science    1
#> 4   ./science/state.Berlin_sex.female/5000/4321 science    1
#> 5     ./science/state.Berlin_sex.male/1000/1234 science    1
#> 6     ./science/state.Berlin_sex.male/1000/4321 science    1
#> 7     ./science/state.Berlin_sex.male/5000/1234 science    1
#> 8     ./science/state.Berlin_sex.male/5000/4321 science    1
#> 9      ./science/state.Berlin_sex.all/1000/1234 science    1
#> 10     ./science/state.Berlin_sex.all/1000/4321 science    1
#> 11     ./science/state.Berlin_sex.all/5000/1234 science    1
#> 12     ./science/state.Berlin_sex.all/5000/4321 science    1
#> 13 ./science/state.Bavaria_sex.female/1000/1234 science    1
#> 14 ./science/state.Bavaria_sex.female/1000/4321 science    1
#> 15 ./science/state.Bavaria_sex.female/5000/1234 science    1
#> 16 ./science/state.Bavaria_sex.female/5000/4321 science    1
#> 17   ./science/state.Bavaria_sex.male/1000/1234 science    1
#> 18   ./science/state.Bavaria_sex.male/1000/4321 science    1
#> 19   ./science/state.Bavaria_sex.male/5000/1234 science    1
#> 20   ./science/state.Bavaria_sex.male/5000/4321 science    1
#> 21    ./science/state.Bavaria_sex.all/1000/1234 science    1
#> 22    ./science/state.Bavaria_sex.all/1000/4321 science    1
#> 23    ./science/state.Bavaria_sex.all/5000/1234 science    1
#> 24    ./science/state.Bavaria_sex.all/5000/4321 science    1
#> 25     ./science/state.all_sex.female/1000/1234 science    1
#> 26     ./science/state.all_sex.female/1000/4321 science    1
#> 27     ./science/state.all_sex.female/5000/1234 science    1
#> 28     ./science/state.all_sex.female/5000/4321 science    1
#> 29       ./science/state.all_sex.male/1000/1234 science    1
#> 30       ./science/state.all_sex.male/1000/4321 science    1
#> 31       ./science/state.all_sex.male/5000/1234 science    1
#> 32       ./science/state.all_sex.male/5000/4321 science    1
#> 33        ./science/state.all_sex.all/1000/1234 science    1
#> 34        ./science/state.all_sex.all/1000/4321 science    1
#> 35        ./science/state.all_sex.all/5000/1234 science    1
#> 36        ./science/state.all_sex.all/5000/4321 science    1
#> 37     ./math/state.Berlin_sex.female/1000/1234    math    1
#> 38     ./math/state.Berlin_sex.female/1000/4321    math    1
#> 39     ./math/state.Berlin_sex.female/5000/1234    math    1
#> 40     ./math/state.Berlin_sex.female/5000/4321    math    1
#> 41       ./math/state.Berlin_sex.male/1000/1234    math    1
#> 42       ./math/state.Berlin_sex.male/1000/4321    math    1
#> 43       ./math/state.Berlin_sex.male/5000/1234    math    1
#> 44       ./math/state.Berlin_sex.male/5000/4321    math    1
#> 45        ./math/state.Berlin_sex.all/1000/1234    math    1
#> 46        ./math/state.Berlin_sex.all/1000/4321    math    1
#> 47        ./math/state.Berlin_sex.all/5000/1234    math    1
#> 48        ./math/state.Berlin_sex.all/5000/4321    math    1
#> 49    ./math/state.Bavaria_sex.female/1000/1234    math    1
#> 50    ./math/state.Bavaria_sex.female/1000/4321    math    1
#> 51    ./math/state.Bavaria_sex.female/5000/1234    math    1
#> 52    ./math/state.Bavaria_sex.female/5000/4321    math    1
#> 53      ./math/state.Bavaria_sex.male/1000/1234    math    1
#> 54      ./math/state.Bavaria_sex.male/1000/4321    math    1
#> 55      ./math/state.Bavaria_sex.male/5000/1234    math    1
#> 56      ./math/state.Bavaria_sex.male/5000/4321    math    1
#> 57       ./math/state.Bavaria_sex.all/1000/1234    math    1
#> 58       ./math/state.Bavaria_sex.all/1000/4321    math    1
#> 59       ./math/state.Bavaria_sex.all/5000/1234    math    1
#> 60       ./math/state.Bavaria_sex.all/5000/4321    math    1
#> 61        ./math/state.all_sex.female/1000/1234    math    1
#> 62        ./math/state.all_sex.female/1000/4321    math    1
#> 63        ./math/state.all_sex.female/5000/1234    math    1
#> 64        ./math/state.all_sex.female/5000/4321    math    1
#> 65          ./math/state.all_sex.male/1000/1234    math    1
#> 66          ./math/state.all_sex.male/1000/4321    math    1
#> 67          ./math/state.all_sex.male/5000/1234    math    1
#> 68          ./math/state.all_sex.male/5000/4321    math    1
#> 69           ./math/state.all_sex.all/1000/1234    math    1
#> 70           ./math/state.all_sex.all/1000/4321    math    1
#> 71           ./math/state.all_sex.all/5000/1234    math    1
#> 72           ./math/state.all_sex.all/5000/4321    math    1
#>                       group Ngroup software     method nodes seed
#> 1   state.Berlin_sex.female      1 conquest montecarlo  1000 1234
#> 2   state.Berlin_sex.female      1 conquest montecarlo  1000 4321
#> 3   state.Berlin_sex.female      1 conquest montecarlo  5000 1234
#> 4   state.Berlin_sex.female      1 conquest montecarlo  5000 4321
#> 5     state.Berlin_sex.male      1 conquest montecarlo  1000 1234
#> 6     state.Berlin_sex.male      1 conquest montecarlo  1000 4321
#> 7     state.Berlin_sex.male      1 conquest montecarlo  5000 1234
#> 8     state.Berlin_sex.male      1 conquest montecarlo  5000 4321
#> 9      state.Berlin_sex.all      1 conquest montecarlo  1000 1234
#> 10     state.Berlin_sex.all      1 conquest montecarlo  1000 4321
#> 11     state.Berlin_sex.all      1 conquest montecarlo  5000 1234
#> 12     state.Berlin_sex.all      1 conquest montecarlo  5000 4321
#> 13 state.Bavaria_sex.female      1 conquest montecarlo  1000 1234
#> 14 state.Bavaria_sex.female      1 conquest montecarlo  1000 4321
#> 15 state.Bavaria_sex.female      1 conquest montecarlo  5000 1234
#> 16 state.Bavaria_sex.female      1 conquest montecarlo  5000 4321
#> 17   state.Bavaria_sex.male      1 conquest montecarlo  1000 1234
#> 18   state.Bavaria_sex.male      1 conquest montecarlo  1000 4321
#> 19   state.Bavaria_sex.male      1 conquest montecarlo  5000 1234
#> 20   state.Bavaria_sex.male      1 conquest montecarlo  5000 4321
#> 21    state.Bavaria_sex.all      1 conquest montecarlo  1000 1234
#> 22    state.Bavaria_sex.all      1 conquest montecarlo  1000 4321
#> 23    state.Bavaria_sex.all      1 conquest montecarlo  5000 1234
#> 24    state.Bavaria_sex.all      1 conquest montecarlo  5000 4321
#> 25     state.all_sex.female      1 conquest montecarlo  1000 1234
#> 26     state.all_sex.female      1 conquest montecarlo  1000 4321
#> 27     state.all_sex.female      1 conquest montecarlo  5000 1234
#> 28     state.all_sex.female      1 conquest montecarlo  5000 4321
#> 29       state.all_sex.male      1 conquest montecarlo  1000 1234
#> 30       state.all_sex.male      1 conquest montecarlo  1000 4321
#> 31       state.all_sex.male      1 conquest montecarlo  5000 1234
#> 32       state.all_sex.male      1 conquest montecarlo  5000 4321
#> 33        state.all_sex.all      1 conquest montecarlo  1000 1234
#> 34        state.all_sex.all      1 conquest montecarlo  1000 4321
#> 35        state.all_sex.all      1 conquest montecarlo  5000 1234
#> 36        state.all_sex.all      1 conquest montecarlo  5000 4321
#> 37  state.Berlin_sex.female      1 conquest montecarlo  1000 1234
#> 38  state.Berlin_sex.female      1 conquest montecarlo  1000 4321
#> 39  state.Berlin_sex.female      1 conquest montecarlo  5000 1234
#> 40  state.Berlin_sex.female      1 conquest montecarlo  5000 4321
#> 41    state.Berlin_sex.male      1 conquest montecarlo  1000 1234
#> 42    state.Berlin_sex.male      1 conquest montecarlo  1000 4321
#> 43    state.Berlin_sex.male      1 conquest montecarlo  5000 1234
#> 44    state.Berlin_sex.male      1 conquest montecarlo  5000 4321
#> 45     state.Berlin_sex.all      1 conquest montecarlo  1000 1234
#> 46     state.Berlin_sex.all      1 conquest montecarlo  1000 4321
#> 47     state.Berlin_sex.all      1 conquest montecarlo  5000 1234
#> 48     state.Berlin_sex.all      1 conquest montecarlo  5000 4321
#> 49 state.Bavaria_sex.female      1 conquest montecarlo  1000 1234
#> 50 state.Bavaria_sex.female      1 conquest montecarlo  1000 4321
#> 51 state.Bavaria_sex.female      1 conquest montecarlo  5000 1234
#> 52 state.Bavaria_sex.female      1 conquest montecarlo  5000 4321
#> 53   state.Bavaria_sex.male      1 conquest montecarlo  1000 1234
#> 54   state.Bavaria_sex.male      1 conquest montecarlo  1000 4321
#> 55   state.Bavaria_sex.male      1 conquest montecarlo  5000 1234
#> 56   state.Bavaria_sex.male      1 conquest montecarlo  5000 4321
#> 57    state.Bavaria_sex.all      1 conquest montecarlo  1000 1234
#> 58    state.Bavaria_sex.all      1 conquest montecarlo  1000 4321
#> 59    state.Bavaria_sex.all      1 conquest montecarlo  5000 1234
#> 60    state.Bavaria_sex.all      1 conquest montecarlo  5000 4321
#> 61     state.all_sex.female      1 conquest montecarlo  1000 1234
#> 62     state.all_sex.female      1 conquest montecarlo  1000 4321
#> 63     state.all_sex.female      1 conquest montecarlo  5000 1234
#> 64     state.all_sex.female      1 conquest montecarlo  5000 4321
#> 65       state.all_sex.male      1 conquest montecarlo  1000 1234
#> 66       state.all_sex.male      1 conquest montecarlo  1000 4321
#> 67       state.all_sex.male      1 conquest montecarlo  5000 1234
#> 68       state.all_sex.male      1 conquest montecarlo  5000 4321
#> 69        state.all_sex.all      1 conquest montecarlo  1000 1234
#> 70        state.all_sex.all      1 conquest montecarlo  1000 4321
#> 71        state.all_sex.all      1 conquest montecarlo  5000 1234
#> 72        state.all_sex.all      1 conquest montecarlo  5000 4321

# Example 7: list elements in cross that contain more than one element need to be
# lists themselves
m07 <- splitModels ( qMatrix=qMatrix, person.groups=NULL ,
       cross = list ( "regression"=list( c("sex") , c("sex","state") ) ,
       "seed"=c(1234,4321) ) )
#> --------------------------------
#> splitModels: generating 8 models
#> ........
#> see <returned>$models
#> number of cores: 8
#> --------------------------------
m07$models
#>   model.no               model.name            model.subpath     dim Ndim group
#> 1        1       science__sex__1234       ./science/sex/1234 science    1  <NA>
#> 2        2       science__sex__4321       ./science/sex/4321 science    1  <NA>
#> 3        3 science__sex.state__1234 ./science/sex.state/1234 science    1  <NA>
#> 4        4 science__sex.state__4321 ./science/sex.state/4321 science    1  <NA>
#> 5        5          math__sex__1234          ./math/sex/1234    math    1  <NA>
#> 6        6          math__sex__4321          ./math/sex/4321    math    1  <NA>
#> 7        7    math__sex.state__1234    ./math/sex.state/1234    math    1  <NA>
#> 8        8    math__sex.state__4321    ./math/sex.state/4321    math    1  <NA>
#>   Ngroup regression seed
#> 1     NA        sex 1234
#> 2     NA        sex 4321
#> 3     NA  sex.state 1234
#> 4     NA  sex.state 4321
#> 5     NA        sex 1234
#> 6     NA        sex 4321
#> 7     NA  sex.state 1234
#> 8     NA  sex.state 4321

# Example 8: create an "empty" model without qMatrix and person.groups
m08 <- splitModels ( qMatrix=NULL, person.groups=NULL )
#> --------------------------------
#> splitModels: generating 1 models
#> .
#> see <returned>$models
#> number of cores: 1
#> --------------------------------
m08$models
#>   model.no model.name model.subpath  dim Ndim group Ngroup
#> 1        1     model1             . <NA>   NA  <NA>     NA
```
