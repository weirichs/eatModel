# Reading and listening achievement test data obtained from large-scale assessment at three times of measurement

This data set contains fictional achievement scores of 13524 students in
two domains (reading, listening) in the long format.

## Usage

``` r
data(trends)
```

## Format

'data.frame': 404281 obs. of 17 variables:

- year:

  Year of evaluation

- idstud:

  Student identifier

- idclass:

  Class identifier

- wgt:

  individual student weight

- jkzone:

  jackknife zone (primary sampling unit)

- jkrep:

  jackknife replicate

- country:

  The country an examinee stems from

- language:

  spoken language at home

- ses:

  student's socio economical status

- sex:

  student's sex

- domain:

  The domain the variable belongs to

- booklet:

  booklet identifier. equal booklet identifiers indicate equal booklets
  across years (assessment cycles)

- block:

  block identifier

- task:

  task identifier

- item:

  item identifier

- format:

  item format

- pos:

  position of the block within the booklet

- value:

  The response of the student to the item (0=incorrect; 1=correct)

## Source

Simulated data

## Examples

``` r
data(trends)
# number of students per year, country and domain
by(data=trends, INDICES = trends[,"year"], FUN = function(x) { tapply(x[,"idstud"], x[,c("country", "domain")], FUN = function(y){length(unique(y))})})
#> trends[, "year"]: 2010
#>           domain
#> country    listening reading
#>   countryA      1598    1598
#>   countryB      1309    1309
#>   countryC      1569    1569
#> ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
#> trends[, "year"]: 2015
#>           domain
#> country    listening reading
#>   countryA      1482    1502
#>   countryB      1220    1237
#>   countryC      1709    1777
#> ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
#> trends[, "year"]: 2020
#>           domain
#> country    listening reading
#>   countryA      1283    1363
#>   countryB      1207    1245
#>   countryC      1830    1865
# number of items per year, country and domain
by(data=trends, INDICES = trends[,"year"], FUN = function(x) { tapply(x[,"item"], x[,c("country", "domain")], FUN = function(y){length(unique(y))})})
#> trends[, "year"]: 2010
#>           domain
#> country    listening reading
#>   countryA        51      80
#>   countryB        51      80
#>   countryC        51      80
#> ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
#> trends[, "year"]: 2015
#>           domain
#> country    listening reading
#>   countryA        96     119
#>   countryB        96     119
#>   countryC        96     119
#> ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
#> trends[, "year"]: 2020
#>           domain
#> country    listening reading
#>   countryA       119     137
#>   countryB       119     137
#>   countryC       119     137

# no overlapping student IDs between assessment cycles
ids <- by(data=trends, INDICES = trends[,"year"], FUN = function (x) {unique(x[,"idstud"])})
length(intersect(ids[[1]], ids[[2]]))
#> [1] 0
length(intersect(ids[[1]], ids[[3]]))
#> [1] 0
length(intersect(ids[[2]], ids[[3]]))
#> [1] 0

# sampling weights substantially differ between countries due to stratified sampling
eatTools::roundDF(do.call("rbind",  by(data=trends, INDICES = trends[,c("year", "country")], FUN = function (x) {data.frame ( x[1,c("year", "country")], eatTools::descr(x[!duplicated(x[,"idstud"]),"wgt"])[,c("Minimum", "Maximum", "Mean", "Median", "SD")], stringsAsFactors = FALSE)})), digits = 3)
#>       year  country Minimum Maximum    Mean Median     SD
#> 577   2010 countryA   1.152   7.489   3.106  3.457  1.077
#> 12660 2015 countryA   1.062  12.738   3.280  3.185  1.526
#> 16615 2020 countryA   1.089  15.680   3.804  3.578  1.670
#> 1889  2010 countryB   8.312  18.480  12.039 11.375  2.716
#> 13428 2015 countryB   2.895  24.221  11.993 12.212  4.134
#> 17111 2020 countryB   5.542  83.014  12.157 11.694  5.402
#> 1     2010 countryC  76.281 337.244 102.684 99.498 24.658
#> 11716 2015 countryC   2.000 397.493  82.224 88.332 47.131
#> 16105 2020 countryC  48.145 249.383  74.794 68.996 22.566

# which booklets occur in which assessment cycles?
# see, for example: Bo01 only occurs 2010; Bo02 occurs 2010, 2015, and 2022; Bo83 occurs 2015 and 2020
reshape2::dcast(do.call("rbind", by(data=trends, INDICES = trends[,"year"], FUN = function (x) {data.frame ( x[1,"year", drop=FALSE], table(x[!duplicated(x[,"idstud"]),"booklet"]), stringsAsFactors = FALSE)})), year~Var1, value.var = "Freq")
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#>   year Bo01 Bo02 Bo03 Bo04 Bo05 Bo06 Bo07 Bo09 Bo13 Bo17 Bo18 Bo19 Bo20 Bo21 Bo22 Bo23 Bo24 Bo25 Bo26 Bo27 Bo28 Bo29 Bo30 Bo31 Bo32 Bo33 Bo34 Bo35 Bo36 Bo37 Bo38 Bo39 Bo40 Bo41 Bo42 Bo43 Bo44 Bo45 Bo46 Bo91 Bo08 Bo10 Bo11 Bo12 Bo14 Bo15 Bo16 Bo47 Bo48 Bo49 Bo50 Bo51 Bo52 Bo53 Bo54 Bo55 Bo56 Bo57 Bo58 Bo59 Bo60 Bo61 Bo62 Bo63 Bo64 Bo83 Bo84 Bo85 Bo86 Bo87 Bo88 Bo89 Bo90 Bo65 Bo66 Bo67 Bo68 Bo69 Bo70 Bo71 Bo72 Bo73 Bo74 Bo75 Bo76 Bo77 Bo78 Bo79 Bo80 Bo81 Bo82
#> 1 2010  183  172  188  168  185  167  184  168  177  162  166  154  169  193  179  180  171  190  177  176   44   43   45   41   47   47   66   58   55   56   55   52   48   46   49   44   42   38   48   43   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> 2 2015   NA  133   NA  136   NA   NA  119  145  126  136   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   94  114  127  153  116  100  123  111  119   79  124  109  106   89  115  104  147  112  142  134  116   77   91   83  119   80  129  106  108  105  153  105  131   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA
#> 3 2020   NA  105   NA  128   NA   NA  128  122  121  126   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA  112  115   97   94  120  137  122   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA  105  130  136  122  129  123   95  129  112   80   93  112  135   87  118  102   98  127  134  106  105  142  130  144  111  100

# which reading tasks occur in which assessment cycles?
# see, for example: T01 occurs 2010, 2015, and 2022; T27 only occurs 2020
reshape2::dcast(do.call("rbind", by(data=subset(trends,domain=="reading"), INDICES = subset(trends,domain=="reading")[,"year"], FUN = function (x) {data.frame ( x[1,"year", drop=FALSE], table(x[!duplicated(x[,"idstud"]),"task"]), stringsAsFactors = FALSE)})), year~Var1, value.var = "Freq")
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#> Warning: row names were found from a short variable and have been discarded
#>   year T01 T02 T03 T04 T05 T06 T07 T08 T09 T10 T11 T22 T23 T24 T25 T26 T27
#> 1 2010 478 600 522 271 237 197 456 444 532 388 351  NA  NA  NA  NA  NA  NA
#> 2 2015 201 321 362 173 316 139 481 277 344 346 190 361 335 376 294  NA  NA
#> 3 2020 206 410 286 132 255 181 391 203 237 325 136 289 303 347 408 183 181

# students nested in classes?
reformulas::isNested(trends[,"idstud"], trends[,"idclass"])
#> [1] TRUE
# items nested in tasks?
reformulas::isNested(trends[,"item"], trends[,"task"])
#> [1] TRUE
# tasks nested in blocks? no, few tasks occur in more than one block
reformulas::isNested(trends[,"task"], trends[,"block"])
#> [1] FALSE
# tasks nested in blocks for specific years?
by(data=trends, INDICES = trends[,"year"], FUN = function (y) {lme4::isNested(y[,"task"], y[,"block"]) })
#> trends[, "year"]: 2010
#> [1] TRUE
#> ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
#> trends[, "year"]: 2015
#> [1] FALSE
#> ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
#> trends[, "year"]: 2020
#> [1] FALSE
# blocks nested in booklets?
lme4::isNested(trends[,"block"], trends[,"booklet"])
#> [1] FALSE
```
