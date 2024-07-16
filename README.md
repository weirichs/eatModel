# eatModel <a href="https://weirichs.github.io/eatModel/"><img src="man/figures/eatModel.svg" align="right" height="139" alt="placeholder logo" /></a>

<!-- badges: start -->
[![R-CMD-check](https://github.com/weirichs/eatModel/workflows/R-CMD-check/badge.svg)](https://github.com/weirichs/eatModel/actions)
<!-- badges: end -->

## Overview 

`eatModel` (Educational Assessment Tools for IRT Modeling) provides functions to compute IRT models using the software `Conquest` or the `R` package `TAM`.

## Installation

```R
# Install eatRep from GitHub via
remotes::install_github("weirichs/eatModel")
```

## View package documentation

```R
library(eatModel)
### View package documentation
package?eatModel
```

## Exemplary analysis

The package comes with an examplary data set containing fictional achievement scores of 13524 students in two domains (reading, listening) from three measurement occasions (2010, 2015, 2020) in the long format. The data are not longitudinal as the individuals stem from different cohorts and do not overlap across measurement occasions. The name of the data set is `trends`. You can find detailed exemplary analyses in the help pages of `defineModel()`. The following example estimates a between-item two-dimensional 2pl model with conditioning (background) model for the student cohort measured in 2020. The `splitModels()` function can be used previously if the corresponding model should be estimated in 2010 and 2015 likewise. 

```R
library(eatModel)
### load exemplary data 
data(trends)

### choose data for 2020 and prepare data in the wide format
dat2020<- reshape2::dcast(subset ( trends, year == 2020),
          idstud+country+sex+ses+language~item, value.var="value")

### create Q matrix indicating which item belongs to which latent dimension. 
qMat   <- unique(trends[ ,c("item","domain")])
qMat   <- data.frame ( qMat[,"item", drop=FALSE], model.matrix(~domain-1, data = qMat))

### specify the model
def    <- defineModel(dat = dat2020, id = "idstud", items = colnames(dat2020)[-c(1:5)],
          qMatrix = qMat, HG.var = c("sex", "ses", "language"), irtmodel = "2PL", software = "tam")

### run the model
run    <- runModel(def)

### collect the results
results<- getResults(run)

### extract item parameters from results
items  <- itemFromRes(results)

### extract WLEs from results
wles   <- wleFromREs(results)
```

