# eatModel <a href="https://weirichs.github.io/eatModel/"><img src="man/figures/logo.png" align="right" height="120" alt="placeholder logo" /></a>

<!-- badges: start -->
[![R-CMD-check](https://github.com/weirichs/eatModel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/weirichs/eatModel/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview 

`eatModel` (Educational Assessment Tools for IRT Modeling) provides functions to compute IRT models (Rasch, 2PL, partial credit, generalized partial credit) using the software `Conquest` or the `R` packages `TAM` and `mirt`.

## Installation

```R
# Install eatRep from GitHub via
remotes::install_github("weirichs/eatModel", upgrade="never")
```

## View package documentation

```R
library(eatModel)
### View package documentation
package?eatModel
```

## Exemplary analysis

The package comes with two examplary data sets, one only for dichotomous items, one also for polytomous items. First data set (`trends`) contains fictional achievement scores of 13524 students in two domains (reading, listening) from three measurement occasions (2010, 2015, 2020) in the long format. The data are not longitudinal as the individuals stem from different cohorts and do not overlap across measurement occasions. The name of the data set is `trends`. You can find detailed exemplary analyses in the help pages of `defineModel()`. The following example estimates a between-item two-dimensional 2pl model with conditioning (background) model for the student cohort measured in 2020. The `splitModels()` function can be used previously if the corresponding model should be estimated in 2010 and 2015 likewise. 

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

## Exemplary partial credit analysis

The second exemplary data set is named `reading` and contains dichotomous as well as polytomous items from a large-scale reading competence assessment. In the following we demonstate a restricted generalized partial credit model (GPCM) using `TAM`. We define arbitrarily some items as dedicated for students with special educational needs (SEN) for which a different slope parameter should be assumed.

```R
# load partial credit long format data
data(reading)

# transform into wide format
datW <- reshape2::dcast(reading[which(reading[,"type"] != "iglu"),],
        uniqueID+sex+language+country~item, value.var = "valueSum")

# combine the two lowest categories for item D205143 to avoid categories with too few observations
datW[,"D205143"] <- car::recode(datW[,"D205143"], "1=0; 2=1; 3=2; 4=3; 5=4")

# restricted generalized partial credit model ... we consider the female subgroup to be the reference population
# we assume that task D223 and D224 are created for students with special educational need and therefore may have
# different slope. Non-SPF items (dichotomous as well as partial credit) should have slope = 1
# data for females
dFema<- datW[which(datW[,"sex"] == "female"),]                               ### reference population mean and SD   
items<- colnames(dFema)[-c(1:4)]

# create an item information data.frame which contain the task an item belongs to, the common slope group,
# and the fixed slope parameter for 1pl items 
items<- data.frame(item = items, task = substr(items,1,4), 
                   slopeGrp = as.numeric(substr(items,1,4) \%in\% c("D223", "D224")), stringsAsFactors = FALSE) |> 
        dplyr::mutate(slope = car::recode(slopeGrp, "0=1; else=NA"))
items[,"dim"] <- car::recode(substr(items[,"item"], 1,2),"'D0'='norm'; 'D2'='pilot'")

# add the group (pilot, norm) each item belongs to 
items<- data.frame ( items, model.matrix(~dim-1, data = items), stringsAsFactors = FALSE)

# define a retsricted generalized partial credit model. Dichotomous items which 
# are constrained to have a slope of 1 are modeled as Rasch items 
def1T<- defineModel(dat=dFema, items = items[,"item"], id=1, irtmodel = "GPCM.groups", 
        software="tam", fixSlopeMat = na.omit(items[,c("item", "slope")]), 
        est.slopegroups = items[,c("item", "slopeGrp")], nodes = 21, 
        fac.oldxsi =0.6, increment.factor=1.05)
run1T<- runModel(def1T)
res1T<- getResults(run1T)
it1T <- itemFromRes(res1T)

# males are focus group: initial free estimation of item parameters for males 
def2T<- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1, 
        irtmodel = "GPCM.groups",software="tam", fixSlopeMat = na.omit(items[,c("item", "slope")]), 
        est.slopegroups = items[,c("item", "slopeGrp")], nodes = 21)
run2T<- runModel(def2T)
res2T<- getResults(run2T)
it2T <- itemFromRes(res2T)

# average discrimination (reg vs. spf) differs for males, but not for females
unique(it1T[,"estSlope"]); unique(it2T[,"estSlope"])                            

# use only slope=1 items for 1pl linking
eq   <- equat1pl(results = res2T, prmNorm = it1T[which(it1T[,"estSlope"] ==1),], 
        item = "item", cat="category", value = "est", difBound=.64, iterativ = TRUE)

# variant 1
# use the DIF-cleaned set of original (female) anchor parameters for calibrat### reference population mean and SDing males on the reference scale
weg  <- eq[["items"]][["not_specified"]][["Dim1"]][["info"]][-1,"itemExcluded"]
weg  <- eatTools::whereAre(weg, paste(it1T[,"item"], it1T[,"category"], sep="_"))
it1C <- subset(it1T[-weg,], estSlope == 1)
def3T<- defineModel(eatTools::na_omit_selection(dat=datW[which(datW[,"sex"] == "male"),],"language"), 
        items = -c(1:4), id=1, irtmodel = "GPCM.groups", anchor = it1C, HG.var = c("language", "country"), 
        itemCol = "item", valueCol = "est", catCol = "category", fixSlopeMat = na.omit(items[,c("item", "slope")]),
        qMatrix = items[,c("item", "dimnorm", "dimpilot")], est.slopegroups = items[,c("item", "slopeGrp")], 
        nodes = 21, software="tam")
run3T<- runModel(def3T)
res3T<- getResults(run3T)
it3T <- itemFromRes(res3T)                                                      

# all items except the ones with linking dif with equal item parameters? check 
comp <- merge(subset(it1T[,c("item", "category", "est", "estSlope")],estSlope == 1), 
              it3T[,c("item", "category", "est", "offset")], by=c("item", "category"), 
              all=TRUE, suffixes = c("_ref", "_foc"))
equal<- na.omit(comp[,c("est_ref", "offset")])
### reference population mean and SD
# all item parameters without linking dif should be equal
stopifnot(all(equal[,1] == equal[,2]))                                          
lDif <- subset(comp, !is.na(est_foc))
stopifnot(all(paste(subset(lDif, !is.na(est_ref))[,"item"], 
          subset(lDif, !is.na(est_ref))[,"category"], sep="_") \%in\% eq$items[["not_specified"]][["Dim1"]][["info"]][,"itemExcluded"]))

# variant 2: use the 1pl item parameters for males (linking dif items excluded), transformed to the metric of females
def4T<- defineModel(eatTools::na_omit_selection(dat=datW[which(datW[,"sex"] == "male"),],"language"), 
        items = -c(1:4), id=1, irtmodel = "GPCM.groups",  
        anchor = eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]],
        HG.var = c("language", "country"), itemCol = "item", valueCol = "est", 
        catCol = "category", fixSlopeMat = na.omit(items[,c("item", "slope")]),
        qMatrix = items[,c("item", "dimnorm", "dimpilot")], 
        est.slopegroups = items[,c("item", "slopeGrp")], nodes = 21, software="tam")
run4T<- runModel(def4T)
res4T<- getResults(run4T)
it4T <- itemFromRes(res4T)                                                      

### all items except the ones with linking dif with equal item parameters? check 
link <- eq[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]][,c("item", "category", "est")]
comp <- merge(link, it4T[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"), all=TRUE)

# all item parameters without linking dif should be equal
equal<- na.omit(comp[,c("est_ref", "offset")])
stopifnot(all(equal[,1] == equal[,2]))                                          

# all items with specific focus paraeter must be included in linking DIF exclusion list 
lDif <- subset(comp, !is.na(est_foc))                                           
stopifnot(all(paste(subset(lDif, !is.na(est_ref))[,"item"], 
          subset(lDif, !is.na(est_ref))[,"category"], sep="_") \%in\% eq$items[["not_specified"]][["Dim1"]][["info"]][,"itemExcluded"]))

# transform to Bista metric
eq4  <- equat1pl(results = res4T)                                               

# reference population mean and SD
refP <- data.frame(domain = c("dimnorm", "dimpilot"), m = 0.0389, sd = 1.07108, stringsAsFactors = FALSE)
cuts <- list ( dimnorm = list(values = 390+0:3*75), dimpilot = list(values = 390+0:3*75))
tf4  <- transformToBista(equatingList=eq4, refPop=refP, cuts=cuts, vera = FALSE)
```

