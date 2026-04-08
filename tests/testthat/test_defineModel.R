### Tests for comparing examples 8 (TAM) and 8.1 (Conquest)

### loading data and edit to fit format ----------------------------------------

# load partial credit long format data
data(reading)
# transform into wide format
datW <- reshape2::dcast(reading[which(reading[,"type"] != "iglu"),],uniqueID+sex+language+country~item, value.var = "valueSum")
# only some items are partial credit, which ones?
pc.it<- reading[which(reading[,"valueSum"] > 1),"item"] |> unique()

# combine the two lowest categories to avoid categories with too few observations
datW[,"D205143"] <- car::recode(datW[,"D205143"], "1=0; 2=1; 3=2; 4=3; 5=4")

# partial credit model ... we consider the female subgroup to be the reference population
# (norm population) that defines the scale and reference item parameters
dFema <- datW[which(datW[,"sex"] == "female"),]                                  ### data for females


################################################################################
###   Example 8: anchored partial credit model excluding linking DIF (TAM)   ###
################################################################################

# partial credit model ... we consider the female subgroup to be the reference population
def1_tam <- defineModel(dat=dFema, items = -c(1:4), id=1, irtmodel = "PCM", software="tam", nodes = 21)
run1_tam <- runModel(def1_tam)
res1_tam <- getResults(run1_tam)
it1_tam  <- itemFromRes(res1_tam)

# males are focus group: initial free estimation of item parameters
def2_tam <- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1, irtmodel = "PCM",software="tam")
run2_tam <- runModel(def2_tam)
res2_tam <- getResults(run2_tam)
it2_tam  <- itemFromRes(res2_tam)

# link males to females ... males perform worse
# 10 items with linking dif identified
eq_tam   <- equat1pl(results = res2_tam, prmNorm = it1_tam, item = "item", cat="category", value = "est", difBound=.64, iterativ = TRUE)

# re-calibrate males with anchoring: use the DIF-cleaned set of anchor parameters for calibrating males on the reference scale
# variant 1: use the original item parameters for females and exclude linking dif items (it1_cleaned)
weg_tam  <- eq_tam[["items"]][["not_specified"]][["Dim1"]][["info"]][-1,"itemExcluded"]
weg_tam  <- eatTools::whereAre(weg_tam, paste(it1_tam[,"item"], it1_tam[,"category"], sep="_"))
it1C_tam <- it1_tam[-weg_tam,]
def3A_tam <- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1, irtmodel = "PCM",
                    anchor = it1C_tam, itemCol = "item", valueCol = "est", catCol = "category", software="tam")
run3A_tam <- runModel(def3A_tam)
res3A_tam <- getResults(run3A_tam)
it3A_tam <- itemFromRes(res3A_tam)                                                      ### all items except the ones with linking dif with equal item parameters? check

comp_tam <- merge(it1_tam[,c("item", "category", "est")], it3A_tam[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"))
equal_tam <- na.omit(comp_tam[,c("est_ref", "offset")])

test_that("all item parameters without linking dif should be equal", {
  expect_equal(equal_tam[,1], equal_tam[,2])
})
lDif <- subset(comp_tam, !is.na(est_foc))                                           ### all items with specific focus parameter must be included in linking DIF exclusion list
stopifnot(all(paste(lDif[,"item"], lDif[,"category"], sep="_") %in% eq_tam$items[["not_specified"]][["Dim1"]][["info"]][,"itemExcluded"]))

# variant 2: use the item parameters for males (linking dif items excluded), transformed to the metric of females
def3B_tam <- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1, irtmodel = "PCM",
                    anchor = eq_tam[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]],
                    itemCol = "item", valueCol = "est", catCol = "category", software="tam")
run3B_tam <- runModel(def3B_tam)
res3B_tam <- getResults(run3B_tam)
it3B_tam <- itemFromRes(res3B_tam)                                                      ### all items except the ones with linking dif with equal item parameters? check
link_tam <- eq_tam[["items"]][["not_specified"]][["Dim1"]][["cleanedLinkItemPars"]][,c("item", "category", "est")]
comp_tam <- merge(link_tam, it3B_tam[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"))
equal_tam <- na.omit(comp_tam[,c("est_ref", "offset")])

test_that("all item parameters without linking dif should be equal", {
  expect_equal(equal_tam[,1], equal_tam[,2])
})
lDif <- subset(comp_tam, !is.na(est_foc))                                           ### all items with specific focus paraeter must be included in linking DIF exclusion list
stopifnot(all(paste(lDif[,"item"], lDif[,"category"], sep="_") %in% eq_tam$items[["not_specified"]][["Dim1"]][["info"]][,"itemExcluded"]))

# transform to Bista metric
eq4_tam  <- equat1pl(results = res3B_tam)                                               ### reference population mean and SD
refP <- data.frame(domain = "Dim1", m = 0.0389, sd = 1.07108, stringsAsFactors = FALSE)
cuts <- list ( Dim1 = list(values = 390+0:3*75))
tf4_tam  <- transformToBista(equatingList=eq4_tam, refPop=refP, cuts=cuts, vera = FALSE)

################################################################################
### Example 8.1: same like example 8 (anchored 1pl PCM), now using Conquest  ###
################################################################################

# partial credit model ... we consider the female subgroup to be the reference population
# (norm population) that defines the scale and reference item parameters
def1_conquest <- defineModel(dat=dFema, items = -c(1:4), id=1, model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest_females", dir=tempdir())
run1_conquest <- runModel(def1_conquest)
res1_conquest <- getResults(run1_conquest) # cannot identify "seed" from cdc file, creates NAs
it1_conquest  <- itemFromRes(res1_conquest)

# males are focus group: initial free estimation of item parameters
def2_conquest <- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1, model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest_males", dir=tempdir())
run2_conquest <- runModel(def2_conquest)
res2_conquest <- getResults(run2_conquest)
it2_conquest  <- itemFromRes(res2_conquest)

# link males to females ... males perform worse
# 11 items with linking dif identified
eq_conquest   <- equat1pl(results = res2_conquest, prmNorm = it1_conquest, item = "item", cat="category", value = "est", difBound=.64, iterativ = TRUE)

# re-calibrate males with anchoring: use the DIF-cleaned set of anchor parameters for calibrating males on the reference scale
# variant 1: use the original item parameters for females and exclude linking dif items (it1_cleaned)
# please note that to date in Conquest only dichotomous items can be used for anchoring
weg_conquest  <- eq_conquest[["items"]][["pcm_conquest_males"]][["Dim1"]][["info"]][-1,"itemExcluded"]
weg_conquest  <- eatTools::whereAre(weg_conquest, paste(it1_conquest[,"item"], it1_conquest[,"category"], sep="_"))
it1C_conquest <- it1_conquest[-weg_conquest,]
pcmIt <- unique(subset(it1C_conquest, category == "Cat2")[,"item"])                       ### additionally exclude pcm items
it1C_conquest <- it1C_conquest[-eatTools::whereAre(pcmIt, it1C_conquest[,"item"]),]
def3A_conquest <- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1,
                    model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest_males_anchored1", dir=tempdir(),
                    anchor = it1C_conquest[,c("item", "est")])
run3A_conquest <- runModel(def3A_conquest)
res3A_conquest <- getResults(run3A_conquest)
it3A_conquest <- itemFromRes(res3A_conquest)                                                      ### all dichotomous items except the ones with linking dif with equal item parameters? check
comp_conquest <- merge(it1_conquest[,c("item", "category", "est")], it3A_conquest[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"))
equal_conquest <- na.omit(comp_conquest[,c("est_ref", "offset")])

test_that("all item parameters without linking dif should be equal", {
  expect_equal(equal_conquest[,1], equal_conquest[,2])
})
### TO DO: all items with specific focus parameter must be included in linking DIF exclusion list or must be pcm items

# variant 2: use the item parameters for males (linking dif items excluded), transformed to the metric of females
ank3B_conquest <- eq_conquest[["items"]][["pcm_conquest_males"]][["Dim1"]][["cleanedLinkItemPars"]]
pcmIt_conquest <- unique(subset(ank3B_conquest, category == "Cat2")[,"item"])                      ### additionally exclude pcm items
ank3B_conquest <- ank3B_conquest[-eatTools::whereAre(pcmIt_conquest, ank3B_conquest[,"item"]),]
def3B_conquest <- defineModel(dat=datW[which(datW[,"sex"] == "male"),], items = -c(1:4), id=1,
                    model.statement = "item+item*step", analysis.name = "pcm_conquest_males_anchored2",
                    dir=tempdir(), anchor = ank3B_conquest[,c("item", "est")])
run3B_conquest <- runModel(def3B_conquest)
res3B_conquest <- getResults(run3B_conquest)
it3B_conquest <- itemFromRes(res3B_conquest)                                                      ### all items except the ones with linking dif with equal item parameters? check
link_conquest <- eq_conquest[["items"]][["pcm_conquest_males"]][["Dim1"]][["cleanedLinkItemPars"]][,c("item", "category", "est")]
comp_conquest <- merge(link_conquest, it3B_conquest[,c("item", "category", "est", "offset")], by=c("item", "category"), suffixes = c("_ref", "_foc"))
equal_conquest <- na.omit(comp_conquest[,c("est_ref", "offset")])

test_that("all item parameters without linking dif should be equal", {
  expect_equal(equal_conquest[,1], equal_conquest[,2])
})
### TO DO: all items with specific focus parameter must be included in linking DIF exclusion list or must be pcm items

# transform to Bista metric
eq4_conquest  <- equat1pl(results = res3B_conquest)                                               ### reference population mean and SD (same lines as in tam model)
#refP <- data.frame(domain = "Dim1", m = 0.0389, sd = 1.07108, stringsAsFactors = FALSE)
#cuts <- list ( Dim1 = list(values = 390+0:3*75))
tf4_conquest  <- transformToBista(equatingList=eq4_conquest, refPop=refP, cuts=cuts, vera = FALSE)

### Tests ----------------------------------------------------------------------

test_that("compare it1", {
  # p-value, Nvalid - exact same
  expect_equal(it1_tam$itemP, it1_conquest$itemP)

  pos <- !is.na(it1_conquest$Nvalid) # positions of rows without NA in `Nvalid`
  expect_equal(it1_tam$Nvalid[pos], it1_conquest$Nvalid[pos]) # conquest has NA

  # est, thurstone - appr. same
  expect_equal(it1_tam$est, it1_conquest$est, tolerance = 0.01)
  expect_equal(it1_tam$thurstone, it1_conquest$thurstone, tolerance = 0.01)
})

test_that("compare it3A", {
  # p-value, Nvalid - exact same
  expect_equal(it3A_tam$itemP, it3A_conquest$itemP)

  pos <- !is.na(it3A_conquest$Nvalid) # positions of rows without NA in `Nvalid`
  expect_equal(it3A_tam$Nvalid[pos], it3A_conquest$Nvalid[pos]) # conquest has NA

  # est, thurstone - appr. same - fix NA problem
  pos <- !is.na(it3A_tam$est & it3A_conquest$est) # positions of rows without NA in `est`
  expect_equal(it3A_tam$est[pos], it3A_conquest$est[pos], tolerance = 0.19) # conquest and tam have NA

  expect_equal(it3A_tam$thurstone, it3A_conquest$thurstone, tolerance = 0.13)
})

test_that("compare it3B", {
  # p-value, Nvalid - exact same
  expect_equal(it3B_tam$itemP, it3B_conquest$itemP)

  pos <- !is.na(it3B_conquest$Nvalid) # positions of rows without NA in `Nvalid`
  expect_equal(it3B_tam$Nvalid[pos], it3B_conquest$Nvalid[pos]) # conquest has NA
  # est, thurstone - appr. same
  pos <- !is.na(it3B_tam$est & it3B_conquest$est) # positions of rows without NA in `est`
  expect_equal(it3B_tam$est[pos], it3B_conquest$est[pos], tolerance = 0.17) # conquest and tam have NA

  expect_equal(it3B_tam$thurstone, it3B_conquest$thurstone, tolerance = 0.15)
})

test_that("compare tf4", {
  # Personpars: plausibel values - appr. same
  expect_equal(tf4_tam$personpars$value, tf4_conquest$personpars$value, tolerance = 1.3)
  # means
  expect_equal(tf4_tam$means[3:6], tf4_conquest$means[3:6], tolerance = 0.01)
})

