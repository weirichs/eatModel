# load partial credit long format data
data(reading)

# create a small wide format exemplary data set, only booklet 0
dw <- reshape2::dcast(subset(reading, bookletID == "TH08"), idstud+sex~item, value.var = "valueSum") |> 
  dplyr::mutate(sex = car::recode(sex, "'male'=0; 'female'=1", as.factor=FALSE)) |>
  eatTools::na_omit_selection("sex")

# define the model: generalized partial credit to get slopes and difficulties. for now, TAM is being used
defT<- defineModel(dat = dw, items = -c(1:2), id = "idstud", irtmodel = "GPCM",software="tam")
runT<- runModel(defT)
resT<- getResults(runT)

# item parameters
it <- itemFromRes(resT)

## the same model in mirt, partially anchored: this is what has to be code reviewed
anchor <- subset(it, item %in% c("D025033", "D204013", "D205143", "D225013"))[,c("item", "category", "est")]
slope <- unique(subset(it, item %in% unique(c("D204013", "D225053", "D025033", "D204013", "D205143", "D225013")))[,c("item", "estSlope")])
items<- colnames(dw)[-c(1:2)]
pcItems <- sapply(items, FUN = function(i) { any(dw[,i] %in% 2)})
items<- data.frame(item = items, irtmod = "2PL")
items[which(items[,"item"] %in% names(pcItems[which(pcItems)])),"irtmod"] <- "gpcm"

## Item parameters we give into the model are in the IRT-parametrisation (with difficulty and discrimination par)
defM<- defineModel(
  dat = dw, 
  items = -c(1:2), 
  id = "idstud",
  irtmodel = items,
  software="mirt", 
  anchor = anchor,  ## where is the discrimination par? 
  itemCol ="item", 
  valueCol = "est", 
  catCol ="category", 
  fixSlopeMat=slope)

## fit the model. In here difficulties should be recalculated to intercepts. Happens in line 166 in mirt.R 
runM<- runModel(defM)
resM<- getResults(runM)
it_mirt <- itemFromRes(resM)

## Anchor items should have been fixed. So they should get the same parameters we gave into the model. 
item_test = "D025033"

## Initial anchoring parameters are found in: 
anchor |> 
  subset(item == item_test)
slope

## output: 
it_mirt |> 
  subset(item == item_test)




coef(runM)
## Hmmn. est is empty for anchors, because it is not estimated there. 
## so maybe better compare the skeleton 

## get the coefficients for the anchor item directly from mirt (withe eatModel prep up front):
coef_mirt <- coef(runM, IRTpars = TRUE)
coef_mirt[[item_test]]

## they should be the same, as when I use the function from mirt to calc slopes and intercept parametrisation from traditional 
## parametrisation, as they are not estimated in the model: 
anch_test <- anchor |> 
  subset(item == item_test)
anch_b <- anch_test[order(anch_test$category), "est"]

## this is IRT parametrisation, so they should be called b1, b2 ...
anchor_vec <- c(
  a = slope[slope$item == item_test, "estSlope"], 
  b1 = anch_vec[1], 
  b2 = anch_vec[2],
  b3 = anch_vec[3],
  b4 = anch_vec[4]
)

## Parametrisierungsvariante: step (d) vs delta (location). 
## d: schwellenwert als Abstand zu mittlerer Schwierigkeit beta. 
## location: Schnittpunkte zwischen Kurven. Das ist das, was wir brauchen. d muss also in delta umgerechnet werden. 
## thurstanian thresholds: ability, wo der Schnittpunkt genau bei 50% liegt, also wo es gleich wahrscheinlich ist, in cat
## 1 oder 2 zu fallen. 
## Wir brauchen wie gesagt delta. Ist das der grund, warum wir das anders berechnen als traditional2mirt? 
expected_mirt <- mirt:::traditional2mirt(anchor_vec, "gpcm", ncat = 5)



#  all D025033  gpcm   a1      9  0.729   -Inf    Inf FALSE  none   none       none     NaN     NaN
# 10   all D025033  gpcm  ak0     10  0.000   -Inf    Inf FALSE  none   none       none     NaN     NaN
# 11   all D025033  gpcm  ak1     11  1.000   -Inf    Inf FALSE  none   none       none     NaN     NaN
# 12   all D025033  gpcm  ak2     12  2.000   -Inf    Inf FALSE  none   none       none     NaN     NaN
# 13   all D025033  gpcm  ak3     13  3.000   -Inf    Inf FALSE  none   none       none     NaN     NaN
# 14   all D025033  gpcm  ak4     14  4.000   -Inf    Inf FALSE  none   none       none     NaN     NaN
# 15   all D025033  gpcm   d0     15  0.000   -Inf    Inf FALSE  none   none       none     NaN     NaN
# 16   all D025033  gpcm   d1     16  1.283   -Inf    Inf FALSE  none   none       none     NaN     NaN
# 17   all D025033  gpcm   d2     17  0.680   -Inf    Inf FALSE  none   none       none     NaN     NaN
# 18   all D025033  gpcm   d3     18  0.968   -Inf    Inf FALSE  none   none       none     NaN     NaN
# 19   all D025033  gpcm   d4     19  0.178   -Inf    Inf FALSE  none   none       none     NaN     NaN



## estimates for the anchor items: This are the ones we want fixed: 
    sub  <- anchor[which(anchor$item == item_test), ]
    b    <- sub[order(sub$category), "est"]
    a    <- slope[which(slope$item == item_test), "estSlope"]



## Now, this should give me the Itemparameters in the slope/intercept form: 
    estimated_eatModel   <- coef(runM, IRTpars = FALSE)[[item_test]]["par", ]

    expect_equal(cf[["a1"]], expected[["a1"]])
    expect_equal(cf[["d1"]], expected[["d1"]])
    expect_equal(cf[["d2"]], expected[["d2"]])
    ## 2) Roundtrip in traditioneller Parametrisierung: die verankerten Schwellen
    ##    und Trennschaerfen muessen unveraendert in itemFromRes(resM) auftauchen
    mrg  <- merge(anchor, it_mirt[, c("item", "category", "est")], by = c("item", "category"), suffixes = c("_ank", "_mirt"))
    expect_equal(nrow(mrg), nrow(anchor))
    expect_equal(mrg[, "est_mirt"], mrg[, "est_ank"], tolerance = 1e-4)
    mrgS <- merge(slope, unique(it_mirt[, c("item", "estSlope")]), by = "item", suffixes = c("_ank", "_mirt"))
    expect_equal(mrgS[, "estSlope_mirt"], mrgS[, "estSlope_ank"], tolerance = 1e-4)




## The problem: mirt needs items in Intercept/Slope form, but Items can also come from the model in difficulty/discrimination parameter form. 
## a: slope/discrimination
## b: difficulty
## d: intercept


# Test run ---------------------------------------------------------------
    data(reading)
    dw   <- reshape2::dcast(subset(reading, bookletID == "TH08"), idstud+sex~item, value.var = "valueSum") |>
            dplyr::mutate(sex = car::recode(sex, "'male'=0; 'female'=1", as.factor=FALSE)) |>
            eatTools::na_omit_selection("sex")
    defT <- defineModel(dat = dw, items = -c(1:2), id = "idstud", irtmodel = "GPCM", software = "tam")
    runT <- runModel(defT)
    it   <- itemFromRes(getResults(runT))
    ankItems <- c("D025033", "D204013", "D205143", "D225013")
    anchor   <- subset(it, item %in% ankItems)[, c("item", "category", "est")]
    slope    <- unique(subset(it, item %in% c(ankItems, "D225053"))[, c("item", "estSlope")])
    pcItems  <- sapply(colnames(dw)[-c(1:2)], FUN = function(i) {any(dw[, i] %in% 2)})
    items    <- data.frame(item = names(pcItems), irtmod = ifelse(pcItems, "gpcm", "2PL"), stringsAsFactors = FALSE)
    defM <- defineModel(dat = dw, items = -c(1:2), id = "idstud", irtmodel = items, software = "mirt",
                        anchor = anchor, itemCol = "item", valueCol = "est", catCol = "category", fixSlopeMat = slope)
    runM <- runModel(defM)
    it_mirt <- itemFromRes(getResults(runM))
    ## 1) interne mirt-Parametrisierung (Intercept/Slope): die fixierten Parameter
    ##    des polytomen Testitems muessen den via traditional2mirt umgerechneten
    ##    Ankerwerten entsprechen, d.h. a1 = a und d_k = -a*(b_1 + ... + b_k)
  