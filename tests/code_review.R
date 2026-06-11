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

# add model information to 'item'

items<- data.frame(item = items, irtmod = "2PL")
items[which(items[,"item"] %in% names(pcItems[which(pcItems)])),"irtmod"] <- "gpcm"
defM<- defineModel(dat = dw, items = -c(1:2), id = "idstud",irtmodel = items,software="mirt", anchor = anchor, itemCol ="item", valueCol = "est", catCol ="category", fixSlopeMat=slope)

runM<- runModel(defM)
resM<- getResults(runM)
it_mirt <- itemFromRes(resM)


item_test = "D204013"

## mirt
coef_mirt <- coef(runM, IRTpars = FALSE)
coef_mirt[[item_test]]

## eatModel
it_mirt |> 
  subset(item == item_test)


## The problem: mirt needs items in Intercept/Slope form, but Items can also come from the model in difficulty/discrimination parameter form. 
## a: slope/discrimination
## b: difficulty
## d: intercept

