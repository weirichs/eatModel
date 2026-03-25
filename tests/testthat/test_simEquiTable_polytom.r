### only 1pl Rasch and partial credit
data(reading)

# tam: To speed up the process, the tests will be administered only for Test
# Booklet 8. This booklet contains questions with some partial credit items.
d   <- subset(reading, bookletID == "TH08" & taskID == "D205")
dw  <- reshape2::dcast(d, uniqueID~item, value.var = "valueSum")
defT<- defineModel(dat = dw, items = -1, id = "uniqueID",  irtmodel = "PCM",software="tam")
runT<- runModel(defT)
resT<- getResults(runT, omitPV=TRUE, Q3 = FALSE)
it  <- itemFromRes(resT)

# conventional (cumbersome) equivalence table
perm<- lapply(colnames(dw)[-1], FUN = function(i) {sort(unique(na.omit(dw[,i])))}) |> setNames(colnames(dw)[-1])
simD<- expand.grid(perm)
def2<- simD |> dplyr::mutate(id = paste0("P", 1:nrow(simD))) |> defineModel(items = colnames(simD), id = "id",  irtmodel = "PCM",software="tam", anchor=it, itemCol = "item", valueCol = "est", catCol = "category")
run2<- runModel(def2)
res2<- getResults(run2, omitPV=TRUE)
wle2<- wleFromRes(res2)
eqt2<- unique(wle2[,c("NitemsSolved", "NitemsTotal", "wle_est")])

test_that("there are exactly as many unique WLEs as there are unique sum values", {
  expect_true(length(unique(eqt2[,"NitemsSolved"])) == nrow(eqt2))
  expect_true(length(unique(eqt2[,"wle_est"])) == nrow(eqt2))
})

test_that("there are exactly as many unique WLEs as there are overall item categories", {
  expect_true(nrow(eqt2) == (nrow(it) + 1))
})

# hotfix equivalence table: Achtung! kommt nicht dasselbe raus!
it[,"newItem"] <- paste0(it[,"item"], it[,"category"], sep="_")
et  <- simEquiTable  ( anchor=it[,c("newItem", "est")], mRef=0.05, sdRef=1.2, addConst = 500, multConst = 100, cutScores = list(Dim1=list(values = c(400, 600))))
comp<- merge(eqt2, et[["complete"]], by.x = "NitemsSolved", by.y = "score")

# A more efficient version for generating the equivalence table 'eqt2'
items <- by(it, INDICES = it[,"item"], FUN = function(i) {sort(unique(c(0, as.numeric(eatTools::removeNonNumeric(i[,"category"])))))})

# erzeuge dichotome pseudo-items, und zwar so viele wie es maximale summenwerte gibt
items1<- sum(unlist(lapply(items, max)))
dtmp  <- data.frame(rbind(1*(lower.tri(matrix(1, nrow = items1, ncol = items1))),1))
pcitem<- items[which(sapply(items, length)>2)]
diItem<- setdiff(names(items), names(pcitem))
for (i in names(pcitem)) {
     ncat <- length(pcitem[[i]])-1
     dtmp <- data.frame(rowSums(dtmp[,(ncol(dtmp)-ncat+1):ncol(dtmp)]), dtmp[,1:(ncol(dtmp)-ncat)], stringsAsFactors = FALSE)
     colnames(dtmp)[1] <- i}
colnames(dtmp)[(1+length(pcitem)):ncol(dtmp)] <- diItem

# 'eqt3' created now should be equivalent to 'eqt2'
def3<- dtmp |> dplyr::mutate(id = paste0("P", 1:nrow(dtmp))) |> defineModel(items = colnames(dtmp), id = "id",  irtmodel = "PCM",software="tam", anchor=it, itemCol = "item", valueCol = "est", catCol = "category", verbose=FALSE)
run3<- runModel(def3)
res3<- getResults(run3, omitPV=TRUE)
wle3<- wleFromRes(res3)
comp<- merge(wle3[,c("NitemsSolved", "wle_est")], eqt2[,c("NitemsSolved", "wle_est")], by = "NitemsSolved")

test_that("comparison object has exactly 15 rows", {
  expect_true(nrow(comp) == 15)
})
test_that("both wle columns are equivaent", {
  expect_true(all(comp[,"wle_est.x"] == comp[,"wle_est.y"]))
})

# creating the same table, now using the function simEquiTable()
eqt4 <- simEquiTable( anchor = it, item = "item", cat="category", value = "est",
       cutScores = list ( values = c(400, 600)), mRef = 0.047, sdRef = 1.181)
comp <- merge(eqt4[["complete"]][,c("score", "est")], wle3[,c("NitemsSolved", "wle_est")], by.x = "score", by.y = "NitemsSolved")

test_that("there are exactly as many unique WLEs as there are overall item categories", {
  expect_true(nrow(comp) == (nrow(it) + 1))
})
test_that("both wle columns are equivaent", {
  expect_true(all(comp[,"est"] == comp[,"wle_est"]))
})

