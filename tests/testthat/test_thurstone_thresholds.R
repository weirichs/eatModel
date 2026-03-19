### only 1pl Rasch and partial credit
data(reading)

# tam: To speed up the process, the tests will be administered only for Test
# Booklet 8. This booklet contains questions with some partial credit items.
d   <- subset(reading, bookletID == "TH08")
dw  <- reshape2::dcast(d, uniqueID~item, value.var = "valueSum")
defT<- defineModel(dat = dw, items = -1, id = "uniqueID",  irtmodel = "PCM",software="tam")
runT<- runModel(defT)
resT<- getResults(runT)

# mirt: It's not very intuitive, but partial credit is also defined here using `irtmod="Rasch"`.
# IRT automatically recognizes that it's partial credit when the items are polytomous rather than dichotomous.
irtM<- data.frame(item = colnames(dw)[-1], irtmod = "Rasch", stringsAsFactors = FALSE)
defM<- defineModel(dat = dw, items = -1, id = "uniqueID",  irtmodel = irtM,software="mirt")
runM<- runModel(defM)
resM<- getResults(runM)

# conquest: To ensure that the tests run on all devices, Conquest is currently being tested only
# on Windows. Conquest does not (yet) run on Linux, so that will have to be left out for now.
sysInfo  <- Sys.info()
if(sysInfo[["sysname"]] != "Linux") {
   defC<- defineModel(dat = dw, items = -1, id = "uniqueID", model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest", dir=tempdir())
   runC<- runModel(defC)
   resC<- getResults(runC)
}

# Unfortunately, the parameters are not exactly the same, but only approximately so.
# This is inherent in probabilistic estimation. The test is therefore also performed
# only on the basis of approximate equality.
itT <- itemFromRes(resT)
itM <- itemFromRes(resM)
merge1 <- mergeAttr(itT[,c("item", "category", "est", "thurstone")], itM[,c("item", "category", "est", "thurstone")], by=c("item", "category"), setAttr=FALSE, all=TRUE, xName="tam", yName = "mirt", suffixes = c("_tam", "_mirt"))
merge1[,"diff_est"] <- merge1[,"est_tam"] - merge1[,"est_mirt"]
merge1[,"diff_thurs"] <- merge1[,"thurstone_tam"] - merge1[,"thurstone_mirt"]
ind <- which(abs(merge1[,"diff_est"]) < 0.01)

# tam vs. mirt
test_that("tf1", {
  expect_true(all(abs(merge1[ind,"diff_thurs"]) < 0.02) )
})

# tam vs. conquest
if(sysInfo[["sysname"]] != "Linux") {
  itC <- itemFromRes(resC)
  merge2 <- mergeAttr(itT[,c("item", "category", "est", "thurstone")], itC[,c("item", "category", "est", "thurstone")], by=c("item", "category"), setAttr=FALSE, all=TRUE, xName="tam", yName = "conquest", suffixes = c("_tam", "_conquest"))
  merge2[,"diff_est"] <- merge2[,"est_tam"] - merge2[,"est_conquest"]
  merge2[,"diff_thurs"] <- merge2[,"thurstone_tam"] - merge2[,"thurstone_conquest"]

  # tam vs. conquest
  test_that("tf2", {
  expect_true(all(abs(merge2[,"diff_thurs"]) < 0.05) )
  })
}











