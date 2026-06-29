### only 1pl Rasch and partial credit
data(reading)
#sysInfo  <- Sys.info()
cf <- system.file("extdata", "console_Feb2007.exe", package = "eatModel")

# tam: To speed up the process, the tests will be administered only for Test
# Booklet 8. This booklet contains questions with some partial credit items.
d   <- subset(reading, bookletID == "TH08")
dw  <- reshape2::dcast(d, idstud~item, value.var = "valueSum")
defT<- suppressWarnings(defineModel(dat = dw, items = -1, id = "idstud",  irtmodel = "PCM",software="tam"))
runT<- runModel(defT)
resT<- suppressWarnings(getResults(runT))

# mirt: It's not very intuitive, but partial credit is also defined here using `irtmod="Rasch"`.
# IRT automatically recognizes that it's partial credit when the items are polytomous rather than dichotomous.
irtM<- data.frame(item = colnames(dw)[-1], irtmod = "Rasch", stringsAsFactors = FALSE)
defM<- suppressWarnings(defineModel(dat = dw, items = -1, id = "idstud",  irtmodel = irtM,software="mirt"))
runM<- runModel(defM)
resM<- suppressWarnings(getResults(runM))

# conquest
if(isTRUE(file.exists(cf))){
    defC<- suppressWarnings(defineModel(dat = dw, items = -1, id = "idstud", model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest", dir=tempdir()))
    runC<- runModel(defC)
    resC<- suppressWarnings(getResults(runC))

    # Unfortunately, the parameters are not exactly the same, but only approximately so.
    # This is inherent in probabilistic estimation. The test is therefore also performed
    # only on the basis of approximate equality.
    itT <- itemFromRes(resT)
    itM <- itemFromRes(resM)
    merge1 <- eatTools::mergeAttr(itT[,c("item", "category", "est", "thurstone")], itM[,c("item", "category", "est", "thurstone")], by=c("item", "category"), setAttr=FALSE, all=TRUE, xName="tam", yName = "mirt", suffixes = c("_tam", "_mirt"))
    merge1[,"diff_est"] <- merge1[,"est_tam"] - merge1[,"est_mirt"]
    merge1[,"diff_thurs"] <- merge1[,"thurstone_tam"] - merge1[,"thurstone_mirt"]
    ind <- which(abs(merge1[,"diff_est"]) < 0.01)

    # tam vs. mirt
    test_that("tf1", {
       expect_true(all(abs(merge1[ind,"diff_thurs"]) < 0.02) )
    })

    # tam vs. conquest
    itC <- itemFromRes(resC)
    merge2 <- eatTools::mergeAttr(itT[,c("item", "category", "est", "thurstone")], itC[,c("item", "category", "est", "thurstone")], by=c("item", "category"), setAttr=FALSE, all=TRUE, xName="tam", yName = "conquest", suffixes = c("_tam", "_conquest"))
    merge2[,"diff_est"] <- merge2[,"est_tam"] - merge2[,"est_conquest"]
    merge2[,"diff_thurs"] <- merge2[,"thurstone_tam"] - merge2[,"thurstone_conquest"]

    # tam vs. conquest
    test_that("tf2", {
    expect_true(all(abs(merge2[,"diff_thurs"]) < 0.05) )
    })
}

# to do: jetzt noch fuer 2pl











