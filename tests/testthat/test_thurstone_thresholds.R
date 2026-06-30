expect_all_abs_lt <- function(x, tolerance) {
    expect_true(all(abs(x) < tolerance), info = paste("Maximum absolute difference:", max(abs(x), na.rm = TRUE)))
}

### only 1pl Rasch and partial credit
data(reading)
#sysInfo  <- Sys.info()
cf <- system.file("extdata", "console_Feb2007.exe", package = "eatModel")
hasConquest <- nchar(cf)>0 && file.exists(cf)

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

# Unfortunately, the parameters are not exactly the same, but only approximately so.
# This is inherent in probabilistic estimation. The test is therefore also performed
# only on the basis of approximate equality.
itT <- itemFromRes(resT)
itM <- itemFromRes(resM)
merge1 <- eatTools::mergeAttr(itT[,c("item", "category", "est", "thurstone")], itM[,c("item", "category", "est", "thurstone")], by=c("item", "category"), setAttr=FALSE, all=TRUE, xName="tam", yName = "mirt", suffixes = c("_tam", "_mirt"))
merge1[,"diff_est"] <- merge1[,"est_tam"] - merge1[,"est_mirt"]
merge1[,"diff_thurs"] <- merge1[,"thurstone_tam"] - merge1[,"thurstone_mirt"]
ind <- which(abs(merge1[,"diff_est"]) < 0.01)

# tam.threshold() is used for TAM, while mirt thresholds are computed in eatModel.
# For PCM, both paths should lead to nearly identical 62.5% Thurstonian thresholds.
test_that("PCM Thurstonian thresholds are comparable for TAM and mirt", {
   expect_gt(nrow(merge1), 0)
   expect_gt(length(ind), 0)
   expect_all_abs_lt(merge1[ind,"diff_thurs"], tolerance = 0.02)
})

# conquest
if(hasConquest){
    defC<- suppressWarnings(defineModel(dat = dw, items = -1, id = "idstud", model.statement = "item+item*step", nodes = 21, analysis.name = "pcm_conquest", dir=tempdir()))
    runC<- runModel(defC)
    resC<- suppressWarnings(getResults(runC))

    # tam vs. conquest
    itC <- itemFromRes(resC)
    merge2 <- eatTools::mergeAttr(itT[,c("item", "category", "est", "thurstone")], itC[,c("item", "category", "est", "thurstone")], by=c("item", "category"), setAttr=FALSE, all=TRUE, xName="tam", yName = "conquest", suffixes = c("_tam", "_conquest"))
    merge2[,"diff_est"] <- merge2[,"est_tam"] - merge2[,"est_conquest"]
    merge2[,"diff_thurs"] <- merge2[,"thurstone_tam"] - merge2[,"thurstone_conquest"]

    # ConQuest is optional on developer machines. If available, the threshold
    # extraction should remain in the same range as the TAM output.
    test_that("PCM Thurstonian thresholds are comparable for TAM and ConQuest", {
       expect_all_abs_lt(merge2[,"diff_thurs"], tolerance = 0.05)
    })
}

test_that("GPCM Thurstonian thresholds and BISTA item parameters are comparable for TAM and mirt", {
   data(data.gpcm, package = "TAM")
   gpcmDat <- data.frame(pid = seq_len(nrow(data.gpcm)), data.gpcm)

   # TAM uses TAM::tam.threshold(prob.lvl = .625) for both PCM and GPCM. The
   # mirt path solves the same cumulative GPCM probability directly in eatModel.
   defTG <- suppressWarnings(defineModel(dat = gpcmDat, items = -1, id = "pid", irtmodel = "GPCM",
                                         software = "tam", n.iterations = 300, nodes = 15,
                                         boundary = 1, remove.no.answers = FALSE, progress = FALSE,
                                         analysis.name = "gpcm_tam", verbose = FALSE))
   runTG <- suppressWarnings(runModel(defTG))
   resTG <- suppressWarnings(getResults(runTG, omitPV = TRUE, omitWle = TRUE, omitFit = TRUE, omitRegr = TRUE))

   irtMG <- data.frame(item = colnames(gpcmDat)[-1], irtmod = "gpcm", stringsAsFactors = FALSE)
   defMG <- suppressWarnings(defineModel(dat = gpcmDat, items = -1, id = "pid", irtmodel = irtMG,
                                         software = "mirt", n.iterations = 300, nodes = 15,
                                         boundary = 1, remove.no.answers = FALSE, progress = FALSE,
                                         analysis.name = "gpcm_mirt", verbose = FALSE))
   runMG <- suppressWarnings(runModel(defMG))
   resMG <- suppressWarnings(getResults(runMG, omitPV = TRUE, omitWle = TRUE, omitFit = TRUE, omitRegr = TRUE))

   itTG <- itemFromRes(resTG)
   itMG <- itemFromRes(resMG)
   mergeG <- eatTools::mergeAttr(itTG[,c("item", "category", "estSlope", "thurstone")],
                                 itMG[,c("item", "category", "estSlope", "thurstone")],
                                 by = c("item", "category"), setAttr = FALSE, all = TRUE,
                                 xName = "tam", yName = "mirt", suffixes = c("_tam", "_mirt"))
   mergeG[,"diff_thurs"] <- mergeG[,"thurstone_tam"] - mergeG[,"thurstone_mirt"]

   expect_gt(nrow(mergeG), 0)
   expect_false(anyNA(mergeG[,c("thurstone_tam", "thurstone_mirt")]))
   expect_all_abs_lt(mergeG[,"diff_thurs"], tolerance = 0.02)

   # transformToBista() should use these Thurstonian thresholds for GPCM. With a
   # common theta reference population, TAM and mirt should therefore produce
   # almost identical BISTA item parameters.
   refPop <- data.frame(domain = "Dim1", m = 0, sd = 1)
   cuts <- list(Dim1 = list(values = c(400, 500, 600)))
   tfTG <- suppressWarnings(transformToBista(equat1pl(resTG), refPop = refPop, cuts = cuts, vera = FALSE))
   tfMG <- suppressWarnings(transformToBista(equat1pl(resMG), refPop = refPop, cuts = cuts, vera = FALSE))
   bistaG <- eatTools::mergeAttr(tfTG[["itempars"]][,c("item", "category", "dimension", "estTransf625", "estTransfBista")],
                                 tfMG[["itempars"]][,c("item", "category", "dimension", "estTransf625", "estTransfBista")],
                                 by = c("item", "category", "dimension"), setAttr = FALSE, all = TRUE,
                                 xName = "tam", yName = "mirt", suffixes = c("_tam", "_mirt"))
   bistaG[,"diff_theta"] <- bistaG[,"estTransf625_tam"] - bistaG[,"estTransf625_mirt"]
   bistaG[,"diff_bista"] <- bistaG[,"estTransfBista_tam"] - bistaG[,"estTransfBista_mirt"]

   expect_gt(nrow(bistaG), 0)
   expect_false(anyNA(bistaG[,c("estTransf625_tam", "estTransf625_mirt", "estTransfBista_tam", "estTransfBista_mirt")]))
   expect_all_abs_lt(bistaG[,"diff_theta"], tolerance = 0.02)
   expect_all_abs_lt(bistaG[,"diff_bista"], tolerance = 2)
})











