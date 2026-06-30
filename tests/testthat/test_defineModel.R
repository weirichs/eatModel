cf <- system.file("extdata", "console_Feb2007.exe", package = "eatModel")
hasConquest <- nchar(cf)>0 && file.exists(cf)
path <- getwd() # defineMOdel aender setwd(), deshalb muss das hier wieder zurueckgesetzt werden

#load(pl2w("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/vera_dat.rda"))
load(test_path("vera_dat.rda"))

if(hasConquest){
    test_that("defineModel works for DIF", {
        itemIndex <- grep("^CE|^EC|^CO|^OC", colnames(findat))
        itemNames <- colnames(findat)[itemIndex]
        qMat     <- data.frame ( item = itemNames, dim = car::recode (substr(itemNames, 1, 2 ), "'CO'='orale'; 'OC'='orale' ; 'CE'='ecrite' ; 'EC'='ecrite'"), stringsAsFactors=FALSE)
        qMat     <- data.frame ( qMat, model.matrix ( ~dim-1, data = qMat ) )
        findat <- findat[-which(findat$Geschl %in% c("mbi","mir", "mbd", "mci", "3")),]
        findat[,"Geschl"] <- car::recode(findat$Geschl, "1=0; 2=1")
        modSplit <- splitModels ( qMatrix = qMat [ c("item", "dimorale", "dimecrite")], nCores = 1 )
        # damit testthat output leichter lesbar, werden warnungen hier ausgeschaltet ... fuer die user sollen sie aber drinbleiben
        modsDef_Geschl_Conquest  <- suppressWarnings(defineModel(dat = findat, splittedModels = modSplit,id = "thnum", HG.var = "thdiff", remove.failures = TRUE, DIF.var = "Geschl", dir = tempdir()))
    })
}

###
load(system.file("extdata", "pcm", "input.rda", package = "eatModel"))

if(hasConquest){
    test_that("defineModel works for three-dimensional model", {
        define3d <- suppressWarnings(defineModel(dat = datAufbereitet, items = qmatrix[,"item_neu"], id = "Pseudonyme",HG.var = "Note.y", model.statement = "item + item*step", nodes = 50,
                    qMatrix = qmatrix[,c("item_neu", "dimensionIn", "dimensionSp","dimensionSt")],  dir = tempdir(), analysis.name = "pcm1_3d"))
    })
}
file.zip <- system.file("extdata", "pcm", "results.zip", package = "eatModel")
utils::unzip(zipfile = file.zip, exdir=tempdir())

test_that("importing results works for three-dimensional model", {
    run3d[["dir"]] <- tempdir()
    results3d <- suppressWarnings(getResults(run3d))
})

# strange item namings ...
# load(pl2w("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/japan.rda"))
setwd(path)
load(test_path("japan.rda"))

if(hasConquest){
    test_that("temporal item renaming works", {
        defDif3 <- suppressWarnings(defineModel(dat = jpn, id = "person", items = grep("^M1", colnames(jpn), value=TRUE), software="conquest", DIF.var = "SEX", dir = tempdir(), analysis.name = "DIF"))
        runDif3 <- runModel(defDif3)
        resultsDif3 <- suppressWarnings(getResults(runDif3))
        itemsDif3   <- itemFromRes(resultsDif3)
    })
}
# example 2
setwd(path)
load(test_path("df_example2.rda"))

if(hasConquest){
    test_that("defineModel works for achoring", {
        def <- suppressWarnings(defineModel(dat = datSel, items = intersect(qmat[,"item"],colnames(datSel)), compute.fit = FALSE, id = "IDSTUD", remove.constant.items = TRUE,method = "montecarlo", nodes = 8000,weight.var = "totwgt_NAW",
            qMatrix = qmat[,c("item", "splitbio_knowledge", "splitbio_procedural")], anchor = anchor[["itempars"]][,c("item", "estTransf")], HG.var = c("designDichotom", grep("^PC", colnames(datSel),value=TRUE)), analysis.name = "withWeights", n.plausible = 15, dir = tempdir()))
    })
}
# lange itemnamen conquest, soll schnell durchlaufen
# single
   data(trends)

if(hasConquest){
    test_that("item renaming (back and forward) works for single model", {
    dat  <- subset(trends, year == 2010 & booklet == "Bo01" & domain == "listening")
    datW <- reshape2::dcast(dat,idstud+sex~item, value.var = "value")
    datW[,"sex"]   <- car::recode(datW[,"sex"], "'male'=0; 'female'=1")
    colnames(datW) <- car::recode(colnames(datW), "'T12_07'='T12_07_zu_lang'")
    datW[,"a11"]   <- rnorm(n = nrow(datW), 0, 1.5)
    defDif4 <- suppressWarnings(defineModel(dat = datW, id = "idstud", items = grep("^T", colnames(datW), value=TRUE), software="conquest", DIF.var = "sex", dir = tempdir(), analysis.name = "DIF"))
    runDif4 <- runModel(defDif4)
    resultsDif4 <- suppressWarnings(getResults(runDif4))
    itemsDif4   <- itemFromRes(resultsDif4)
    })
}
# multiple
if(hasConquest){
    test_that("item renaming (back and forward) works for multiple models", {
    dat  <- subset(trends, year == 2010 & booklet == "Bo01")
    datW <- reshape2::dcast(dat,idstud+sex~item, value.var = "value")
    datW[,"sex"]   <- car::recode(datW[,"sex"], "'male'=0; 'female'=1", as.factor=FALSE)
    colnames(datW) <- car::recode(colnames(datW), "'T12_07'='T12_07_zu_lang'; 'T02_01'='T02_01_auch_zu_lang'")
    datW[,"a11"]   <- rnorm(n = nrow(datW), 0, 1.5)
    qma  <- unique(dat[,c("item", "domain")])
    qma  <- data.frame(qma, model.matrix(~domain-1,data=qma), stringsAsFactors = FALSE)
    split<- splitModels( qMatrix = qma[,c("item", "domainlistening", "domainreading")], nCores = 1)
    defDif4 <- suppressWarnings(defineModel(dat = datW, id = "idstud", splittedModels = split, software="conquest", DIF.var = "sex", dir = tempdir(), analysis.name = "DIF"))
    runDif4 <- runModel(defDif4)
    resultsDif4 <- suppressWarnings(getResults(runDif4))
    itemsDif4   <- itemFromRes(resultsDif4)
    })
}
