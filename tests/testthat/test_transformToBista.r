compareObj <- function (obj1, obj2, by) {
  vgl <- merge(obj1, obj2, by = by, all=TRUE, suffixes = c("_x", "_y"))
  pre <- unique(eatTools::halveString(grep("_x$", colnames(vgl), value=TRUE),"_", first=FALSE)[,1])
  ret <- unique(unlist(lapply(pre, FUN = function (i) {
         cols <- grep(paste0("^", i, "_"), colnames(vgl), value=TRUE)
         stopifnot(length(cols)==2)
         vgl.i<- vgl[,cols]
         return(expect_true(all.equal(vgl.i[,1],vgl.i[,2])))})))
  return(ret)}

# load(pl2w("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/tf_example1.rda"))
load(test_path("tf_example1.rda"))

test_that("tf1", {
  tf2b.neu <- suppressWarnings(transformToBista ( equatingList = eq2b, refPop = ref, cuts = cuts))
  expect_true(compareObj( tf2b$itempars, tf2b.neu$itempars, by=c("item", "dimension")))
  expect_true(all.equal(tf2b$personpars , tf2b.neu$personpars))
  expect_true(all.equal(tf2b$linkingErrors , tf2b.neu$linkingErrors))
  expect_true(compareObj(tf2b$means,tf2b.neu$means, by = c("model", "domain")))
  expect_true(compareObj(tf2b$refPop , tf2b.neu$refPop , by="domain"))
})

# load(pl2w("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/tf_example2.rda"))
   load(test_path("tf_example2.rda"))

test_that("tf2", {
  dfr.neu   <- suppressWarnings(transformToBista ( equatingList = anch3, refPop = refPop, cuts=cuts , vera=FALSE))
  expect_true(compareObj( dfr$itempars, dfr.neu$itempars, by=c("item", "dimension", "model")))
  expect_true(all.equal(dfr$personpars , dfr.neu$personpars))
  expect_true(all.equal(dfr$linkingErrors , dfr.neu$linkingErrors))
  expect_true(compareObj(dfr$means, dfr.neu$means, by = c("model", "domain")))
  expect_true(compareObj(dfr$refPop , unique(dfr.neu$refPop[,-2]) , by="domain"))
})

# load(pl2w("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/tf_example3.rda"))
load(test_path("tf_example3.rda"))
test_that("tf3", {
  dfr.neu3  <- transformToBista ( equatingList = anch3, refPop = refPop, cuts=cuts )
  expect_true(compareObj( dfr$itempars, dfr.neu3$itempars, by=c("item", "dimension", "model")))
  expect_true(all.equal(dfr$personpars , dfr.neu3$personpars))
  expect_true(all.equal(dfr$linkingErrors , dfr.neu3$linkingErrors))
  expect_true(compareObj(dfr$refPop , unique(dfr.neu3$refPop[,-2]) , by="domain"))
})


# load(pl2w("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/tf_example3a.rda"))
# load(pl2w("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/tf_example3b.rda"))
   load(test_path("tf_example3a.rda"))
   load(test_path("tf_example3b.rda"))
test_that("tf3ab", {
  tfRef.neu <- suppressWarnings(transformToBista ( equatingList = eqRef, cuts=cuts, vera=FALSE,weights = datT1[,c("idstud", "wgt")] ))
  expect_true(compareObj( tfRef$itempars, tfRef.neu$itempars, by=c("item", "dimension", "model")))
  expect_true(all.equal(tfRef$personpars , tfRef.neu$personpars))
  expect_true(all.equal(tfRef$linkingErrors , tfRef.neu$linkingErrors))
  expect_true(compareObj(tfRef$refPop , tfRef.neu$refPop, by="domain"))
  expect_true(compareObj(tfRef$means , tfRef.neu$means , by=c("model", "domain")))
})

# load(pl2w("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/tf_example4.rda"))
   load(test_path("tf_example4.rda"))
test_that("tf4", {
  dfrT1P.neu <- suppressWarnings(transformToBista ( equatingList = ankT1P, refPop = tfRef[["refPop"]],cuts=cuts, vera=FALSE ))
  expect_true(compareObj(dfrT1P$itempars, dfrT1P.neu$itempars, by=c("item", "dimension", "model")))
  expect_true(all.equal(dfrT1P$personpars, dfrT1P.neu$personpars))
  expect_true(all.equal(dfrT1P$linkingErrors , dfrT1P.neu$linkingErrors))
  expect_true(compareObj(dfrT1P$means, dfrT1P.neu$means, by=c("model","domain")))
  expect_true(compareObj(dfrT1P$refPop, dfrT1P.neu$refPop, by="domain"))
dfrT2P.neu<- suppressWarnings(transformToBista ( equatingList = ankT2P, refPop=tfRef[["refPop"]], cuts=cuts, vera=FALSE))
  expect_true(compareObj(dfrT2P$itempars, dfrT2P.neu$itempars, by=c("item", "dimension", "model")))
  expect_true(all.equal(dfrT2P$personpars, dfrT2P.neu$personpars))
  expect_true(all.equal(dfrT2P$linkingErrors , dfrT2P.neu$linkingErrors))
  expect_true(compareObj(dfrT2P$means, dfrT2P.neu$means, by=c("model","domain")))
  expect_true(compareObj(dfrT2P$refPop, dfrT2P.neu$refPop, by="domain"))
dfrT3P.neu<- suppressWarnings(transformToBista ( equatingList = ankT3P, refPop=tfRef[["refPop"]], cuts=cuts, vera=FALSE))
  expect_true(compareObj(dfrT3P$itempars, dfrT3P.neu$itempars, by=c("item", "dimension", "model")))
  expect_true(all.equal(dfrT3P$personpars, dfrT3P.neu$personpars))
  expect_true(all.equal(dfrT3P$linkingErrors , dfrT3P.neu$linkingErrors))
  expect_true(compareObj(dfrT3P$means, dfrT3P.neu$means, by=c("model","domain")))
  expect_true(compareObj(dfrT3P$refPop, dfrT3P.neu$refPop, by="domain"))
tle.neu   <- suppressWarnings(transformToBista ( equatingList = L.t1t3, refPop=tfRef[["refPop"]], cuts = cuts, vera=FALSE, years = c(2010,2020)))
  expect_true(compareObj(tle$itempars, tle.neu$itempars, by=c("item", "dimension", "model")))
  expect_true(all.equal(tle$personpars, tle.neu$personpars))
  expect_true(all.equal(tle$linkingErrors , tle.neu$linkingErrors))
  expect_true(compareObj(tle$means, tle.neu$means, by=c("model","domain")))
  expect_true(compareObj(tle$refPop, tle.neu$refPop, by="domain"))
})

test_that("2PL transformation to 62.5% respects software parameterization", {
  makeEq <- function(source) {
    items <- c("I1", "I2")
    est <- c(-0.6, 0.8)
    slope <- c(0.75, 1.5)
    rows <- rbind(
      data.frame(model = "m1", source = source, var1 = items, var2 = "Cat1", type = "fixed", indicator.group = "items", group = "Dim1", par = "est", derived.par = NA, value = est),
      data.frame(model = "m1", source = source, var1 = items, var2 = "Cat1", type = "fixed", indicator.group = "items", group = "Dim1", par = "estSlope", derived.par = NA, value = slope),
      # Include slope standard errors because itemFromRes() creates estSlope_se;
      # transformToBista() must still select the actual slope column.
      data.frame(model = "m1", source = source, var1 = items, var2 = "Cat1", type = "fixed", indicator.group = "items", group = "Dim1", par = "estSlope", derived.par = "se", value = c(0.1, 0.2)),
      data.frame(model = "m1", source = source, var1 = c("p1", "p2"), var2 = NA, type = "fixed", indicator.group = "persons", group = "Dim1", par = "pv", derived.par = "pv1", value = c(-0.2, 0.3)),
      data.frame(model = "m1", source = source, var1 = NA, var2 = NA, type = "tech", indicator.group = NA, group = NA, par = "ID", derived.par = "id", value = NA)
    )
    ret <- list(
      items = list(m1 = list(Dim1 = list(eq = list(B.est = c(Mean.Mean = 0, Haebara = 0, Stocking.Lord = 0), descriptives = c(N.Items = 0, SD = NA, Var = NA, linkerror = NA)), items = NULL, method = "Mean.Mean"))),
      results = rows
    )
    class(ret) <- c("eq2tom", "list")
    ret
  }
  refPop <- data.frame(domain = "Dim1", m = 0, sd = 1)
  cuts <- list(Dim1 = list(values = c(400, 500, 600)))
  logit625 <- log(0.625/(1 - 0.625))

  # TAM stores the 2PL difficulty in the slope-scaled logit, so the whole
  # threshold expression has to be divided by the slope.
  tam <- suppressWarnings(transformToBista(makeEq("tam"), refPop = refPop, cuts = cuts, vera = FALSE, idVarName = "id"))
  expect_equal(tam$itempars$estTransf625, with(tam$itempars, (estTransf + logit625) / estSlope), tolerance = 1e-10)
  pTam <- with(tam$itempars, plogis(estSlope * estTransf625 - estTransf))
  expect_equal(pTam, rep(0.625, 2), tolerance = 1e-10)

  # mirt stores the 2PL difficulty on the theta scale; only the .625 logit
  # offset is divided by the slope.
  mirt <- suppressWarnings(transformToBista(makeEq("mirt"), refPop = refPop, cuts = cuts, vera = FALSE, idVarName = "id"))
  expect_equal(mirt$itempars$estTransf625, with(mirt$itempars, estTransf + logit625 / estSlope), tolerance = 1e-10)
  pMirt <- with(mirt$itempars, plogis(estSlope * (estTransf625 - estTransf)))
  expect_equal(pMirt, rep(0.625, 2), tolerance = 1e-10)
})
