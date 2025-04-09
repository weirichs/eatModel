# ref1 <- readRDS("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/lErrors.rds")
ref1 <- readRDS("lErrors.rds")

# example 1
txt <- capture.output(results <- by(data = trends, INDICES = trends[,"year"], FUN = function (y){
           dat <- reshape2::dcast(subset ( y, domain == "reading"), idstud~item, value.var="value")
           def <- defineModel(dat=dat, items= -1, id="idstud", software="tam")
           run <- runModel(def)
           res <- getResults(run)
           return(res)})) |> suppressMessages()
items   <- lapply(results, itemFromRes)
txt <- capture.output(eq.1_2  <- equat1pl(items[[1]][,c("item", "est")], items[[2]][,c("item", "est")],difBound = 0.64, iterativ = TRUE))
txt <- capture.output(eq.2_3  <- equat1pl(items[[2]][,c("item", "est")], items[[3]][,c("item", "est")],difBound = 0.64, iterativ = TRUE))
txt <- capture.output(eq.1_3  <- equat1pl(items[[1]][,c("item", "est")], items[[3]][,c("item", "est")],difBound = 0.64, iterativ = TRUE))


test_that("chained linking errors of example 1 are equal between package versions", {
  lErrors <- multiEquatError(eq.1_2, eq.2_3, eq.1_3, verbose=FALSE)
  expect_equal(lErrors[[1]][["le3221"]], ref1[[1]][["le3221"]])
})

# ref2 <- readRDS("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/lErrDif.rds")
ref2 <- readRDS("lErrDif.rds")

test_that("chained linking errors of example 1 (dependentDIF) are equal between package versions", {
  lErrDif <- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE, verbose=FALSE)
  expect_equal(lErrDif[[1]][["le3221"]], ref2[[1]][["le3221"]])
})

# ref3 <- readRDS("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/lErrors2.rds")
ref3 <- readRDS("lErrors2.rds")

# example 2
txt <- capture.output(results2<- by(data = trends, INDICES = trends[,"year"], FUN = function (y){
           qmat<- unique(y[ ,c("item","domain")])
           qmat<- data.frame ( qmat[,"item", drop=FALSE], model.matrix(~domain-1, data = qmat))
           spl <- splitModels ( qMatrix = qmat,  nCores = 1)
           dat <- reshape2::dcast(y, idstud~item, value.var="value")
           def <- defineModel(dat=dat, splittedModels = spl, id="idstud", software="tam")
           run <- runModel(def)
           res <- getResults(run)
           return(res)}))
txt <- capture.output(eq.1_2  <- equat1pl(results2[[1]], itemFromRes(results2[[2]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE))
txt <- capture.output(eq.2_3  <- equat1pl(results2[[2]], itemFromRes(results2[[3]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE))
txt <- capture.output(eq.1_3  <- equat1pl(results2[[1]], itemFromRes(results2[[3]])[,c("item", "est")],difBound = 0.64, iterativ = TRUE))

test_that("chained linking errors of example 1 (dependentDIF) are equal between package versions", {
  lErrors2<- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)
  lapply(1:2, FUN = function (i) { expect_equal(lErrors2[[i]][["le3221"]], ref3[[i]][["le3221"]])})
})

# ref4 <- readRDS("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/jkErr2.rds")
ref4 <- readRDS("jkErr2.rds")

# example 3
txt <- capture.output(eq.1_2  <- equat1pl(results2[[1]], itemFromRes(results2[[2]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE))
txt <- capture.output(eq.2_3  <- equat1pl(results2[[2]], itemFromRes(results2[[3]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE))
txt <- capture.output(eq.1_3  <- equat1pl(results2[[1]], itemFromRes(results2[[3]])[,c("item", "est")] |> dplyr::mutate(testlet = substr(item,1,3)), item = "item", testlet = "testlet", value = "est", difBound = 0.64, iterativ = TRUE))

test_that("chained linking errors of example 3 (dependentDIF) are approximately equal between package versions", {
  lErrors3<- multiEquatError(eq.1_2, eq.2_3, eq.1_3, dependentDIF =TRUE)
  lapply(1:2, FUN = function (i) { expect_equal(round(lErrors3[[i]][["le3221"]],digits = 2), round(ref4[[i]][["le3221"]], digits = 2))})
})


# example 4
# ref5 <- readRDS("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatModel/tests/testthat/ex4.rds")
ref5 <- readRDS("ex4.rds")

data(data.pars1.rasch, package="sirt")

# use the first three item parameter lists
e1 <- subset ( data.pars1.rasch, study == "study1")[,c("item", "b")]
e2 <- subset ( data.pars1.rasch, study == "study2")[,c("item", "b")]
e3 <- subset ( data.pars1.rasch, study == "study3")[,c("item", "b")]

# pairwise linking
txt <- capture.output(eq.1_2 <- equat1pl(e1, e2,difBound = 0.64, iterativ = TRUE))
txt <- capture.output(eq.2_3 <- equat1pl(e2, e3,difBound = 0.64, iterativ = TRUE))
txt <- capture.output(eq.1_3 <- equat1pl(e1, e3,difBound = 0.64, iterativ = TRUE))

test_that("chained linking errors of example 4 are equal between package versions", {
  lErrors4<- multiEquatError(eq.1_2, eq.2_3, eq.1_3)
  expect_equal(lErrors4[[1]][["le3221"]], ref5[["le3221"]])
})



