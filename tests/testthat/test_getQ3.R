q3MinReference <- function(dat, q3MinType) {
  itemPairs <- utils::combn(seq_len(ncol(dat)), 2)
  minValue <- apply(itemPairs, MARGIN = 2, FUN = function(pair) {
    tab <- table(dat[, pair, drop = FALSE])
    if(any(dim(tab) < 2)) {return(0)}
    if(q3MinType == "singleObs") {
      min(tab)
    } else {
      min(c(rowSums(tab), colSums(tab)))
    }
  })
  data.frame(Var1 = colnames(dat)[itemPairs[1,]], Var2 = colnames(dat)[itemPairs[2,]],
             minValue = as.numeric(minValue), stringsAsFactors = FALSE)
}

test_that("nObsItemPairs handles item pairs without pairwise complete observations", {
  dat <- data.frame(v1 = c(NA, NA, 1, 2),
                    v2 = c(1, 1, NA, NA),
                    stringsAsFactors = FALSE)

  expect_equal(table(dat), structure(c(0L, 0L), dim = c(2L, 1L),
                                     dimnames = list(v1 = c("1", "2"), v2 = "1"),
                                     class = "table"))
  expect_equal(nObsItemPairs(dat, q3MinType = "singleObs"),
               data.frame(Var1 = "v1", Var2 = "v2", minValue = 0, stringsAsFactors = FALSE))
  expect_equal(nObsItemPairs(dat, q3MinType = "marginalSum"),
               data.frame(Var1 = "v1", Var2 = "v2", minValue = 0, stringsAsFactors = FALSE))
})

test_that("nObsItemPairs matches table-based reference with missing values", {
  set.seed(20260817)
  dat <- data.frame(v1 = sample(c(0:2, NA), 80, replace = TRUE),
                    v2 = sample(c(0:1, NA), 80, replace = TRUE),
                    v3 = sample(c(1:3, NA), 80, replace = TRUE),
                    v4 = sample(c(NA, NA, 1), 80, replace = TRUE),
                    stringsAsFactors = FALSE)

  for(q3MinType in c("singleObs", "marginalSum")) {
    expect_equal(nObsItemPairs(dat, q3MinType = q3MinType),
                 q3MinReference(dat, q3MinType = q3MinType))
  }
})

test_that("nObsItemPairs respects factor levels like table", {
  dat <- data.frame(v1 = factor(c(1, 1, NA), levels = c(1, 2)),
                    v2 = factor(c(1, 2, 2), levels = c(1, 2)))

  expect_equal(nObsItemPairs(dat, q3MinType = "singleObs"),
               q3MinReference(dat, q3MinType = "singleObs"))
  expect_equal(nObsItemPairs(dat, q3MinType = "marginalSum"),
               q3MinReference(dat, q3MinType = "marginalSum"))
})
