makeMiniTrafo <- function() {
  domainEntry <- function(domain, item, parameter, extra = NULL) {
    anchor <- data.frame(item = item, parameter = parameter, stringsAsFactors = FALSE)
    if(!is.null(extra)) {anchor[["extra"]] <- extra}
    list(
      refPop = data.frame(domain = domain, m = parameter, sd = 1, stringsAsFactors = FALSE),
      cuts = stats::setNames(list(list(values = c(400, 500), labels = c("I", "II", "III"))), domain),
      anchor = anchor,
      info = stats::setNames(paste("info", domain), domain)
    )
  }

  list(
    paper = list(
      primary = list(
        mat = list(vera = list(
          GL = domainEntry("GL", "m1", 0),
          MS = domainEntry("MS", "m2", 1)
        )),
        deu = list(vera = list(
          lesen = domainEntry("lesen", "d1", 2, extra = "kept")
        ))
      )
    )
  )
}

test_that("getTrafo combines transformation entries across subjects and domains", {
  out <- getTrafo(
    dataBase = makeMiniTrafo(), mode = "paper", grade = "primary",
    subject = c("math", "deu"), domain = "all", study = "vera"
  )

  expect_setequal(out$refPop$domain, c("GL", "MS", "lesen"))
  expect_setequal(names(out$cuts), c("GL", "MS", "lesen"))
  expect_equal(length(out$info), 3)
  expect_true(all(c("item", "parameter", "extra") %in% colnames(out$anchor)))
  expect_true(is.na(out$anchor[match("m1", out$anchor$item), "extra"]))
  expect_equal(out$anchor[match("d1", out$anchor$item), "extra"], "kept")
})

test_that("getTrafo skips subjects without the requested domain", {
  out <- getTrafo(
    dataBase = makeMiniTrafo(), mode = "paper", grade = "primary",
    subject = c("math", "deu"), domain = "GL", study = "vera"
  )

  expect_equal(out$refPop$domain, "GL")
  expect_equal(names(out$cuts), "GL")
  expect_equal(out$anchor$item, "m1")
})
