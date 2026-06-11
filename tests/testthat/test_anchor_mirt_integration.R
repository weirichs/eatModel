# Integrationstest nach dem Beispiel in tests/code_review.R: GPCM-Kalibrierung in
# TAM, vier Items werden in mirt verankert (Anker in traditioneller b-Parametri-
# sierung), danach wird das gefittete mirt-Modell geprueft. Langsam (zwei Fits)!

test_that("Integration: TAM-Ankerwerte werden im verankerten mirt-Modell reproduziert", {
    skip_on_cran()
    data(reading)
    dw   <- reshape2::dcast(subset(reading, bookletID == "TH08"), idstud+sex~item, value.var = "valueSum") |>
            dplyr::mutate(sex = car::recode(sex, "'male'=0; 'female'=1", as.factor=FALSE)) |>
            eatTools::na_omit_selection("sex")
    defT <- defineModel(dat = dw, items = -c(1:2), id = "idstud", irtmodel = "GPCM", software = "tam")
    runT <- runModel(defT)
    it   <- itemFromRes(getResults(runT))
    ankItems <- c("D025033", "D204013", "D205143", "D225013")
    anchor   <- subset(it, item %in% ankItems)[, c("item", "category", "est")]
    slope    <- unique(subset(it, item %in% c(ankItems, "D225053"))[, c("item", "estSlope")])
    pcItems  <- sapply(colnames(dw)[-c(1:2)], FUN = function(i) {any(dw[, i] %in% 2)})
    items    <- data.frame(item = names(pcItems), irtmod = ifelse(pcItems, "gpcm", "2PL"), stringsAsFactors = FALSE)
    defM <- defineModel(dat = dw, items = -c(1:2), id = "idstud", irtmodel = items, software = "mirt",
                        anchor = anchor, itemCol = "item", valueCol = "est", catCol = "category", fixSlopeMat = slope)
    runM <- runModel(defM)
    it_mirt <- itemFromRes(getResults(runM))
    ## 1) interne mirt-Parametrisierung (Intercept/Slope): die fixierten Parameter
    ##    des polytomen Testitems muessen den via traditional2mirt umgerechneten
    ##    Ankerwerten entsprechen, d.h. a1 = a und d_k = -a*(b_1 + ... + b_k)
    item_test<- "D204013"
    sub  <- anchor[which(anchor$item == item_test), ]
    b    <- sub[order(sub$category), "est"]
    a    <- slope[which(slope$item == item_test), "estSlope"]
    expected <- mirt:::traditional2mirt(c(a, b), "gpcm", ncat = length(b) + 1)
    cf   <- coef(runM, IRTpars = FALSE)[[item_test]]["par", ]
    expect_equal(cf[["a1"]], expected[["a1"]])
    expect_equal(cf[["d1"]], expected[["d1"]])
    expect_equal(cf[["d2"]], expected[["d2"]])
    ## 2) Roundtrip in traditioneller Parametrisierung: die verankerten Schwellen
    ##    und Trennschaerfen muessen unveraendert in itemFromRes(resM) auftauchen
    mrg  <- merge(anchor, it_mirt[, c("item", "category", "est")], by = c("item", "category"), suffixes = c("_ank", "_mirt"))
    expect_equal(nrow(mrg), nrow(anchor))
    expect_equal(mrg[, "est_mirt"], mrg[, "est_ank"], tolerance = 1e-4)
    mrgS <- merge(slope, unique(it_mirt[, c("item", "estSlope")]), by = "item", suffixes = c("_ank", "_mirt"))
    expect_equal(mrgS[, "estSlope_mirt"], mrgS[, "estSlope_ank"], tolerance = 1e-4)
})
