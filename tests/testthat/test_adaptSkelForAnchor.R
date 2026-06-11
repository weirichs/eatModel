# minimaler Unit-Test fuer adaptSkelForAnchor(): die Funktion wird direkt aufgerufen,
# der skeleton kommt aus mirt(..., pars = "values"), d.h. es wird nichts geschaetzt.
# Verankerung erfolgt in traditioneller IRT-Parametrisierung (b), mirt erwartet
# Intercepts (d). Transformation dichotom: d = -a*b (vgl. mirt:::traditional2mirt)

buildSkelInput <- function() {
    set.seed(42)
    dat <- data.frame(I1 = rbinom(200, 1, 0.5), I2 = rbinom(200, 1, 0.5),
                      I3 = rbinom(200, 1, 0.5), I4 = sample(0:2, 200, replace = TRUE))
    irtmodel <- data.frame(item = paste0("I", 1:4), type = c("2PL", "2PL", "2PL", "gpcm"), stringsAsFactors = FALSE)
    skel <- mirt::mirt(dat, 1, itemtype = irtmodel[, "type"], pars = "values", verbose = FALSE)
    list(skel = skel, irtmodel = irtmodel,
         allNam = list(slopeMatItemCol = "item", slopeMatValueCol = "estSlope"),
         qmat = data.frame(item = paste0("I", 1:4), dim1 = 1, stringsAsFactors = FALSE))
}

test_that("anchored dichotomous item: slope and intercept fixed, d = -a*b", {
    ## Itempars in mirt can be called as irt parameters via  coef(..., IRTpars = TRUE))
    ## However, the intercept and slope are needed for anchoring.  
    
    ## pretend I have called the irt parameters before. This are the parameters for the anchor item:
    b = 0.5
    a = 1.2
    inp  <- buildSkelInput()
    anch <- list(ank = data.frame(item = "I1", parameter = b, stringsAsFactors = FALSE))
    slop <- list(ori = data.frame(item = "I1", estSlope = a, stringsAsFactors = FALSE))
    out  <- adaptSkelForAnchor(allNam = inp$allNam, skel = inp$skel, anch = anch,
                               qmat = inp$qmat, slope = slop, irtmodel = inp$irtmodel,
                               est.slopegroups = NULL)
    rowA <- which(out$item == "I1" & out$name == "a1")
    rowD <- which(out$item == "I1" & out$name == "d")
    expect_false(out[rowA, "est"])
    ## a stays as it is
    expect_equal(out[rowA, "value"], a)
    expect_false(out[rowD, "est"])
    ## We should get the transformed parameter d: 
    expect_equal(out[rowD, "value"], -1*(a * b)) ## formula: -1(b * a) = d
    # Mittelwert der Dimension mit Ankeritem sowie Varianz (Dimension enthaelt
    # Item mit fixierter Trennschaerfe) muessen frei geschaetzt werden
    expect_true(out[which(out$name == "MEAN_1"), "est"])
    expect_true(out[which(out$name == "COV_11"), "est"])
    # nicht verankerte Items bleiben frei
    expect_true(out[which(out$item == "I2" & out$name == "d"), "est"])
})


## also check against coef(itempars = TRUE )



test_that("anchored gpcm item: skeleton values match mirt:::traditional2mirt", {
    # Achtung: kodiert die mirt-Konvention (mirt:::traditional2mirt, itemclass 'gpcm');
    # schlaegt derzeit fehl, weil adaptSkelForAnchor d_k = -a*b_k ohne cumsum setzt
    a    <- 0.8                                  # Trennschaerfe des Ankeritems
    b    <- c(-0.5, 1.0)                         # traditionelle Schwellen b1, b2 (aus coef(..., IRTpars = TRUE))
    inp  <- buildSkelInput()
    anch <- list(ank = data.frame(item = c("I4", "I4"), parameter = b,
                                  category = c("Cat1", "Cat2"), stringsAsFactors = FALSE))
    slop <- list(ori = data.frame(item = "I4", estSlope = a, stringsAsFactors = FALSE))
    out  <- adaptSkelForAnchor(allNam = inp$allNam, skel = inp$skel, anch = anch,
                               qmat = inp$qmat, slope = slop, irtmodel = inp$irtmodel,
                               est.slopegroups = NULL)
    # traditional2mirt erwartet x = c(a, b1, ..., bK): Trennschaerfe zuerst, dann
    # die Schwellen; die Reihenfolge zaehlt, etwaige Namen des Vektors werden ignoriert.
    # ncat = Anzahl der Antwortkategorien (hier 3: 0/1/2), NICHT Anzahl der Schwellen.
    # Rueckgabe ist benannt: a1, ak0..ak2, d0..d2 (d0 ist immer 0)
    mirt_pars <- mirt:::traditional2mirt(c(a, b), "gpcm", ncat = length(b) + 1)
    rowA  <- which(out$item == "I4" & out$name == "a1")
    rowD1 <- which(out$item == "I4" & out$name == "d1")
    rowD2 <- which(out$item == "I4" & out$name == "d2")
    # verankert werden die Trennschaerfe (a1) und alle Intercepts ausser d0
    expect_false(out[rowA,  "est"])
    expect_false(out[rowD1, "est"])
    expect_false(out[rowD2, "est"])
    expect_equal(out[rowA,  "value"], mirt_pars[["a1"]])    # = a
    expect_equal(out[rowD1, "value"], mirt_pars[["d1"]])    # = -a*b1
    expect_equal(out[rowD2, "value"], mirt_pars[["d2"]])    # = -a*(b1+b2)
})



