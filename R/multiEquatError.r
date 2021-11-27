# wird auf namespace exportiert, soll sowohl data.frames als auch results-objekte verarbeiten koennen
# inputobjekte heissen unspezifisch x1, x2, x3, weil es ja sowohl data.frames als auch results-objekte sein koennen
multiEquatError <- function(x1, x2, x3, difBound = 1, dependentDIF = FALSE, testletStr = NULL) {
  ### unspecific checks
  liste<- list(x1, x2, x3)
  chk1 <- checkInput(liste) #chk1 is defined but not used, brauch wir ggf. nicht als Rückgabe?
  if(!is.null(testletStr)) {
    # check input for equalting.rasch.jackknife: hier muss noch gecheckt werden, ob alle Items eine Unit haben und ausgegeben werden welche nicht, bzw. muss dann der Itemname als Unitname benutzt werden

    # dann muss einiges anders werden
    stopifnot(ncol(x1) == 3) # etc.
  }
  ### what kind of input?
  if ( "derived.par" %in% colnames(x1) ) {                                 ### specific checks for 'eatModelResults' object
    obj  <- prepareAndCheckEatModelObject(liste, difBound=difBound)     ### diese Funktion checkt das 'defModelObjList' Objekt und passt es so an, dass damit 'tripleEquatError' ausgefuehrt werden kann
    link <- lapply(names(obj), FUN = function (dim ) { tripleEquatError(e1=obj[[dim]][[1]][,c("item", "est")], e2=obj[[dim]][[2]][,c("item", "est")], e3= obj[[dim]][[3]][,c("item", "est")], dependentDIF=dependentDIF)})
    names(link) <- names(obj)
  } else {                                                                 ### specific checks for 'data.frame' object
    link <- tripleEquatError(e1=x1, e2=x2, e3=x3, dependentDIF=dependentDIF, testletStr=testletStr)
  }
  return(link) }


# hilfsfunktion (nicht auf NAMESPACE exportieren)
tripleEquatError <- function(e1, e2, e3, dependentDIF, testletStr) {

  if(is.null(testletStr)) {
      el12 <- sirt::equating.rasch(e1,e2)
      el13 <- sirt::equating.rasch(e1,e3)
      el23 <- sirt::equating.rasch(e2,e3)
      l12 <- el12$descriptives
      l13 <- el13$descriptives
      l23 <- el23$descriptives
      trend1223 <- -el12$B.est$Mean.Mean - el23$B.est$Mean.Mean
      trend13 <- -el13$B.est$Mean.Mean
     if(dependentDIF) {
       is1223 <- intersect(el12$anchor$item, el23$anchor$item)
       dif12 <- el12$anchor$TransfItempar.Gr1[match(is1223, el12$anchor$item)] - el12$anchor$Itempar.Gr2[match(is1223, el12$anchor$item)]
       dif23 <- el23$anchor$TransfItempar.Gr1[match(is1223, el23$anchor$item)] - el23$anchor$Itempar.Gr2[match(is1223, el23$anchor$item)]
        cov1223a <- cov(dif12, dif23)
        if(cov1223a < 0) cov1223a <- 0
        cov1223 <- sqrt(cov1223a)/sqrt(length(dif12))
        le1223 <- sqrt(l12$linkerror^2 + l23$linkerror^2 - 2*cov1223^2)
      } else {
        le1223 <- sqrt(l12$linkerror^2 + l23$linkerror^2)
      }
      le13 <- l13$linkerror
  } else {
      e12a <- merge(e1, e2, by = "item", all=TRUE)
      e12 <- merge(testletStr, e12a, by="item", all.x=FALSE, all.y=TRUE)
      e12 <- e12[,c(2,3,4,1)]

      e13a <- merge(e1, e3, by = "item", all=TRUE)
      e13 <- merge(testletStr, e13a, by="item", all.x=FALSE, all.y=TRUE)
      e13 <- e13[,c(2,3,4,1)]

      e23a <- merge(e2, e3, by = "item", all=TRUE)
      e23 <- merge(testletStr, e23a, by="item", all.x=FALSE, all.y=TRUE)
      e23 <- e23[,c(2,3,4,1)]

      eli12 <- sirt::equating.rasch.jackknife(e12)
      eli13 <- sirt::equating.rasch.jackknife(e13)
      eli23 <- sirt::equating.rasch.jackknife(e23)
      l12 <- eli12$descriptives$linkerror.jackknife
      l13 <- eli13$descriptives$linkerror.jackknife
      l23 <- eli23$descriptives$linkerror.jackknife
      if(dependentDIF) {
        n1 <- na.omit(e12)
        n1$est.x <- n1$est.x - mean(n1$est.x)
        n1$est.y <- n1$est.y - mean(n1$est.y)
        n2 <- na.omit(e23)
        n2$est.x <- n2$est.x - mean(n2$est.x)
        n2$est.y <- n2$est.y - mean(n2$est.y)
        e1223 <- na.omit(merge(n1, n2, by=c("unit","item")))
        dif12 <-(e1223[,4]-mean(e1223[,4])) - (e1223[,3]-mean(e1223[,3]))
        dif23 <- e1223[,6] - e1223[,5]
# alles falsch hier, haut überhaupt nicht hin, arbeite gleich morgen früh dran weiter ;)
        cov1223a <- tcrossprod(dif12, dif23)
        if(cov1223a < 0) cov1223a <- 0
        cov1223 <- sqrt(cov1223a)/sqrt(length(dif12))
        le1223 <- sqrt(l12$linkerror^2 + l23$linkerror^2 - 2*cov1223^2)
      } else {
        le1223 <- sqrt(l12$linkerror^2 + l23$linkerror^2)
      }
  }


  res <- data.frame(trend13 = trend13, trend1223 = trend1223, le13 = le13, le1223 = le1223)
  return(res)
}

checkInput <- function(inputlist) {
     cls <- lapply(inputlist, class)
     if(!all.equal(cls[[1]], cls[[2]], cls[[3]])) {stop("'x1', 'x2', and 'x3' must have the same class.")}
     if(!is.data.frame(inputlist[[1]])) {stop("'x1', 'x2', and 'x3' must be of class 'data.frame'.")} }

#multiEquatError <- function (defModelObjList, nCores=1 ){
#    ### checks
#       obj  <- checkInput(defModelObjList)                                      ### diese Funktion checkt das 'defModelObjList' Objekt und aendert es geringfuegig
#    ### linkingfehlerbestimmung separat fuer alle domaenen: erstmal alle Modelle mit ltm skalieren
#       les  <- lapply(names(obj[[1]]), FUN = function ( d ) {                   ### aeussere Schleife: Kompetenzbereiche
#               runs <- lapply(1:length(obj), FUN = function ( mzp ) {           ### innere Schleife: Zeitpunkte
#                       run <- ltm::rasch(obj[[mzp]][[d]][["daten"]][,-1], constraint = cbind(ncol(obj[[mzp]][[d]][["daten"]][,-1]) + 1, 1))
#                       txt <- capture.output(matr<- equateIRT::import.ltm(run)) ### ggf. option fuer multicore
#                       return(matr)})
#    ### Listenstruktur umdrehen
#               neu  <- list()
#               for ( i in 1:length(runs[[1]])) {
#                     for ( j in 1:length(runs)) {
#                           if ( j == 1) { neu[[i]] <- list()}
#                           neu[[i]][[j]] <- runs[[j]][[i]] }}
#    ### direktes Equating
#               nams <- paste("mzp", 1:length(neu[[1]]), sep = "")
#               mod  <- equateIRT::modIRT(coef = neu[[1]], var = neu[[2]], names = nams,	display = FALSE)
#               direc<- equateIRT::alldirec(mods = mod, method = "mean-mean", all = TRUE)
#    # indirektes Equating
#               chain<- equateIRT::chainec(direclist = direc, pths = nams)
#               summary(chain)
#               })
#}


#1 vs. 2
#r1 <- runModel(obj[[1]][[1]])
#r2 <- runModel(obj[[2]][[1]])

#res1 <- getResults(r1, omitPV=TRUE, omitFit = TRUE, omitWle=TRUE)
#res2 <- getResults(r2, omitPV=TRUE, omitFit = TRUE, omitWle=TRUE)

#prm1 <- itemFromRes(res1)
#eq <- equat1pl(res2,prm1[,c("item","est")])


# fuer ein einziges modell siehts so aus:
#Browse[1]> str(est.mod1pl) #  <- import.ltm(mod1pl)
#List of 2
# $ coef: num [1:5, 1:2] 2.73 0.999 0.24 1.306 2.099 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:5] "Item 1" "Item 2" "Item 3" "Item 4" ...
#  .. ..$ : chr [1:2] "beta.i" "beta"
# $ var : num [1:10, 1:10] 0.017015 0.001269 0.000747 0.001462 0.001882 ...

# fuer mehrere modelle soll es so aussehen
#Browse[1]> str(estrasch)
#List of 2
# $ coef:List of 5
#  ..$ : num [1:20, 1:2] 0.0118 -0.2343 0.7033 -0.0949 0.5893 ...
#  .. ..- attr(*, "dimnames")=List of 2
#  .. .. ..$ : chr [1:20] "I1" "I2" "I3" "I4" ...
#  .. .. ..$ : chr [1:2] "beta.i" "beta"
#  ..$ : num [1:20, 1:2] 0.1798 -0.0458 0.9153 0.0621 0.788 ...
#  .. ..- attr(*, "dimnames")=List of 2
#  .. .. ..$ : chr [1:20] "I1" "I2" "I3" "I4" ...
#  .. .. ..$ : chr [1:2] "beta.i" "beta"
#  ..$ : num [1:20, 1:2] -0.035 -0.4318 0.7366 0.0254 0.3023 ...
#  .. ..- attr(*, "dimnames")=List of 2
#  .. .. ..$ : chr [1:20] "I11" "I12" "I13" "I14" ...
#  .. .. ..$ : chr [1:2] "beta.i" "beta"
#  ..$ : num [1:20, 1:2] 0.83 -0.311 0.643 0.45 0.251 ...
#  .. ..- attr(*, "dimnames")=List of 2
#  .. .. ..$ : chr [1:20] "I21" "I22" "I23" "I24" ...
#  .. .. ..$ : chr [1:2] "beta.i" "beta"
#  ..$ : num [1:20, 1:2] -0.103 0.379 0.401 0.306 -0.283 ...
#  .. ..- attr(*, "dimnames")=List of 2
#  .. .. ..$ : chr [1:20] "I31" "I32" "I33" "I34" ...
#  .. .. ..$ : chr [1:2] "beta.i" "beta"
# $ var :List of 5
#  ..$ : num [1:20, 1:20] 0.001169 0.000202 0.000202 0.000202 0.000202 ...
#  ..$ : num [1:20, 1:20] 0.001174 0.000202 0.000202 0.000202 0.000202 ...
#  ..$ : num [1:20, 1:20] 0.00117 0.000202 0.000202 0.000202 0.000202 ...
#  ..$ : num [1:20, 1:20] 0.001278 0.000201 0.000204 0.000203 0.000203 ...
#  ..$ : num [1:20, 1:20] 0.001172 0.000202 0.000202 0.000202 0.000202 ...

prepareAndCheckEatModelObject <- function ( liste, difBound ) {
    ### all models unidim?
       uni  <- sapply(liste, FUN = function ( x ) { "correlation" %in% x[,"par"] })
       chk1 <- which(uni == TRUE)
       if ( length(chk1)>0) {stop("Following model(s) do not seem to be unidimensional: ",paste(chk1, collapse = ", "),".")}
    ### all models equal domains?
       doms <- lapply(liste, FUN = function (o) {na.omit(setdiff(unique(o[,"group"]), "persons"))})
       comps<- data.frame(combinat::combn(1:length(doms),2))
       chk2 <- sapply(comps, FUN = function ( col ) { if(!all.equal(doms[[col[1]]], doms[[col[2]]])) {stop(paste0("Domains between ",col[1]," and ",col[2]," do not match."))}})
    ### common items between adjacent measurement time points?
       its  <- lapply(liste, itemFromRes)
       dims <- unique(its[[1]][,"dimension"])
       its  <- lapply(dims, FUN = function (d) {  list(t1 = its[[1]][which(its[[1]][,"dimension"] == d),], t2 = its[[2]][which(its[[2]][,"dimension"] == d),], t3 = its[[3]][which(its[[3]][,"dimension"] == d),]) })
       names(its) <- dims
       chk3 <- lapply(dims, FUN = function ( d ) {
               bridges <- data.frame(combinat::combn(1:length(its[[d]]),2), stringsAsFactors = FALSE)
               bridges <- bridges[,which(sapply(bridges, diff)==1)]
               chk4    <- lapply(bridges, FUN = function (b ) {
                          com<- length(intersect( its[[d]][[b[1]]][,"item"], its[[d]][[b[2]]][,"item"]))
                          if ( com == 0 ) {stop(paste0("No common items for linking mzp ",b[1]," to mzp ", b[2]," for domain '",d,"'."))}
                          if ( com < 10 ) {cat(paste0("Warning: Only ",com," common items for linking mzp ",b[1]," to mzp ", b[2]," for domain '",d,"'.\n"))}
                          })})
    ### link adjacent measurement time points to exclude linking dif items
       ref <- itemFromRes(liste[[1]])
    ### 2 vs. 1
       eq  <- equat1pl(liste[[2]], prmNorm = ref[,c("item", "est")], difBound = difBound, iterativ = TRUE)
       stru<- lapply(eq[["items"]], names)
       for ( d in 1:length(stru)) {
           weg <- setdiff(eq[["items"]][[names(stru)[d]]][[stru[[d]]]][["info"]][["itemExcluded"]], "")
           if ( length(weg)>0) {
                its[[stru[[d]]]][[1]] <- its[[stru[[d]]]][[1]][-match(weg, its[[stru[[d]]]][[1]][,"item"]),]
           }
       }
    ### 3 vs. 2
       ref2<- itemFromRes(liste[[2]])
       eq  <- equat1pl(liste[[3]], prmNorm = ref2[,c("item", "est")], difBound = difBound, iterativ = TRUE)
       stru<- lapply(eq[["items"]], names)
       for ( d in 1:length(stru)) {
           weg <- setdiff(eq[["items"]][[names(stru)[d]]][[stru[[d]]]][["info"]][["itemExcluded"]], "")
           if ( length(weg)>0) {
                its[[d]][[3]] <- its[[d]][[3]][-match(weg, its[[d]][[3]][,"item"]),]
           }
       }
       return(its)}

