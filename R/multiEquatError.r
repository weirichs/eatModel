# wird auf namespace exportiert, soll sowohl data.frames als auch results-objekte verarbeiten koennen
# inputobjekte heissen unspezifisch x1, x2, x3, weil es ja sowohl data.frames als auch results-objekte sein koennen
multiEquatError <- function (x1, x2, x3, difBound = 1, dependentDIF =FALSE, testletStr = NULL, verbose=TRUE ){
### unspecific checks
  liste<- lapply(list(x1, x2, x3), eatTools::makeDataFrame, verbose=FALSE) ### specific checks for 'eatModelResults' object
### what kind of input?                                                     ### diese Funktion checkt das 'defModelObjList' Objekt und passt es so an, dass damit 'tripleEquatError' ausgefuehrt werden kann
  if("derived.par" %in% colnames(x1) ) {                                 ### wenn der Input das Rueckgabeobjekt von eatModel ist, wird es hier in drei data.frames mit je zwei Spalten umgewandelt. Das ist der Input, den 'tripleEquatError' braucht
    obj  <- prepareAndCheckEatModelObject(liste, difBound=difBound, verbose=verbose)
    link <- lapply(names(obj), FUN = function (dim ) {
            tripleEquatError(e1=obj[[dim]][[1]][,c("item", "est")], e2=obj[[dim]][[2]][,c("item", "est")],
               e3= obj[[dim]][[3]][,c("item", "est")], dependentDIF=dependentDIF, testletStr=testletStr,
               difBound=difBound, dimName = dim, verbose=verbose)})
    names(link) <- names(obj)
  } else {                                                                 ### specific checks for 'data.frame' object
    link <- tripleEquatError(e1=liste[[1]], e2=liste[[2]], e3=liste[[3]],
            dependentDIF=dependentDIF, testletStr=testletStr, difBound=difBound,
            dimName = "Dim1", verbose=verbose)
  }
  return(link) }


# hilfsfunktion (nicht auf NAMESPACE exportieren)
tripleEquatError <- function(e1, e2, e3, dependentDIF, testletStr, difBound, dimName, verbose) {
### checks
  chk1 <- checkInputConsistency(e1=e1,e2=e2,e3=e3,testletStr=testletStr)    ### Die Rueckgabe von 'checkInputConsistency' ist NUR die Testletstruktur.
  liste<- list(e1, e2, e3)                                                  ### Aber Achtung! Noch nicht fuer Testletstruktur! Achtung! bei 'checkIPD()' wird 'chk1' reingeschrieben statt 'testletStr', da das Objekt ggf. von der check-Funktion veraendert wird
  kombs<- combinat::combn(1:3, 2, simplify=FALSE)
  eqs  <- lapply(kombs, FUN = function (k) {
  txt  <- capture.output(eqx <- equat1pl(liste[[k[1]]], liste[[k[2]]],  difBound = difBound, iterativ = TRUE))
    return(eqx[["items"]][["global"]][["global"]][["eq"]][["anchor"]][,c(1,2,4)])})
    names(eqs) <- unlist(lapply(kombs, FUN = paste, collapse="_"))
  el21 <- equ.rasch(eqs[["1_2"]][,c(1,3)], eqs[["1_2"]][,c(1,2)])
  el31 <- equ.rasch(eqs[["1_3"]][,c(1,3)], eqs[["1_2"]][,c(1,2)])
  el32 <- equ.rasch(eqs[["2_3"]][,c(1,3)], eqs[["1_2"]][,c(1,2)])
  l21 <- el21$descriptives
  l31 <- el31$descriptives
  l32 <- el32$descriptives
  trend3221 <- el32$B.est$Mean.Mean + el21$B.est$Mean.Mean
  trend31 <- el31$B.est$Mean.Mean
### nicht nur wenn testletStr NULL ist, soll auf jackknife verzichtet werden, sondern auch wenn testlets missings haben oder fehlen
  if(is.null(chk1)) {                                                         ### hier wird 'chk1' reingeschrieben statt 'testletStr', weil die oben aufgerufene
    if(dependentDIF) {                                                         ### 'checkInputConsistency'-Funktion das Testletobjekt ggf. aendert oder auf NULL
      is3221 <- intersect(el21$anchor$item, el32$anchor$item)                  ### setzt, wenn es bspw. missings auf der Testletvariablen gibt
      dif21 <- el21$anchor$TransfItempar.Gr1[match(is3221, el21$anchor$item)] - el21$anchor$Itempar.Gr2[match(is3221, el21$anchor$item)]
      dif32 <- el32$anchor$TransfItempar.Gr1[match(is3221, el32$anchor$item)] - el32$anchor$Itempar.Gr2[match(is3221, el32$anchor$item)]
      cov1223a <- cov(dif21, dif32)
      if(cov1223a < 0) {cov1223a <- 0}
      cov1223 <- sqrt(cov1223a)/sqrt(length(dif21))
      le3221 <- sqrt(l21$linkerror^2 + l32$linkerror^2 - 2*cov1223^2)
    } else {
      le3221 <- sqrt(l21$linkerror^2 + l32$linkerror^2)
    }
    le31 <- l31$linkerror
  } else {
    e12a <- merge(e1, e2, by = "item", all=TRUE)
### Achtung! hier wird 'chk1' reingeschrieben statt 'testletStr', da das Objekt ggf. von der check-Funktion veraendert wird ('item'-Spalte wird
### in 'item' umbenannt, falls sie nicht schon so heisst), und das veraenderte Objekt hier benutzt werden soll!
    e12 <- merge(chk1, e12a, by="item", all.x=FALSE, all.y=TRUE)
    e12 <- e12[,c(2,3,4,1)]
    e13a <- merge(e1, e3, by = "item", all=TRUE)
    e13 <- merge(chk1, e13a, by="item", all.x=FALSE, all.y=TRUE)
    e13 <- e13[,c(2,3,4,1)]
    e23a <- merge(e2, e3, by = "item", all=TRUE)
    e23 <- merge(chk1, e23a, by="item", all.x=FALSE, all.y=TRUE)
    e23 <- e23[,c(2,3,4,1)]
    eli12 <- equa.rasch.jk(e12)
    eli13 <- equa.rasch.jk(e13)
    eli23 <- equa.rasch.jk(e23)
    l12 <- eli12$descriptives$linkerror.jackknife
    le31 <- eli13$descriptives$linkerror.jackknife
    l23 <- eli23$descriptives$linkerror.jackknife
    if(dependentDIF) {
      eli12$pars.data[,2] <- eli12$pars.data[,2]+eli12$descriptives$shift
      eli23$pars.data[,2] <- eli23$pars.data[,2]+eli23$descriptives$shift
      e1223 <- stats::na.omit(merge(eli12$pars.data, eli23$pars.data, by=c("unit","item"),all=TRUE))
      e1223$dif12 <- e1223[,3] - e1223[,4]
      e1223$dif23 <- e1223[,5] - e1223[,6]
      res1 <- NULL
      for(un in e1223$unit) {
        p1223 <- e1223[e1223[,1] != un,]
        res1 <- c(cov(p1223$dif12,p1223$dif23),res1)
        }
      cov1223a <- mean(res1, na.rm=TRUE)
      if(cov1223a < 0) {cov1223a <- 0}
      cov1223 <- sqrt(cov1223a)/sqrt(length(p1223$dif12))
      if((l12^2 + l23^2 - 2*cov1223^2) < 0) {
        le3221 <- sqrt(l12^2 + l23^2)
      } else {
        le3221 <- sqrt(l12^2 + l23^2 - 2*cov1223^2)
       }
    } else {
      le3221 <- sqrt(l12^2 + l23^2)
    }
  }                                                                             ### untere Zeile: Anzeige der direkten Linkingfehler
  message(paste0("\nDimension '",dimName,"': Direct linking error of the three combinations of measurements: \n",eatTools::print_and_capture(eatTools::roundDF(plyr::rbind.fill(data.frame ( mzp1 = 1, mzp2 = 2, l21), data.frame ( mzp1 = 1, mzp2 = 3, l31, chained = le3221, mzp_1.vs.3= trend31, mzp_1.vs.2.vs.3 = trend3221), data.frame ( mzp1 = 2, mzp2 = 3, l32)), digits = 3), spaces = 5)))
  res <- data.frame(trend31 = trend31, trend3221 = trend3221, le31 = le31, le3221 = le3221)
  return(res)}

prepareAndCheckEatModelObject <- function ( liste, difBound, verbose ) {
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
       its  <- lapply(dims, FUN = function (d) {  list(t1 = its[[1]][which(its[[1]][,"dimension"] == d),c("item", "est")], t2 = its[[2]][which(its[[2]][,"dimension"] == d),c("item", "est")], t3 = its[[3]][which(its[[3]][,"dimension"] == d),c("item", "est")]) })
       names(its) <- dims
       chk3 <- lapply(dims, FUN = function ( d ) {
               bridges <- data.frame(combinat::combn(1:length(its[[d]]),2), stringsAsFactors = FALSE)
               bridges <- bridges[,which(sapply(bridges, diff)==1)]
               chk4    <- lapply(bridges, FUN = function (b ) {
                          com<- length(intersect( its[[d]][[b[1]]][,"item"], its[[d]][[b[2]]][,"item"]))
                          if ( com == 0 ) {stop(paste0("No common items for linking mzp ",b[1]," to mzp ", b[2]," for domain '",d,"'."))}
                          if ( com < 10 ) {cat(paste0("Warning: Only ",com," common items for linking mzp ",b[1]," to mzp ", b[2]," for domain '",d,"'.\n"))}
                          })})
       return(its)}
       
### checkfunktion
checkInputConsistency <- function(e1,e2,e3,testletStr) {
    ### alle data.frames zwei Spalten?
       il  <- list(e1, e2, e3)
       if(!all(sapply(il, ncol) == 2)) {stop("All three data.frames must have two columns.")}
    ### erste Spalte muss 'item' heissen, zweite Spalte egal, aber heisst einfach mal 'est'
       if(!all(lapply(il, FUN = function ( x ) { all(colnames(x) == c("item", "est")) })==TRUE)) {
          il <- lapply(il, FUN = function ( x ) { colnames(x) <- c("item", "est"); return(x)})
       }
    ### unique item identifiers and common item pars between measurement occasions?
       it  <- lapply(il, FUN = function (x) {
                     stopifnot ( length(x[,"item"]) == length(unique(x[,"item"])))
                     return(x[,"item"])})
       comp<- as.data.frame(combinat::combn(length(it),2))
       chk1<- lapply(comp, FUN = function ( com ){ if ( length(intersect(it[[com[1]]], it[[com[2]]])) < 2) { stop(paste0(length(intersect(it[[com[1]]], it[[com[2]]])) ," common items between measurement occasions ",com[1], " and ", com[2], "."))} })
    ### testlet checken
       if (!is.null(testletStr)) {
            testletStr <- eatTools::makeDataFrame(testletStr, verbose=FALSE)
            if ( ncol(testletStr) != 2) {stop("'testletStr' must have two columns.")}
            if(!"item" %in% colnames(testletStr)) { stop("One column in 'testletStr' must be named 'item'.")}
            if (setdiff(colnames(testletStr), "item") != "unit") {
                cat(paste0("The 'unit' column in 'testletStr' seemed to be called '",setdiff(colnames(testletStr), "item"),"'. Rename this column into 'unit'.\n"))
                colnames(testletStr)[setdiff(1:2, match("item", colnames(testletStr)))] <- "unit"
            }
    ### pruefen, ob alle Item, die es in e1, e2, e3 gibt, auch eine NICHT-FEHLENDE Testletzuordnung haben!
            allI<- unique(unlist(lapply(il, FUN = function (x){x[,"item"]})))
            if(!all(allI %in% testletStr[,"item"])) {
                message(paste0("Following items without testlet definition: '",paste(setdiff(allI,testletStr[,"item"]) , collapse="', '"), "'.\n    Linking errors will be compute without considering testlets."))
                testletStr <- NULL
            }  else  {
                ntl <- length(unique(testletStr[,"unit"]))
                if(ntl < 2) {stop(paste0("Number of testlets (",ntl,") is too low."))}
                if(ntl < 5) {warning(paste0("Only ",ntl," testlets defined."))}
            }
            tls <- testletStr[which(testletStr[,"item"] %in% allI),]
            isna<- length(which(is.na(tls[,setdiff(colnames(tls), "item")])))
            if(isna>0) {
                message(paste0("Following ",isna," items without testlet definition: '",paste(tls[which(is.na(tls[,setdiff(colnames(tls), "item")])),"item"] , collapse="', '"), "'.\n    Linking errors will be compute without considering testlets."))
                testletStr <- NULL
            }
    ### wenn es testlets gibt, soll das testlet-Objekt zurueckgegeben werden, entweder unveraendert oder veraendert (Spaltenbezeichnung angepasst)
            return(testletStr)
       }  else  {
            return(NULL)
       } }

equa.rasch.jk <- function (pars.data, se.linkerror = FALSE, alpha1 = 0, alpha2 = 0) {
  pars.data <- as.data.frame(stats::na.omit(pars.data))
  itemunits <- unique(pars.data[, 1])
  N.units <- length(itemunits)
  N.items <- nrow(pars.data)
  mod1 <- equ.rasch(pars.data[, c(4, 2)], pars.data[, c(4, 3)])
  res1 <- data.frame(unit = itemunits, shift = 0, SD = 0, linkerror = 0)
  for (nn in 1:N.units) {
    pars.data1 <- pars.data[pars.data[, 1] != itemunits[nn], ]
    mod.nn <- equ.rasch(x = pars.data1[, c(4, 2)], y = pars.data1[, c(4, 3)])
    res1[nn, "shift"] <- mod.nn$B.est$Mean.Mean
    res1[nn, "SD"] <- mod.nn$descriptives$SD
    if (se.linkerror) {
      itemunits.nn <- itemunits[-nn]
      l1 <- NULL
      for (ii in itemunits.nn) {
        pars.data1.ii <- pars.data1[paste(pars.data1[, 1]) != ii, ]
        mod.ii <- equ.rasch(x = pars.data1.ii[,  c(4, 2)], y = pars.data1.ii[, c(4, 3)], alpha1 = alpha1, alpha2 = alpha2)
        l1 <- c(l1, mod.ii$B.est$Mean.Mean)
      }
      res1[nn, "linkerror"] <- sqrt((N.units - 2)/(N.units -  1) * sum((l1 - res1[nn, "shift"])^2))
    }
  }
  linkerror <- sqrt((N.units - 1)/N.units * sum((res1[, 2] - mod1$B.est$Mean.Mean)^2))
  se.sd <- sqrt((N.units - 1)/N.units * sum((res1[, 3] - mod1$descriptives$SD)^2))
  if (se.linkerror) {
    se.linkerror <- sqrt((N.units - 1)/N.units * sum((res1[, 4] - linkerror)^2))
  }
  else {
    se.linkerror <- NA
  }
  descriptives <- data.frame(N.items = N.items, N.units = N.units,
                             shift = mod1$B.est$Mean.Mean, SD = mod1$descriptives$SD,
                             linkerror.jackknife = linkerror, SE.SD.jackknife = se.sd,
                             se.linkerror.jackknife = se.linkerror)
  res <- list(pars.data = pars.data, itemunits = itemunits, descriptives = descriptives)
  return(res)
}

equ.rasch <- function (x, y, theta = seq(-4, 4, len = 100), alpha1 = 0, alpha2 = 0) {
  x[, 1] <- gsub(" ", "", paste(x[, 1]))
  y[, 1] <- gsub(" ", "", paste(y[, 1]))
  b.xy <- stats::na.omit(merge(x, y, by = "item"))
  colnames(b.xy) <- c("item", "Itempar.Gr1", "Itempar.Gr2")
  B.mm <- mean(b.xy[, 3]) - mean(b.xy[, 2])
  g1   <- prob_raschtype_genlogis(theta = theta, b = b.xy[, 2], alpha1 = 0, alpha2 = 0)
  opt_interval <- 10 * c(-1, 1)
  ha <- function(B) {
    fct1 <- prob_raschtype_genlogis(theta = theta, b = b.xy[, 2], alpha1 = alpha1, alpha2 = alpha2)
    fct2 <- prob_raschtype_genlogis(theta = theta, b = b.xy[, 3] - B, alpha1 = alpha1, alpha2 = alpha2)
    sum((fct1 - fct2)^2)
  }
  B.ha <- stats::optimize(f = ha, interval = opt_interval)$minimum
  sl <- function(B) {
    fct1 <- prob_raschtype_genlogis(theta = theta, b = b.xy[, 2], alpha1 = alpha1, alpha2 = alpha2)
    fct2 <- prob_raschtype_genlogis(theta = theta, b = b.xy[, 3] - B, alpha1 = alpha1, alpha2 = alpha2)
    sum((rowSums(fct1 - fct2))^2)
  }
  B.sl <- stats::optimize(f = sl, interval = opt_interval)$minimum
  B.est <- data.frame(B.mm, B.ha, B.sl)
  colnames(B.est) <- c("Mean.Mean", "Haebara","Stocking.Lord")
  b.xy$TransfItempar.Gr1 <- b.xy[, 2] + B.est[1, "Mean.Mean"]
  x[, 2] <- x[, 2] + B.est[1, "Mean.Mean"]
  transf.par <- merge(x = x, y = y, by.x = 1, by.y = 1, all = TRUE)
  colnames(transf.par) <- c("item", "TransfItempar.Gr1", "Itempar.Gr2")
  transf.par <- transf.par[order(paste(transf.par$item)), ]
  des <- data.frame(N.Items = nrow(b.xy), SD = stats::sd(b.xy$TransfItempar.Gr1 -  b.xy$Itempar.Gr2), stringsAsFactors = FALSE)
  des$Var <- des$SD^2
  des$linkerror <- sqrt(des[,"SD"]^2/des[,"N.Items"])
  res <- list(B.est = B.est, descriptives = des, anchor = b.xy[,  c(1, 2, 4, 3)], transf.par = transf.par)
  return(res)}

prob_raschtype_genlogis <- function( theta, b, alpha1, alpha2, fixed.a=1+0*b,  Qmatrix=NULL, dimensions=NULL ) {
  LT <- length(theta)
  if (is.matrix(theta)){
    LT <- nrow(theta)
  }

  if ( is.null(Qmatrix) ){
    XX <- TAM::tam_outer(x=theta, y=b, op="-")
    XX <- as.vector(XX)
    aM <- sirt::sirt_matrix2(x=fixed.a, nrow=LT)
    XX <- aM * XX
  }
  if ( ! is.null(Qmatrix) ){
    XX0 <- tcrossprod( as.matrix(theta), Qmatrix )
    XX <- XX0 - outer( rep(1,nrow(theta)), b )
    XX <- as.vector(XX)
    aM <- sirt::sirt_matrix2(x=fixed.a, nrow=LT)
    XX <- aM * XX
  }
  pm <- sirt::pgenlogis(x=XX, alpha1=alpha1, alpha2=alpha2 )
  pm <- matrix( pm, ncol=length(b))
  return(pm)}
  
#checkIPD <- function(e1,e2,e3,testletStr, difBound, verbose){
#    ### mzp1 vs. mzp2
#       if (verbose) {
#           eq1  <- equat1pl(results=e2, itemF = colnames(e2)[1], valueF = colnames(e2)[2], prmNorm = e1, difBound=difBound, iterativ=TRUE)
#           eq2  <- equat1pl(results=e3, itemF = colnames(e3)[1], valueF = colnames(e3)[2], prmNorm = e2, difBound=difBound, iterativ=TRUE)
#       }  else  {
#           txt1 <- capture.output(eq1  <- equat1pl(results=e2, itemF = colnames(e2)[1], valueF = colnames(e2)[2], prmNorm = e1, difBound=difBound, iterativ=TRUE))
#           txt2 <- capture.output(eq2  <- equat1pl(results=e3, itemF = colnames(e3)[1], valueF = colnames(e3)[2], prmNorm = e2, difBound=difBound, iterativ=TRUE))
#       }
#       if ( "itemExcluded" %in% colnames(eq1[["items"]][["global"]][["global"]][["info"]])) {
#            e1_weg <- setdiff(eq1[["items"]][["global"]][["global"]][["info"]][,"itemExcluded"], "")
#       }  else  {
#            e1_weg <- NULL
#       }
#       if ( "itemExcluded" %in% colnames(eq2[["items"]][["global"]][["global"]][["info"]])) {
#            e3_weg <- setdiff(eq2[["items"]][["global"]][["global"]][["info"]][,"itemExcluded"], "")
#       }  else  {
#            e3_weg <- NULL
#       }
#       return(list(e1_weg = e1_weg, e3_weg = e3_weg))}

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

# equatingList: The object returned by \code{\link{equat1pl}}, for measurement occasion 1 vs. 3!
# multiEquatError_output: The object returned by \code{\link{multiEquatError}}
replaceLinkingError <- function(equatingList, multiEquatError_output, verbose = TRUE, digits = 4) {
    ### checks
       if ( !is.list(multiEquatError_output) && is.null(dim(multiEquatError_output)) ) {stop("'multiEquatError_output' must be a list of data.frames.")}
    ### dimensionen in der Equatingliste. Diese Dimensionen muessen auch in dem 'multiEquatError_output'-Objekt vorhanden sein
       dims <- lapply(equatingList[["items"]], names)
       drin <- unlist(dims) %in% names(multiEquatError_output)
       if(!all(drin)) { stop(paste0("Dimension(s) '",paste(dims[which(drin==FALSE)],collapse="', '"),"' are missing in the 'multiEquatError_output'."))}
       for (d in 1:length(dims)) {
            oldLE <- equatingList[["items"]][[names(dims)[d]]][[dims[[d]]]][["eq"]][["descriptives"]][["linkerror"]]
            if ( verbose) {message(paste0("Dimension '",dims[[d]],"': Replace old linking error ",round(oldLE, digits = digits), " with ", round(multiEquatError_output[[dims[[d]]]][["le3221"]], digits = digits)))}
            equatingList[["items"]][[names(dims)[d]]][[dims[[d]]]][["eq"]][["descriptives"]][["linkerror"]] <- multiEquatError_output[[dims[[d]]]][["le3221"]]
       }
       return(equatingList)}


