### getTamResult() is called by getResults()
### the other functions are called by getTamResult()

getTamResults <- function(runModelObj, omitFit, omitRegr, omitWle, omitPV, nplausible , ntheta , normal.approx, samp.regr, theta.model, np.adj, Q3=Q3, q3MinObs =  q3MinObs, q3MinType = q3MinType,
                          pvMethod , group, beta_groups , level , n.iter , n.burnin, adj_MH , adj_change_MH , refresh_MH, accrate_bound_MH,	sample_integers, theta_init, print_iter , verbose, calc_ic, seed) {
  qMatrix<- attr(runModelObj, "defineModelObj")[["qMatrix"]]
  qL     <- reshape2::melt(qMatrix, id.vars = colnames(qMatrix)[1], variable.name = "dimensionName", na.rm=TRUE)
  qL     <- qL[which(qL[,"value"] != 0 ) , ]
  varName<- colnames(qMatrix)[1]
  if( omitRegr == FALSE && !inherits(runModelObj, "tamBayes")) {
    beg <- Sys.time()
    txt <- capture.output ( regr <- tam.se(runModelObj))
    stopifnot ( nrow(regr$beta) == ncol(attr(runModelObj, "Y") )+1)
    rownames(regr$beta) <- c("(Intercept)", colnames(attr(runModelObj, "Y")))
  } else {
    regr <- NULL
  }
  if ( !inherits(runModelObj, "tamBayes") ) {leseAlles <- TRUE} else {leseAlles <- FALSE}
  ret    <- NULL
  resItem<- getTamItempars(runModelObj=runModelObj, qL=qL, qMatrix=qMatrix, leseAlles = leseAlles)
  ret    <- rbind(ret, resItem[["shw1"]], resItem[["shw2"]])
  ret    <- rbind(ret, getTamDescriptives(runModelObj=runModelObj, qL=qL, qMatrix=qMatrix, leseAlles = leseAlles))
  ret    <- rbind(ret, getTamDiscrim(runModelObj=runModelObj, qL=qL, qMatrix = qMatrix, leseAlles = leseAlles))
  ret    <- rbind(ret, getTam2plDiscrim(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, regr = regr, omitRegr=omitRegr))
  beg    <- Sys.time()
  ret    <- rbind(ret, getTamInfit(runModelObj=runModelObj, qL=qL, qMatrix = qMatrix, leseAlles = leseAlles, omitFit = omitFit, seed=seed))
  ret    <- rbind(ret, getTamPopPar(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles))
  ret    <- rbind(ret, getTamRegPar(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, omitRegr = omitRegr, regr=regr))
  ret    <- rbind(ret, getTamModInd(runModelObj=runModelObj, leseAlles = leseAlles))
  beg    <- Sys.time()
  ret    <- rbind(ret, getTamWles(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, omitWle = omitWle))
  beg    <- Sys.time()
  tamArg <- as.list(match.call(definition = getTamResults))
  weg    <- which(names(tamArg) %in% c("tamobj", "Y", "runModelObj", "qMatrix", "leseAlles", "omitPV", "pvMethod", "omitFit", "omitRegr", "omitWle"))
  if ( length(weg)>0) {tamArg <- tamArg[-weg]}
  tamarg <- list()
  for ( i in 2:length(tamArg)) {
    tamarg[[names(tamArg)[i]]] <- eval(tamArg[[i]])
  }
  retPVs <- getTamPVs ( runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, omitPV = omitPV, pvMethod = pvMethod, tam.pv.arguments = tamarg)
  ret    <- rbind(ret, retPVs)
  ret    <- rbind(ret, getTamEAPs(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles))
  ret    <- rbind(ret, getTamQ3(runModelObj=runModelObj, leseAlles = leseAlles, shw1 = resItem[["shw1"]], Q3=Q3, q3MinObs=q3MinObs, q3MinType=q3MinType))
  return(ret)}

### ----------------------------------------------------------------------------

getTamItempars    <- function(runModelObj, qL, qMatrix, leseAlles) {
  if(leseAlles == FALSE) {return(NULL)}
  if ( is.null(attr(runModelObj, "defineModelObj")[["all.Names"]][["DIF.var"]])) {
    xsis <- merge(data.frame ( item = rownames(runModelObj[["xsi"]]), runModelObj[["xsi"]], stringsAsFactors = FALSE), qL[,-match("value", colnames(qL))],  by.x = "item", by.y = colnames(qMatrix)[1], all = TRUE)
  }  else  {
    xsis <- mergeDimensionIfDIF (dat = data.frame ( item = rownames(runModelObj[["xsi"]]), runModelObj[["xsi"]], stringsAsFactors = FALSE), qmatLong = qL[,-match("value", colnames(qL))], datMergingVar="item", remove = "toMerge")
  }
  shw1 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = xsis[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = xsis[,"dimensionName"], par = "est",  derived.par = NA, value = xsis[,"xsi"], stringsAsFactors = FALSE)
  shw2 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = xsis[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = xsis[,"dimensionName"], par = "est",  derived.par = "se", value = xsis[,"se.xsi"], stringsAsFactors = FALSE)
  if ( !is.null(attr(runModelObj, "defineModelObj")[["all.Names"]][["DIF.var"]])) {
    shw1 <- renameDifParameters (dat=shw1, qmatLong = qL[,-match("value", colnames(qL))])
    shw2 <- renameDifParameters (dat=shw2, qmatLong = qL[,-match("value", colnames(qL))])
  }
  toOff<- shw2[ which(shw2[,"value"] == 0 ), "var1"]
  if(length(toOff)>0) {
    shw1[match(toOff, shw1[,"var1"]), "par"] <- "offset"
    shw2  <- shw2[-which(shw2[,"value"] == 0 ),] }                      ### entferne Zeilen aus shw2, die in der "value"-Spalte NA haben, danach: p-Werte einfuegen
  return(list ( shw1=shw1, shw2=shw2))}

### called by getTamItempars() and getTamInfit ---------------------------------

mergeDimensionIfDIF <- function(dat, qmatLong, datMergingVar, remove) {
  dat[,datMergingVar]    <- as.character(dat[,datMergingVar])
  dat[,"toMerge"] <- eatTools::halveString(dat[,datMergingVar], ":", first=TRUE)[,1]
  dat  <- merge(dat, qmatLong[,c("item", "dimensionName")], by.x = "toMerge", by.y = "item", all.x = TRUE)
  dat  <- dat[,-match(remove, colnames(dat))]
  return(dat)}

### called by getTamItempars()  ------------------------------------------------

renameDifParameters <- function(dat, qmatLong) {
  indD5<- setdiff( 1:nrow(dat), grep(":DIF", dat[,"var1"]))
  indD5<- setdiff( dat[indD5,"var1"], qmatLong[,"item"])
  to   <- eatTools::removePattern(string = indD5, pattern = "DIF_")
  to1  <- eatTools::removeNumeric(to)
  to2  <- eatTools::removeNonNumeric(to)
  indD5<- data.frame ( from = indD5, to = paste(to1, to2, sep="_"), stringsAsFactors = FALSE)
  recSt<- paste("'",indD5[,"from"] , "' = '" , indD5[,"to"],"'", collapse="; ",sep="")
  indD <- grep(":DIF", dat[,"var1"])
  indD2<- eatTools::halveString(dat[indD,"var1"], pattern = ":DIF_", first=TRUE)
  indD3<- eatTools::removeNumeric(indD2[,2])
  indD4<- eatTools::removeNonNumeric(indD2[,2])
  dat[indD,"var1"] <- paste("item_", indD2[,1], "_X_", paste(indD3, indD4, sep="_"), sep="")
  dat[,"var1"]     <- car::recode(dat[,"var1"], recSt)
  return(dat)}

### ----------------------------------------------------------------------------

getTamDescriptives    <- function(runModelObj, qL, qMatrix, leseAlles) {
  if(leseAlles == FALSE || is.null ( attr(runModelObj, "defineModelObj")[["deskRes"]] )) {return(NULL)}
  deskR<- merge(attr(runModelObj, "defineModelObj")[["deskRes"]], qL[,-match("value", colnames(qL))],  by.x = "item.name", by.y = colnames(qMatrix)[1], all = TRUE)
  shw3 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = as.character(deskR[,"item.name"]), var2 = NA , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "itemP",  derived.par = NA, value = deskR[,"item.p"], stringsAsFactors = FALSE)
  shw4 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = as.character(deskR[,"item.name"]), var2 = NA , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "Nvalid",  derived.par = NA, value = deskR[,"valid"], stringsAsFactors = FALSE)
  cols <- setdiff ( colnames(deskR)[grep("^item.p", colnames(deskR))], "item.p")
  if ( length ( cols ) > 0 ) {
    colsR <- data.frame ( original = cols, reduziert = eatTools::removePattern ( string = cols, pattern = "item.p.") , stringsAsFactors = FALSE)
    shw31 <- do.call("rbind", apply ( colsR, MARGIN = 1, FUN = function ( zeile ) { data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = as.character(deskR[,"item.name"]), var2 = zeile[["reduziert"]] , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "itemP",  derived.par = NA, value = deskR[,zeile[["original"]]], stringsAsFactors = FALSE) }))
    return(rbind(shw3, shw4, shw31))
  }  else  {
    return(rbind(shw3, shw4))
  } }

### ----------------------------------------------------------------------------

getTamDiscrim    <- function(runModelObj, qL, qMatrix, leseAlles) {
  if(leseAlles == FALSE || is.null ( attr(runModelObj, "defineModelObj")[["discrim"]] )) {return(NULL)}
  discR<- merge( attr(runModelObj, "defineModelObj")[["discrim"]] , qL[,-match("value", colnames(qL))],  by.x = "item.name", by.y = colnames(qMatrix)[1], all = TRUE)
  shw5 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = discR[,"item.name"], var2 = NA , type = "fixed", indicator.group = "items", group = discR[,"dimensionName"], par = "itemDiscrim",  derived.par = NA, value = discR[,"item.diskrim"], stringsAsFactors = FALSE)
  return(shw5)}

### ----------------------------------------------------------------------------

getTam2plDiscrim <- function(runModelObj, qMatrix, leseAlles, regr, omitRegr) {
  if(leseAlles == FALSE || !attr(runModelObj, "defineModelObj")[["irtmodel"]] %in% c("2PL", "2PL.groups", "GPCM", "3PL") ) {return(NULL)}
  shw6 <- do.call("rbind", lapply (  1 : length ( colnames( qMatrix ) [-1] ) , FUN = function ( dims ) {
    if ( isFALSE(omitRegr) ) {
      obj <- regr[["B"]]
    } else {
      obj <- as.data.frame ( runModelObj[["B"]])
      colnames(obj) <- paste0("B.", gsub("Dim0", "Dim", colnames(obj)))
      obj[,"item"]  <- rownames(obj)
      isNull        <- which(sapply(obj, FUN = function ( x ) { all(x==0)})==TRUE)
      if (length (isNull)>0) {
        obj <- obj[,-isNull]
      }
    }
    cols  <- grep(paste0(".Dim",dims,"$" ), colnames(obj), value=TRUE)
    tamMat<- obj[,c("item",cols)]
    weg   <- which(tamMat[,2] == 0)
    if(length(weg)>0) {tamMat <- tamMat[-weg,]}
    shw6D <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = tamMat[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = colnames(qMatrix)[dims+1], par = "estSlope",  derived.par = NA, value = tamMat[,2], stringsAsFactors = FALSE)
    if (ncol(tamMat) == 3 ) {
      shw6se<- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = tamMat[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = colnames(qMatrix)[dims+1], par = "estSlope",  derived.par = "se", value = tamMat[,3], stringsAsFactors = FALSE)
    }  else  {
      shw6se<- NULL
    }
    return(rbind(shw6D, shw6se)) }))
  return(shw6)}

### ----------------------------------------------------------------------------

getTamInfit    <- function(runModelObj, qL, qMatrix, leseAlles, omitFit, seed) {
  if(leseAlles == FALSE || omitFit == TRUE ) {return(NULL)}
  infit<- tam.fit(runModelObj, progress=FALSE, seed=seed)
  fits <- merge(infit[["itemfit"]], qL[,-match("value", colnames(qL))],  by.x = "parameter", by.y = colnames(qMatrix)[1], all = TRUE)
  if ( is.null(attr(runModelObj, "defineModelObj")[["all.Names"]][["DIF.var"]])) {
    ret  <- rbind(data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = fits[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = fits[,"dimensionName"], par = "est",  derived.par = "infit", value = fits[,"Infit"], stringsAsFactors = FALSE),
                  data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = fits[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = fits[,"dimensionName"], par = "est",  derived.par = "outfit", value = fits[,"Outfit"], stringsAsFactors = FALSE))
  }  else  {
    ret  <- rbind(data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = infit$itemfit[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = NA, par = "est",  derived.par = "infit", value = infit$itemfit[,"Infit"], stringsAsFactors = FALSE),
                  data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = infit$itemfit[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = NA, par = "est",  derived.par = "outfit", value = infit$itemfit[,"Outfit"], stringsAsFactors = FALSE) )
    ret  <- mergeDimensionIfDIF(dat=ret, qmatLong=qL[,-match("value", colnames(qL))], datMergingVar="var1", remove = c("group", "toMerge"))
    colnames(ret) <- car::recode(colnames(ret), "'dimensionName'='group'")
    ret  <- renameDifParameters(dat=ret, qmatLong=qL[,-match("value", colnames(qL))])
  }
  return(ret)}

### ----------------------------------------------------------------------------

getTamPopPar    <- function(runModelObj, qMatrix, leseAlles) {
  if(leseAlles == FALSE ) {return(NULL)}
  if(ncol(qMatrix) == 2) {
    ret  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = colnames(qMatrix)[2], var2 = NA , type = "distrpar", indicator.group = NA, group = "persons", par = "var",  derived.par = NA, value = runModelObj[["variance"]][1,1] , stringsAsFactors = FALSE)
  }  else  {
    cov1 <- runModelObj[["variance"]]
    colnames(cov1) <- colnames(qMatrix)[-1]
    rownames(cov1) <- colnames(qMatrix)[-1]
    cor1 <- cov2cor(cov1)
    for (ii in 1:nrow(cor1))   {
      cor1[ii,ii:ncol(cor1)] <- NA}
    cor1 <- reshape2::melt(cor1, measure.vars = colnames(cor1), na.rm=TRUE)
    vars <- Matrix::diag(cov1)
    ret  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = c(names(vars),as.character(cor1[,"Var1"])) , var2 = c(rep(NA, length(vars)), as.character(cor1[,"Var2"])) , type = "random", indicator.group = NA, group = "persons", par = c(rep("var",length(vars)), rep("correlation", nrow(cor1))) ,  derived.par = NA, value = c(unlist(vars), cor1[,"value"]), stringsAsFactors = FALSE)
  }
  return(ret)}

### ----------------------------------------------------------------------------

getTamRegPar    <- function(runModelObj, leseAlles, qMatrix, omitRegr, regr) {
  if(leseAlles == FALSE || omitRegr == TRUE ) {return(NULL)}
  if( !isTRUE(all.equal ( dim(runModelObj$beta) , c(1,1))))  {
    regr <- data.frame ( reg.var = rownames(regr$beta), regr$beta, stringsAsFactors = FALSE)
    regr <- reshape2::melt(regr, id.vars = "reg.var", na.rm=TRUE)
    regr2<- data.frame ( par = "est", derived.par = car::recode(unlist(lapply(strsplit(as.character(regr[,"variable"]),"\\."), FUN = function ( l ) {l[1]})), "'se'='se'; else=NA"), group = colnames(qMatrix)[as.numeric(eatTools::removePattern( string = unlist(lapply(strsplit(as.character(regr[,"variable"]),"\\."), FUN = function ( l ) {l[2]})), pattern = "Dim")) + 1], regr, stringsAsFactors = FALSE)
    regr3<- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = regr2[,"reg.var"], var2 = NA , type = "regcoef", indicator.group = NA, group = regr2[,"group"], par = regr2[,"par"],  derived.par = regr2[,"derived.par"], value = regr2[,"value"] , stringsAsFactors = FALSE)
  }  else {
    return(NULL)
  }
  return(regr3)}

### ----------------------------------------------------------------------------

getTamModInd    <- function(runModelObj, leseAlles) {
  if(leseAlles == FALSE ) {return(NULL)}
  return(data.frame ( model = attr(runModelObj,"defineModelObj")[[ "analysis.name"]], source = "tam", var1 = NA, var2 = NA , type = "model", indicator.group = NA, group = NA, par = c("deviance", "Npar", "AIC", "BIC"), derived.par = NA, value = unlist(runModelObj[["ic"]][c("deviance", "Npars", "AIC", "BIC")]), stringsAsFactors = FALSE))}

### ----------------------------------------------------------------------------

getTamWles    <- function(runModelObj, qMatrix, leseAlles, omitWle) {
  if(leseAlles == FALSE || omitWle == TRUE ) {return(NULL)}
  txt  <- capture.output(wle  <- tam.wle(runModelObj, progress = FALSE))
  eind1<- ncol(wle) == 7
  if(isTRUE(eind1)) {
    cols <- grep("^PersonScores$|^PersonMax$|^theta$|^error$|^WLE.rel$", colnames(wle))
    stopifnot(length(cols) == 5)
    colnames(wle)[cols] <- paste(colnames(wle)[cols], ".Dim01", sep="")
  }
  weg1 <- grep("WLE.rel", colnames(wle))
  wleL <- reshape2::melt(wle, id.vars = "pid", measure.vars = colnames(wle)[-c(1:2,weg1)], na.rm=TRUE)
  wleL[,"group"] <- colnames(qMatrix)[as.numeric(eatTools::removePattern(string = unlist(lapply(strsplit(as.character(wleL[,"variable"]),"\\."), FUN = function (l) {l[2]})), pattern = "Dim"))+1]
  trans<- na.omit(unique(data.frame ( original = unlist(lapply(strsplit(as.character(wleL[,"variable"]),"\\."), FUN = function (l) {l[2]})), uebersetzt = wleL[,"group"], stringsAsFactors = FALSE)))
  wleL[,"par"]   <- car::recode(unlist(lapply(strsplit(as.character(wleL[,"variable"]),"\\."), FUN = function (l) {l[1]})), "'PersonScores'='NitemsSolved'; 'PersonMax'='NitemsTotal'; 'theta'='wle'; 'error'='wle'")
  wleL[,"derived.par"] <- car::recode(unlist(lapply(strsplit(as.character(wleL[,"variable"]),"\\."), FUN = function (l) {l[1]})), "'theta'='est'; 'error'='se';else=NA")
  rel  <- reshape2::melt(as.data.frame ( wle)[1,weg1], na.rm = TRUE)
  if ( "variable" %in% colnames(rel)) {
    rel[,"original"] <- unlist(lapply(strsplit(as.character(rel[,"variable"]),"\\."), FUN = function (l) {l[length(l)]}))
    rel  <- merge(rel, trans, by="original", all=TRUE)
  }  else  {
    rel[,"uebersetzt"] <- colnames(qMatrix)[-1]
  }
  res  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = c(wleL[,"pid"],rep(NA,nrow(rel))), var2 = NA , type = c(rep("indicator",nrow(wleL)), rep("tech", nrow(rel))), indicator.group = "persons", group = c(wleL[,"group"],rel[,"uebersetzt"]), par = c(wleL[,"par"],rep("wle", nrow(rel))),  derived.par = c(wleL[,"derived.par"],rep("rel",nrow(rel))), value = c(wleL[,"value"] ,rel[,"value"]), stringsAsFactors = FALSE)
  return(res)}

### ----------------------------------------------------------------------------

getTamPVs <- function ( runModelObj, qMatrix, leseAlles, omitPV, pvMethod, tam.pv.arguments) {
  if(omitPV == TRUE ) {return(NULL)}
  for ( i in names( tam.pv.arguments )) { assign(i, tam.pv.arguments[[i]]) }
  if(leseAlles == TRUE ) {
    if ( pvMethod == "regular" ) {
      do   <- paste ( "tam.pv ( ", paste(names(formals(tam.pv)), car::recode ( names(formals(tam.pv)), "'tamobj'='runModelObj'"), sep =" = ", collapse = ", "), ")",sep="")
    } else {
      if ( is.null ( attr(runModelObj, "Y") ) ) {
        warning("Conditioning model was not defined ('Y' is NULL).")
        Y1 <- NULL
      } else {
        Y1 <- data.frame ( intercpt = 1, attr(runModelObj, "Y"))
      }
      do   <- paste ( "tam.pv.mcmc ( ", paste(names(formals(tam.pv.mcmc)), car::recode ( names(formals(tam.pv.mcmc)), "'tamobj'='runModelObj'; 'Y'='Y1'"), sep =" = ", collapse = ", "), ")",sep="")
    }
  }  else  {
    class(runModelObj) <- "list"
    stopifnot ( pvMethod == "bayesian")
    do   <- paste ( "tam.pv.mcmc ( ", paste(names(formals(tam.pv.mcmc)), car::recode ( names(formals(tam.pv.mcmc)), "'tamobj'='runModelObj'; 'Y'='runModelObj[[\"Y\"]]'; 'nplausible'='attr(runModelObj, \"defineModelObj\")[[\"n.plausible\"]]'"), sep =" = ", collapse = ", "), ")",sep="")
  }
  pv   <- eval(parse(text=do))
  pvL  <- reshape2::melt(pv$pv, id.vars = "pid", na.rm=TRUE)
  pvL[,"PV.Nr"] <- as.numeric(eatTools::removePattern(string = unlist(lapply(strsplit(as.character(pvL[,"variable"]),"\\."), FUN = function (l) {l[1]})), pattern = "PV"))
  pvL[,"group"] <- colnames(qMatrix)[as.numeric(eatTools::removePattern(string = unlist(lapply(strsplit(as.character(pvL[,"variable"]),"\\."), FUN = function (l) {l[2]})), pattern = "Dim"))+1]
  res <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = pvL[,"pid"], var2 = NA , type = "indicator", indicator.group = "persons", group = pvL[,"group"], par = "pv",  derived.par = paste("pv", pvL[,"PV.Nr"],sep=""), value = pvL[,"value"] , stringsAsFactors = FALSE)
  return(res)}

### ----------------------------------------------------------------------------

getTamEAPs <- function ( runModelObj, qMatrix, leseAlles = leseAlles) {
  if(leseAlles == FALSE ) {return(NULL)}
  eaps <- runModelObj[["person"]]
  eind1<- ncol(eaps) == 7
  if(eind1 == TRUE) {
    cols <- grep("EAP$", colnames(eaps))
    stopifnot(length(cols) == 2)
    colnames(eaps)[cols] <- paste(colnames(eaps)[cols], ".Dim1", sep="")
  }
  eaps <- reshape2::melt(eaps, id.vars = "pid", measure.vars = grep("EAP", colnames(eaps)), na.rm=TRUE)
  eaps[,"tam"]      <- eatTools::halveString(string = as.character(eaps[,"variable"]), pattern = "\\.", first = FALSE)[,"X2"]
  eaps[,"dimnumber"]<- as.numeric(eatTools::removePattern ( eaps[,"tam"], "Dim"))
  eaps[,"group"]    <- colnames(qMatrix)[eaps[,"dimnumber"] + 1]
  checkmate::assert_numeric(unique(eaps[,"dimnumber"]), len = ncol(qMatrix)-1, any.missing = FALSE)
  eaps[,"par"]      <- "est"
  eaps[grep("^SD.",as.character(eaps[,"variable"])),"par"]   <- "se"
  res  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = eaps[,"pid"], var2 = NA , type = "indicator", indicator.group = "persons", group = eaps[,"group"], par = "eap", derived.par = eaps[,"par"], value = eaps[,"value"] , stringsAsFactors = FALSE)
  if (ncol(qMatrix)>2) {
    grp <-  eatTools::recodeLookup(names(runModelObj[["EAP.rel"]]), unique(eaps[,c("tam", "group")]))
  } else {
    stopifnot(length(unique(eaps[,"group"])) == 1)
    grp <- unique(eaps[,"group"])
  }
  res  <- rbind(res, data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = NA, var2 = NA , type = "tech", indicator.group = "persons", group = grp, par = "eap", derived.par = "rel", value = runModelObj[["EAP.rel"]] , stringsAsFactors = FALSE) )
  return(res)}

### ----------------------------------------------------------------------------

getTamQ3 <- function(runModelObj, leseAlles, shw1, Q3, q3MinObs, q3MinType){
  if(leseAlles == FALSE || Q3 == FALSE) {return(NULL)}
  nObs <- NULL
    if ( q3MinObs > 1 ) {
      nObs <- nObsItemPairs ( responseMatrix = runModelObj[["resp"]], q3MinType = q3MinType )
    }
    mat  <- tam.modelfit ( tamobj = runModelObj, progress = FALSE )
    matL <- reshapeQ3 (mat = mat$Q3.matr, q3MinObs = q3MinObs, nObs = nObs)
    if( nrow(matL)>0) {
      res  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = matL[,"Var1"], var2 = matL[,"Var2"] , type = "fixed",indicator.group = "items", group = paste(names(table(shw1[,"group"])), collapse="_"), par = "q3", derived.par = NA, value = matL[,"value"] , stringsAsFactors = FALSE)
    }  else  {
      res  <- NULL
    }
  return(res)}
