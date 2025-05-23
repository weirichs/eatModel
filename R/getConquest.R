### getConquestResult() is called by getResults()
### the other functions are called by getConquestResult()

getConquestResults<- function(path, analysis.name, model.name, qMatrix, all.Names, abs.dif.bound , sig.dif.bound, p.value, deskRes, discrim, omitFit, omitRegr, omitWle, omitPV, daten, Q3=Q3, q3theta=q3theta, q3MinObs =  q3MinObs, q3MinType = q3MinType, omitUntil) {
  allFiles <- list.files(path=path, pattern = analysis.name, recursive = FALSE)
  qL       <- reshape2::melt(qMatrix, id.vars = colnames(qMatrix)[1], variable.name = "dimensionName", na.rm=TRUE)
  qL       <- qL[which(qL[,"value"] != 0 ) , ]
  varName  <- colnames(qMatrix)[1]
  ret      <- NULL
  logFile  <- paste(analysis.name, "log", sep=".")
  isConv   <- converged ( dir = path, logFile = logFile )
  isPoly   <- length(unique(deskRes[,"Codes"]))>1
  plotPdf  <- getConquestDeviance(path=path, analysis.name = analysis.name, omitUntil = omitUntil)
  ret      <- rbind(ret, getConquestItn (model.name=model.name, analysis.name=analysis.name, qMatrix=qMatrix, qL=qL, allFiles=allFiles, isPoly=isPoly, path=path))
  ret      <- rbind(ret, getConquestDesc (model.name=model.name, deskRes = deskRes, qMatrix=qMatrix, qL = qL, isPoly=isPoly))
  ret      <- rbind(ret, getConquestDiscrim (model.name=model.name, discrim = discrim, qMatrix=qMatrix, qL = qL))
  shwFile  <- paste(analysis.name, "shw", sep=".")
  if (!shwFile %in% allFiles) {
    cat("Cannot find Conquest showfile.\n")
  } else {
    fle  <- file.path(path, shwFile)
    attr(fle, "allNames") <- all.Names
    shw  <- get.shw( file = fle )
    if(is.null( dim(shw$cov.structure) )) {from <- NA} else { from <- shw$cov.structure[-ncol(shw$cov.structure),1]}
    altN <- data.frame ( nr = 1:(ncol(qMatrix)-1), pv = paste("dim", 1:(ncol(qMatrix)-1),sep="."), from = from ,  to = colnames(qMatrix)[-1], stringsAsFactors = FALSE)
    shw[["item"]]  <- merge(shw[["item"]], qL[,-match("value", colnames(qL))], by.x = "item", by.y = colnames(qMatrix)[1], all=TRUE)
    shw12<- getConquestShw (model.name=model.name, qMatrix=qMatrix, qL=qL, shw=shw, altN=altN)
    ret  <- rbind(ret, shw12[["shw1"]], shw12[["shw2"]])
    ret  <- rbind(ret, getConquestInfit (model.name=model.name, shw=shw))
    ret  <- rbind(ret, getConquestAdditionalTerms (model.name=model.name, qMatrix=qMatrix, shw=shw, shwFile = shwFile))
    ret  <- rbind(ret, data.frame ( model = model.name, source="conquest", var1=NA, var2=NA,type="tech", indicator.group="persons", group = colnames(qMatrix)[-1], par="eap", derived.par = "rel", value = shw[["reliability"]][,"eap.rel"], stringsAsFactors=FALSE))
    ret  <- rbind(ret, getConquestPopPar (model.name=model.name, qMatrix=qMatrix, shw=shw))
    ret  <- rbind(ret, getConquestRegPar (model.name=model.name, shw=shw, altN = altN))
    ret  <- rbind(ret, data.frame ( model = model.name, source = "conquest", var1 = NA, var2 = NA , type = "model", indicator.group = NA, group = NA, par = c("deviance", "Npar"),  derived.par = NA, value = shw$final.deviance , stringsAsFactors = FALSE))
    wles <- getConquestWles (model.name=model.name, analysis.name=analysis.name, qMatrix=qMatrix, allFiles=allFiles, omitWle = omitWle, altN = altN, path=path)
    ret  <- rbind(ret, wles[["res"]])
    pvs  <- getConquestPVs (model.name=model.name, analysis.name=analysis.name, omitPV = omitPV, altN = altN, path=path, allFiles=allFiles)
    ret  <- rbind(ret, pvs[["res"]])
    ret  <- rbind(ret, getConquestQ3 (model.name=model.name, shw=shw,Q3=Q3, q3theta=q3theta, omitWle=omitWle, omitPV=omitPV, pv=pvs[["pv"]],wle=wles[["wle"]],daten=daten,all.Names=all.Names, q3MinObs=q3MinObs, q3MinType=q3MinType, shw1 = shw12[["shw1"]]))
  }
  if(!is.null(ret)) {
    attr(ret, "isConverged") <- isConv
    attr(ret, "available")   <- list ( itn =  paste(analysis.name, "itn", sep=".") %in% allFiles, shw =  paste(analysis.name, "shw", sep=".") %in% allFiles, wle = ( paste(analysis.name, "wle", sep=".") %in% allFiles) & (omitWle == FALSE), pv = ( paste(analysis.name, "pvl", sep=".") %in% allFiles) & (omitPV == FALSE))
  }
  return(ret)}

### ----------------------------------------------------------------------------

converged<- function (dir, logFile) {
  isConv <- TRUE
  if (!file.exists(file.path ( dir, logFile ))) {
    cat(paste0("Warning: Model seems not to have converged. Cannot find log file '",file.path ( dir, logFile ),"'.\n"))
    isConv <- FALSE
  }  else  {
    logF  <- scan(file = file.path ( dir, logFile ), what="character",sep="\n",quiet=TRUE)
    if(length(logF) == 0 ) {
      cat(paste0("Warning: Model seems not to have converged. Log file '",file.path ( dir, logFile ),"' is empty.\n"))
      isConv <- FALSE
    }  else  {
      last  <- logF[length(logF)]
      if ( ! eatTools::crop(last) == "=>quit;" ) {
        if ( length( grep("quit;" , last)) == 0 ) {
          cat(paste0("Warning: Model seems not to have converged. Log file unexpectedly finishs with '",last,"'.\nReading in model output might fail.\n"))
          isConv <- FALSE
        }  }  }  }
  return(isConv)  }

### ----------------------------------------------------------------------------

getConquestItn <- function (model.name, analysis.name, qMatrix, qL, allFiles, isPoly, path){
  itnFile  <- paste(analysis.name, "itn", sep=".")
  if (!itnFile %in% allFiles) {
    cat("Cannot find Conquest itn-file.\n")
    return(NULL)
  } else {
    itn  <- get.itn( file.path(path, itnFile) )
    allID<- c("dif.name", "dif.value", "item.name", "Label")
    drin <- allID[which(allID %in% colnames(itn))]
    itnL <- reshape2::melt(itn, id.vars = drin, measure.vars = "pt.bis", value.name = "ptBis", variable.name = "pointBiserialCorrelation", na.rm=FALSE)
    both <- merge(qL, itnL, by.x = colnames(qMatrix)[1], by.y = "item.name", all=TRUE)
    drin2<- setdiff ( drin, "item.name")
    both[,"var2"] <- apply(X = both, MARGIN = 1, FUN = function ( zeile ) { paste( names ( zeile[drin2]), zeile[drin2], sep="=", collapse= ", ") })
    itn3 <- data.frame ( model = model.name, source = "conquest", var1 = both[,colnames(qMatrix)[1]], var2 = NA , type = "fixed", indicator.group = "items", group = both[,"dimensionName"], par = "ptBis",  derived.par = both[,"var2"], value = as.numeric(both[,"ptBis"]), stringsAsFactors = FALSE)
    if ( isPoly == TRUE ) {
      pval<- reshape2::melt(itn, id.vars = drin, measure.vars = "Rel.Freq", variable.name = " itemP", value.name = "pval", na.rm=FALSE)
      both<- merge(qL, pval, by.x = colnames(qMatrix)[1], by.y = "item.name", all=TRUE)
      dri <- setdiff ( drin, "item.name")
      both[,"var2"] <- apply(X = both, MARGIN = 1, FUN = function ( zeile ) { paste( names ( zeile[dri]), zeile[dri], sep="=", collapse= ", ") })
      itn4 <- data.frame ( model = model.name, source = "conquest", var1 = both[,colnames(qMatrix)[1]], var2 = NA , type = "fixed", indicator.group = "items", group = both[,"dimensionName"], par = "itemP",  derived.par = both[,"var2"], value = as.numeric(both[,"pval"])/100, stringsAsFactors = FALSE)
      itn3 <- rbind(itn3, itn4)
    }
  }
  return(itn3)}

### ----------------------------------------------------------------------------

getConquestShw <- function (model.name, qMatrix, qL, shw, altN){
  shw1 <- data.frame ( model = model.name, source = "conquest", var1 = shw$item[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = shw$item[,"dimensionName"], par = "est",  derived.par = NA, value = as.numeric(shw$item[,"ESTIMATE"]), stringsAsFactors = FALSE)
  shw2 <- data.frame ( model = model.name, source = "conquest", var1 = shw$item[,"item"], var2 = NA , type = "fixed", indicator.group = "items",group = shw$item[,"dimensionName"], par = "est",  derived.par = "se", value = as.numeric(shw$item[,"ERROR"]), stringsAsFactors = FALSE)
  toOff<- shw2[ which(is.na(shw2[,"value"])), "var1"]
  if(length(toOff)>0) {
    shw1[match(toOff, shw1[,"var1"]), "par"] <- "offset"
    shw2  <- shw2[-which(is.na(shw2[,"value"])),]                       ### entferne Zeilen aus shw2, die in der "value"-Spalte NA haben
  }
  return(list(shw1=shw1, shw2=shw2))}

### ----------------------------------------------------------------------------

getConquestDesc <- function ( model.name, deskRes, qMatrix, qL, isPoly){
  shw3 <- shw31 <- NULL
  if(is.null ( deskRes ) ) { return(NULL)}
  deskR<- merge(deskRes, qL[,-match("value", colnames(qL))], by.x = "item.name", by.y = colnames(qMatrix)[1], all=TRUE)
  if ( isPoly == FALSE ) {
    shw3 <- data.frame ( model = model.name, source = "conquest", var1 = deskR[,"item.name"], var2 = NA , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "itemP",  derived.par = NA, value = deskR[,"item.p"], stringsAsFactors = FALSE)
  }
  shw4 <- data.frame ( model = model.name, source = "conquest", var1 = deskR[,"item.name"], var2 = NA , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "Nvalid",  derived.par = NA, value = deskR[,"valid"], stringsAsFactors = FALSE)
  shw4 <- shw4[!duplicated(shw4[,"var1"]),]
  cols <- setdiff ( colnames(deskR)[grep("^item.p", colnames(deskR))], "item.p")
  if ( length ( cols ) > 0 ) {
    colsR <- data.frame ( original = cols, reduziert = eatTools::removePattern ( string = cols, pattern = "item.p.") , stringsAsFactors = FALSE)
    shw31 <- do.call("rbind", apply ( colsR, MARGIN = 1, FUN = function ( zeile ) { data.frame ( model = model.name, source = "conquest", var1 = deskR[,"item.name"], var2 = zeile[["reduziert"]] , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "itemP",  derived.par = NA, value = deskR[,zeile[["original"]]], stringsAsFactors = FALSE) }))
  }
  return(rbind(shw3, shw31, shw4))}

### ----------------------------------------------------------------------------

getConquestDiscrim <- function (model.name, discrim , qMatrix, qL){
  if( is.null(discrim) )  {return(NULL)}
  discR<- merge(discrim, qL[,-match("value", colnames(qL))], by.x = "item.name", by.y = colnames(qMatrix)[1], all=TRUE)
  shw5 <- data.frame ( model = model.name, source = "conquest", var1 = discR[,"item.name"], var2 = NA , type = "fixed", indicator.group = "items", group = discR[,"dimensionName"], par = "itemDiscrim",  derived.par = NA, value = discR[,"item.diskrim"], stringsAsFactors = FALSE)
  return(shw5)}

### ----------------------------------------------------------------------------

getConquestInfit <- function (model.name,  shw){
  res <- rbind(data.frame ( model = model.name, source = "conquest", var1 = shw[["item"]][,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = shw$item[,"dimensionName"], par = "est",  derived.par = "infit", value = as.numeric(shw$item[,"MNSQ.1"]), stringsAsFactors = FALSE),
               data.frame ( model = model.name, source = "conquest", var1 = shw[["item"]][,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = shw$item[,"dimensionName"], par = "est",  derived.par = "outfit", value = as.numeric(shw$item[,"MNSQ"]), stringsAsFactors = FALSE) )
  return(res)}

### ----------------------------------------------------------------------------

getConquestAdditionalTerms <- function(model.name, qMatrix, shw, shwFile){
  if(length(shw) <= 4 )  {  return(NULL)}
  res   <- NULL
  read  <- 2 : (length(shw) - 3)
  for ( i in names(shw)[read] ) {
    cols <- unlist(isLetter(i))
    if( !all(cols %in% colnames(shw[[i]])) ) {
      cat(paste("Cannot identify variable identifier for additional term '",i,"' in file '",shwFile,"'. Skip procedure.\n",sep=""))
    }  else  {
      if(length(cols) == 1 ) {
        var1 <- paste( cols, shw[[i]][,cols],sep="_")
      } else {
        var1 <- unlist(apply(shw[[i]][,cols], MARGIN=1, FUN = function ( y ) {paste ( unlist(lapply ( 1:length(y), FUN = function ( yy ) { paste(names(y)[yy], y[yy],sep="_")})), sep="", collapse = "_X_")  }))
      }
      if(ncol(qMatrix) != 2 ){
        cat(paste0("Warning: Cannot identify the group the term '",i,"' in file '",shwFile,"' belongs to. Insert 'NA' to the 'group' column.\n"))
        gr <- NA
      }  else {
        gr <- colnames(qMatrix)[2]
      }
      vars<- c("ESTIMATE", "MNSQ", "MNSQ.1", "ERROR")
      cls <- sapply(shw[[i]][,vars], inherits, what=c("numeric", "integer"))
      if(!all(cls) ) {
         cat(paste0("Expect column(s) '",paste(vars[which(cls==FALSE)],collapse= "', '"), "' in file '",shwFile,"' (statement '",i,"') to be numeric. Current column format is: '",paste(sapply(shw[[i]][,vars[which(cls==FALSE)]],class), collapse="', '"),"'. Column will be transformed.\n"))
         shw[[i]] <- dplyr::mutate_at(shw[[i]], .vars = names(cls[which(cls==FALSE)]), .funs = eatTools::asNumericIfPossible, maintain.factor.scores = TRUE)
      }
      shwE <- data.frame ( model = model.name, source = "conquest", var1 = var1, var2 = NA , type = "fixed", indicator.group = "items", group = gr, par = "est",  derived.par = NA, value = shw[[i]][,"ESTIMATE"], stringsAsFactors = FALSE)
      shwE2<- data.frame ( model = model.name, source = "conquest", var1 = var1, var2 = NA , type = "fixed", indicator.group = "items", group = gr, par = "est",  derived.par = "infit", value = shw[[i]][,"MNSQ.1"], stringsAsFactors = FALSE)
      shwE3<- data.frame ( model = model.name, source = "conquest", var1 = var1, var2 = NA , type = "fixed", indicator.group = "items", group = gr, par = "est",  derived.par = "outfit", value = shw[[i]][,"MNSQ"], stringsAsFactors = FALSE)
      shwSE<- data.frame ( model = model.name, source = "conquest", var1 = var1, var2 = NA , type = "fixed", indicator.group = "items", group = gr, par = "est",  derived.par = "se", value = shw[[i]][,"ERROR"], stringsAsFactors = FALSE)
      toOff<- shwSE[ which(is.na(shwSE[,"value"])), "var1"]
      if(length(toOff)>0) {
        shwE[match(toOff, shwE[,"var1"]), "par"] <- "offset"
        shwSE <- shwSE[-which(is.na(shwSE[,"value"])),]
      }
      res  <- rbind ( res, shwE, shwE2, shwE3, shwSE)
    }
  }
  return(res)}

### ----------------------------------------------------------------------------

getConquestPopPar <- function(model.name, qMatrix, shw){
  if(ncol(qMatrix) == 2) {
    res  <- data.frame ( model = model.name, source = "conquest", var1 = colnames(qMatrix)[2], var2 = NA , type = "distrpar", indicator.group = NA, group = "persons", par = "var",  derived.par = NA, value = shw$cov.structure, stringsAsFactors = FALSE)
  }  else  {
    stopifnot(nrow(shw$cov.structure) == ncol(qMatrix))
    shw$cov.structure[-nrow(shw$cov.structure),1] <- colnames(qMatrix)[-1]
    cov1 <- shw$cov.structure[,-1]
    cov1[upper.tri(shw$cov.structure[,-1])] <- NA
    cov1 <- data.frame ( shw$cov.structure[,1,drop=FALSE], cov1, stringsAsFactors = FALSE)
    colnames(cov1)[-1] <- cov1[-nrow(cov1),1]
    cov2 <- eatTools::facToChar( dataFrame = reshape2::melt(cov1[-nrow(cov1),], id.vars = colnames(cov1)[1], na.rm=TRUE))
    res  <- data.frame ( model = model.name, source = "conquest", var1 = c(colnames(qMatrix)[-1], cov2[,1]), var2 = c(rep(NA, ncol(qMatrix)-1), cov2[,2]) , type = "random", indicator.group = NA, group = "persons", par = c(rep("var",ncol(qMatrix)-1), rep("correlation", nrow(cov2))) ,  derived.par = NA, value = unlist(c(cov1[nrow(cov1),-1], cov2[,3])) , stringsAsFactors = FALSE)
  }
  return(res)}

### ----------------------------------------------------------------------------

getConquestRegPar <- function ( model.name, shw, altN){
  if(nrow(shw$regression)<=1) {return(NULL)}
  reg  <- shw$regression
  if(!is.null( dim(shw$cov.structure) )) {
    for ( i in 1:nrow(altN)) { colnames(reg) <- gsub(altN[i,"from"], altN[i,"to"], colnames(reg))}
  }  else  {
    index  <- grep("_$", colnames(reg))
    colnames(reg)[index] <- paste(colnames(reg)[index], altN[,"to"], sep="")
  }
  regL <- reshape2::melt(reg, id.vars = colnames(reg)[1], measure.vars = colnames(reg)[-c(1, ncol(reg))], na.rm=TRUE)
  foo  <- data.frame ( do.call("rbind", strsplit(as.character(regL[,"variable"]), "_")), stringsAsFactors = FALSE)
  colnames(foo) <- c("par", "group")
  foo[,"derived.par"] <- car::recode(foo[,"par"], "'error'='se'; else = NA")
  foo[,"par"] <- "est"
  regL <- data.frame ( regL[,-match("variable", colnames(regL)), drop=FALSE], foo, stringsAsFactors = FALSE)
  regL[,"reg.var"] <- car::recode(regL[,"reg.var"], "'CONSTANT'='(Intercept)'")
  res  <- data.frame ( model = model.name, source = "conquest", var1 = regL[,"reg.var"], var2 = NA , type = "regcoef", indicator.group = NA, group = regL[,"group"], par = regL[,"par"],  derived.par = regL[,"derived.par"], value = regL[,"value"] , stringsAsFactors = FALSE)
  return(res)}

### ----------------------------------------------------------------------------

getConquestWles <- function ( model.name, analysis.name, qMatrix, allFiles, omitWle, altN, path){
  wleFile  <- paste(analysis.name, "wle", sep=".")
  if ( omitWle == TRUE ) {return(NULL)}
  if (!wleFile %in% allFiles) {
    cat("Cannot find Conquest WLE file.\n")
    return(NULL)
  }
  wle  <- get.wle( file.path(path, wleFile) )
  res  <- NULL
  for ( i in 1:nrow(altN)) { colnames(wle) <- gsub(  paste(".",altN[i,"nr"],"$",sep=""), paste("_", altN[i,"to"],sep="") , colnames(wle))}
  wleL <- reshape2::melt(wle, id.vars = "ID", measure.vars = colnames(wle)[-c(1:2)], na.rm=TRUE)
  foo  <- data.frame ( eatTools::halveString( as.character(wleL[,"variable"]), pattern = "_"), stringsAsFactors=FALSE)
  colnames(foo) <- c("par", "group")
  foo[,"derived.par"] <- car::recode(foo[,"par"], "'wle'='est'; 'std.wle'='se'; else=NA")
  foo[,"par"]         <- car::recode(foo[,"par"], "'wle'='wle'; 'std.wle'='wle'; 'n.solved'='NitemsSolved'; 'n.total'='NitemsTotal'")
  wleL <- data.frame ( wleL[,-match("variable", colnames(wleL)), drop=FALSE], foo, stringsAsFactors = FALSE)
  wleW <- reshape2::dcast(wleL[which(wleL[,"par"] == "wle"),], ID+group~derived.par, value="value")
  rels <- do.call("rbind", by(wleW, INDICES = wleW[,"group"], FUN = function ( g ) { data.frame (dim = g[1,"group"], rel = 1 - mean(g[,"se"]^2)/var(g[,"est"]), stringsAsFactors = FALSE)}))
  res  <- rbind ( res, data.frame ( model = model.name, source = "conquest", var1 = c(wleL[,"ID"],rep(NA,nrow(rels))), var2 = NA , type = c(rep("indicator",nrow(wleL)), rep("tech",nrow(rels))), indicator.group = "persons", group = c(wleL[,"group"],rels[,"dim"]), par = c(wleL[,"par"],rep("wle",nrow(rels))),  derived.par = c(wleL[,"derived.par"],rep("rel", nrow(rels))), value = c(wleL[,"value"] ,rels[,"rel"]) , stringsAsFactors = FALSE))
  return(list(res=res, wle=wle))   }

### ----------------------------------------------------------------------------

getConquestPVs <- function ( model.name, analysis.name, omitPV, altN, path, allFiles){
  pvFile<- paste(analysis.name, "pvl", sep=".")
  if ( omitPV == TRUE ) {return(NULL)}
  if (!pvFile %in% allFiles) {
    cat("Cannot find Conquest PV file.\n")
    return(NULL)
  }
  pv    <- get.plausible( file.path(path, pvFile), forConquestResults = TRUE )
  rec   <- paste("'",altN[,"pv"] , "' = '" , altN[,"to"], "'" ,sep = "", collapse="; ")
  pv$pvLong[,"variable"] <- car::recode( pv$pvLong[,"variable"], rec)
  res   <- data.frame ( model = model.name, source = "conquest", var1 = pv$pvLong[,"ID"], var2 = NA , type = "indicator", indicator.group = "persons", group = pv$pvLong[,"variable"], par = "pv",  derived.par = paste("pv", as.numeric(pv$pvLong[,"PV.Nr"]),sep=""), value = as.numeric(pv$pvLong[,"value"]) , stringsAsFactors = FALSE)
  eaps  <- reshape2::melt ( data.frame ( pv$pvWide[,"ID", drop=FALSE], pv$eap, stringsAsFactors = FALSE), id.vars = "ID", na.rm=TRUE)
  foo   <- data.frame ( do.call("rbind", strsplit(as.character(eaps[,"variable"]), "_")), stringsAsFactors = FALSE)
  colnames(foo) <- c("par", "group")
  foo[,"derived.par"] <- car::recode(foo[,"par"], "'eap'='est'; 'se.eap'='se'; else=NA")
  foo[,"par"]         <- "eap"
  foo[,"group"]       <- car::recode(tolower(foo[,"group"]), rec)
  res   <- rbind(res, data.frame ( model = model.name, source = "conquest", var1 = eaps[,"ID"], var2 = NA , type = "indicator", indicator.group = "persons", group = foo[,"group"], par = "eap",  derived.par = foo[,"derived.par"], value = eaps[,"value"] , stringsAsFactors = FALSE))
  return(list(res=res, pv=pv))}

### ----------------------------------------------------------------------------

getConquestQ3 <- function(model.name, shw,Q3, q3theta, omitWle, omitPV, pv,wle,daten,all.Names, q3MinObs, q3MinType, shw1){
  if ( Q3 == FALSE ) {return(NULL)}
  if ( q3theta == "pv") {
    if ( omitPV == TRUE ) {
      cat("Cannot compute Q3 if 'omitPV == TRUE' and 'q3theta == \"pv\"'. Skip computation.\n")
      return(NULL)
    }
    theta <- pv[["pvWide"]][,2:3]
  }
  if ( q3theta == "wle") {
    if ( omitWle == TRUE ) {
      cat("Cannot compute Q3 if 'omitWle == TRUE' and 'q3theta == \"wle\"'. Skip computation.\n")
      return(NULL)
    }
    colW  <- grep("^wle", colnames(wle))[1]
    theta <- wle[,c(2,colW)]
  }
  if ( q3theta == "eap") {
    if ( omitPV == TRUE ) {
      cat("Cannot compute Q3 if 'omitPV == TRUE' and 'q3theta == \"eap\"'. Skip computation.\n")
      return(NULL)
    }
    colEAP<- grep("^eap", colnames(pv[["pvWide"]]))[1]
    theta <- pv[["pvWide"]][,c(2,colEAP)]
  }
  drinI <- match( shw[["item"]][,"item"], colnames(daten))
  drinP <- match(theta[,1], daten[,"ID"])
  stopifnot(length(which(is.na(drinP))) == 0 , length(which(is.na(drinI))) == 0 )
  q3.res<- sirt::Q3(dat = daten[drinP,drinI], theta = theta[,2], b = shw[["item"]][,"ESTIMATE"], progress = FALSE)
  nObs  <- NULL
    if ( q3MinObs > 1 ) { nObs <- nObsItemPairs ( responseMatrix = daten[,all.Names[["variablen"]]], q3MinType = q3MinType ) }
  matL  <- reshapeQ3 (mat = q3.res$q3.matrix, q3MinObs = q3MinObs, nObs = nObs)
  if( nrow(matL)== 0) { return(NULL)}
  res   <- data.frame ( model = model.name, source = "conquest", var1 = matL[,"Var1"],  var2 = matL[,"Var2"] , type = "fixed",indicator.group = "items",group = paste(names(table(shw1[,"group"])), collapse="_"), par = "q3", derived.par = NA, value = matL[,"value"] , stringsAsFactors = FALSE)
  return(res)}

### ----------------------------------------------------------------------------

getConquestDeviance <- function ( path, analysis.name, omitUntil = omitUntil) {
    ### erstmal zusaetzliche Informationen (Anzahl nodes etc.) gewinnen
  cqc  <- scan(file.path ( path, paste0(analysis.name, ".cqc")),what="character",sep="\n",quiet=TRUE)
  such <- c("method", "nodes", "converge", "seed")
  ret  <- lapply(such, FUN = function ( su ) {
    indm <- grep(paste0(su, "="), cqc)
    if ( length(indm)>1) {                                         ### schlechter Hotfix, f_nodes entfernen
       hf   <- grep("f_nodes", cqc)
       indm <- setdiff(indm, hf)
    }
    if(length(indm) != 1) {
       cat(paste("Cannot identify '",su,"' from cqc file.\n",sep=""))
       met <- NULL
    }  else  {
       pos1<- nchar(unlist(strsplit(cqc[indm], su))[1])            ### position finden, an der 'method' steht
       pos2<- which(sapply(1:nchar(cqc[indm]), FUN = function(x){ substr(cqc[indm],x,x) == ","}))
       if(length(pos2)==0) {
          pos2<- nchar(cqc[indm])
       } else {
          pos2<- min(pos2[which(pos2>pos1)])
       }
       met <- eatTools::removePattern(substr(cqc[indm], pos1+1, pos2-1), paste0(su,"="))
    }
    return(met)})                                                  ### Zeit als Differenz von cqc und shw file
  names(ret) <- such
  ret  <- Filter(Negate(is.null), ret)
  tme  <- file.info ( file.path ( path, paste0(analysis.name, ".shw")))[["mtime"]] - file.info ( file.path ( path, paste0(analysis.name, ".cqc")))[["mtime"]]
  grDevices::pdf(file = file.path ( path, paste0(analysis.name, "_dev.pdf")), width = 10, height = 7.5)
  plotDevianceConquest ( logFile = list ( path=path, analysis.name=analysis.name, ret=ret, tme=tme), omitUntil = omitUntil)
  grDevices::dev.off() }

### called by getConquestDeviance() --------------------------------------------

plotDevianceConquest <- function ( logFile, omitUntil = 1, reverse = TRUE, change = TRUE ) {
  if ( inherits(logFile, "character")) {lf <- logFile}  else  { lf <- file.path(logFile[["path"]], paste0(logFile[["analysis.name"]], ".log"))}
  input<- scan(lf,what="character",sep="\n",quiet=TRUE)
  ind  <- grep("eviance=", input)
  dev  <- unlist(lapply(input[ind], FUN = function (x) {               ### bei negativer deviance veraenderung erzeugt conquest klammern mit
          brace <- grep("\\(", x)                                      ### niedrigster bisheriger deviance, die muessen jetzt mitsamt der werte darinnen entfernt werden
          if(length(brace)>0) {
              weg <- grep("\\(", unlist(strsplit(x, "")))              ### stelle mit aufgehender klammer finden
              x   <- substr(x, 1, weg-1)
          }
          return(x)}))
  dev  <- data.frame ( lapply(data.frame ( eatTools::halveString(dev, "\\."), stringsAsFactors = FALSE), eatTools::removeNonNumeric), stringsAsFactors = FALSE)
  mat  <- data.frame ( iter = 1:length(ind), dev = as.numeric(paste(dev[,1], dev[,2], sep=".")), stringsAsFactors = FALSE)
  if(omitUntil>0)  {
     dc<- mat[-c(1:omitUntil),2]                                       ### 'dc' = 'deviance chance'
  } else {
     dc<- mat[,2]
  }
  if ( change ){
     dc<- diff(dc)
     yl<- "Deviance Change"                                            ### labels der y-Achse definieren
  } else {
     yl<- "Deviance"
  }
  if(reverse){
     dc<- -1 * dc
  }
  dc   <- data.frame ( nr=omitUntil + 1:length(dc), dc)
  xm   <- ceiling( max(dc[,1])/10 )*10
  xt   <- NULL
  for ( i in c( 1:30 ) ){
        xt <- c ( xt, (xm/10) %% i==0 )
  }
  xt   <- max ( which ( xt ) )
  cex  <- 0.85 - ( length(dc[,1]) / 1000 )
  if ( cex < 0.40 ) {
       cex <- 0.40
  }
  if (inherits(logFile,"list")) {
      titel <- paste0("Deviance Change Plot for model '",logFile[["analysis.name"]],"'\n")
  }  else  {
      titel <- "Deviance Change Plot\n"
  }
  par(mar = c(5.1, 4.1, 8.1, 2.1))                                     ### mehr Abstand zwischen titel und plot, damit mehrzeilige ueberschriften reinpassen
  plot ( dc[,1], dc[,2], type="o",
       main=titel,  xlab="Iteration",
       xlim=c(min(dc[,1]),max(dc[,1])),  xaxp=c(0,xm,xt),
       ylab=yl, pch=20, cex=cex, lwd=0.75, mar = c(5, 4, 10, 2) + 0.1 )
  si   <- devtools::session_info(pkgs = "eatModel")
  si   <- si[["packages"]][which(si[["packages"]][,"package"] == "eatModel"),]
  inf  <- Sys.getenv()
  sysi <- Sys.info()
  sys  <- sessionInfo()
  cpu  <- benchmarkme::get_cpu()
  ram  <- benchmarkme::get_ram()
  stri <- paste0("'eatModel', version ", si[["loadedversion"]], ", build ",si[["date"]], ", user: ",sysi[["user"]], ", computername: ", ifelse(sysi[["sysname"]] == "Linux", sysi[["nodename"]], inf["COMPUTERNAME"]), "\nsystem: ", sys[["running"]], ", cpu: ", cpu[["model_name"]], ", cores: ",cpu[["no_of_cores"]], ", RAM: ",capture.output(ram))
  if (inherits(logFile,"list")) {
      stri <- paste0(paste(names(logFile[["ret"]]), logFile[["ret"]], sep=" = ", collapse = "  |  "), "  |  elapsed time: ", timeFormat(logFile[["tme"]]), "\n" , stri, "\n")
  } else {
      stri <- paste0(stri, "\n")
  }
  graphics::mtext(stri)
  graphics::abline( a=0, b=0 )
  dcr  <- dc[dc[,2]<0,]
  graphics::points( dcr[,1], dcr[,2], pch=20, cex=cex, col="red") }

# funktion soll mal nach CRAN -> eatTools
timeFormat <- function(timediff, digits) {
          if(missing(digits)) {
             if(as.numeric(timediff) < 0.05 ) {
                digits <- 3
             }  else  {
                digits <- 1
             }
          }
          paste(round(as.numeric(timediff), digits = digits), attr(timediff, "units"), sep=" ")}
