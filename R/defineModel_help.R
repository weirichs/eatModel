### diverse functions called by defineModel()

doAufb <- function ( m, matchCall, anf, verbose, dir, multicore ) {
  matchL <- match(m, unlist(lapply(matchCall[["splittedModels"]][["models.splitted"]], FUN = function ( l ) { l[["model.no"]] } )))
  mess1  <- NULL
  if(!is.null(matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["qMatrix"]])) {
    if ( !is.null(matchCall[["items"]]) )  {
      if(m == anf) { mess1 <- c(mess1, cat("Warning: 'defineModel' was called using 'splitModels' argument. Model split according to item groups is intended. Item selection is defined \n    via 'splittedModels' object. Hence, 'items' argument is expected to be missed in 'defineModel()' and will be ignored.\n")) }
    }
    itemMis<- setdiff ( matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["qMatrix"]][,1], colnames(matchCall[["dat"]]))
    if( length ( itemMis ) > 0) {
      mess1 <- c(mess1, paste( "Warning! Model No. ",matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["model.no"]], ", model name: '",matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["model.name"]],"': ", length(itemMis) ," from ",nrow(matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["qMatrix"]])," items listed the Q matrix not found in data:\n    ", paste(itemMis,collapse=", "),"\n",sep=""))
    }
    itemSel<- intersect ( matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["qMatrix"]][,1], colnames(matchCall[["dat"]]))
    qMatrix<- matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["qMatrix"]]
  }  else  {
    if ( is.null(matchCall[["items"]]) )  { stop(paste0("Model no. ",m," ('",matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["model.name"]],"'): no items defined.\n"))}
    itemSel<- matchCall[["items"]]                                     ### itemSel = "items selected"
  }
  if(!is.null(matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["person.grouping"]])) {
    persMis<- setdiff ( matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["person.grouping"]][,1], matchCall[["dat"]][,matchCall[["id"]]])
    if( length ( persMis ) > 0) {
      mess1 <- c(mess1, paste0( "Warning: ",length(persMis) ," from ",nrow(matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["person.grouping"]])," persons not found in data.\n"))
    }
    persons<- intersect ( matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["person.grouping"]][,1], matchCall[["dat"]][,matchCall[["id"]]])
    datSel <- matchCall[["dat"]][match(persons, matchCall[["dat"]][,matchCall[["id"]]]),]
  }  else  { datSel <- matchCall[["dat"]] }
  if(is.null(matchCall[["dir"]])) { dirI <- NULL }  else  { dirI   <- file.path(dir, substring(matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["model.subpath"]],3)) }
  nameI  <- matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["model.name"]]
  if ( !is.null(matchCall[["qMatrix"]]) ) {qMatrix <- matchCall[["qMatrix"]]}
  if(!exists("qMatrix") && is.null(matchCall[["qMatrix"]]) ) {
    nDim    <- 1
    qMatrix <- NULL
  }   else   {
    if(!exists("qMatrix") ) { qMatrix <- matchCall[["qMatrix"]]}
    nDim <- ncol(qMatrix)-1
  }
  overwr1<- list( dat=datSel, items = itemSel, qMatrix = qMatrix, analysis.name = nameI, dir = dirI)
  overwrF<- setdiff ( colnames(matchCall[["splittedModels"]][["models"]]), c("model.no", "model.name", "model.subpath", "dim", "Ndim", "group", "Ngroup"))
  if(length(overwrF)>0) {
    notAllow <- setdiff ( overwrF, names(formals(defineModel)))
    if ( length ( notAllow ) > 0 ) {
      if ( m == anf ) {
        mess1 <- c(mess1, paste("Column(s) '",paste(notAllow, collapse = "', '"),"' of 'splittedModels' definition frame do not match arguments of 'defineModel()'. Columns will be ignored.\n", sep=""))
      }
      overwrF <- setdiff (overwrF, notAllow)
    }
    notAllow2<- intersect ( overwrF, names(overwr1))
    if ( length ( notAllow2 ) > 0 ) {
      if ( m == anf ) {
        mess1 <- c(mess1, paste("Column(s) '",paste(notAllow2, collapse = "', '"),"' of 'splittedModels' definition frame are not allowed to be modified by user. Columns will be ignored.\n", sep=""))
      }
      overwrF <- setdiff (overwrF, notAllow2)
    }
    notAllow3<- intersect ( overwrF, names(matchCall))
    if ( length ( notAllow3 ) > 0 ) {
      if ( m == anf ) {
        mess1 <- c(mess1, paste("Column(s) '",paste(notAllow3, collapse = "', '"),"' were defined twice, in <models>$models and 'defineModel'. The latter one will be ignored.\n", sep=""))
      }
    }
    if ( length ( overwrF ) > 0 ) {
      for ( hh in overwrF ) {
        overwr1[[hh]] <- matchCall[["splittedModels"]][["models"]][which(matchCall[["splittedModels"]][["models"]][,"model.no"] == m),hh]
        matchCall[["splittedModels"]][["models.splitted"]][[matchL]][[hh]] <- matchCall[["splittedModels"]][["models"]][which(matchCall[["splittedModels"]][["models"]][,"model.no"] == m),hh]
      }
    }
  }  else  { overwrF <-  NULL }
  zusatz <- setdiff ( setdiff ( names(matchCall), "splittedModels"), names( overwr1))
  if ( length ( zusatz ) > 0 ) { overwr1 <- c(overwr1, matchCall[zusatz]) }
  if(!is.null(matchCall[["items"]])) {allVars<- list(variablen=matchCall[["items"]])}
  if(exists("itemSel"))              {allVars<- list(variablen=itemSel)}
  allNams<- lapply(allVars, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = matchCall[["dat"]], variable=ii)})
  overwr3<- data.frame ( arg = c("Model name", "Number of items", "Number of persons", "Number of dimensions"), eval = as.character(c(matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["model.name"]],length(allNams[["variablen"]]), nrow(datSel) , nDim)), stringsAsFactors = FALSE)
  if ( length ( overwrF) > 0 )  {
    zusatz <- lapply ( overwrF, FUN = function ( y ) { matchCall[["splittedModels"]][["models"]][which(matchCall[["splittedModels"]][["models"]][,"model.no"] == m),y] })
    zusatz <- data.frame ( arg = overwrF, eval = as.character ( zusatz ) )
    overwr3<- rbind ( overwr3, zusatz)
  }
  overwr3[,"leerz"] <- max (nchar(overwr3[,"arg"])) - nchar(overwr3[,"arg"]) + 1
  txt    <- apply(overwr3, MARGIN = 1, FUN = function ( j ) { paste("\n    ", j[["arg"]], ":", paste(rep(" ", times = j[["leerz"]]), sep="", collapse=""), j[["eval"]], sep="")})
  nDots  <- max(nchar(overwr3[,"arg"])) + max(na.omit(nchar(overwr3[,"eval"]))) + 6
  if(verbose == TRUE ) {
    cat(paste("\n\n",paste(rep("=",times = nDots), sep="", collapse=""),"\nModel No. ",m, paste(txt,sep="", collapse=""), "\n",paste(rep("=",times = nDots), sep="", collapse=""),"\n\n", sep=""))
    if(!is.null(mess1)) { cat(mess1)}
  }
  if(is.null ( matchCall[["splittedModels"]][["nCores"]] ) | matchCall[["splittedModels"]][["nCores"]] == 1 ) {
    ret    <- do.call("defineModel", args = overwr1)
  }  else  {
    attr(overwr1[["dat"]], "multicore") <- TRUE
    ret    <- overwr1
  }
  return(ret) }

### ----------------------------------------------------------------------------

identifyConquestFolder <- function () {
  cf <- system.file("extdata", "console_Feb2007.exe", package = "eatModel")
  if ( nchar(cf)==0) {cf <- system.file("exec", "console_Feb2007.exe", package = "eatModel")} else{return(cf)}
  if ( nchar(cf)==0) {
    root <- system.file(package = "eatModel")
    if ( !file.exists(file.path(root, "exec"))) {dir.create(file.path(root, "exec"))}
    if ( !file.exists( system.file("exec", "console_Feb2007.exe", package = "eatModel") )) {
      if ( !file.exists("i:/Methoden/00_conquest_console/console_Feb2007.exe") ) {
        packageStartupMessage("Cannot find conquest 2007 executable file. Please choose manually.")
        fname <- file.choose()
      }  else  {
        fname <- "i:/Methoden/00_conquest_console/console_Feb2007.exe"
      }
      if ( nchar(fname)>0) { foo <- file.copy(from = fname, to = file.path(root, "exec", "console_Feb2007.exe") ) }
    }
  } else {
    return(cf)
  }
  return(file.path(root, "exec", "console_Feb2007.exe"))}

### ----------------------------------------------------------------------------

prepGuessMat <- function(guessMat, allNam){
  if(!is.null(guessMat)) {
    weg1 <- setdiff(allNam[["variablen"]], guessMat[,1])
    if(length(weg1)>0) {cat(paste(length(weg1), " item(s) in dataset which are not defined in guessing matrix. No guessing parameter will be estimated for these/this item(s).\n",sep="")) }
    weg2 <- setdiff(guessMat[,1], allNam[["variablen"]])
    if(length(weg2)>0) {
      cat(paste(length(weg2), " item(s) in guessing matrix missing in dataset. Remove these items from guessing matrix.\n",sep=""))
      guessMat <- guessMat[-match( weg2, guessMat[,1])  ,]
    }
    gues <- guessMat[ match( allNam[["variablen"]], guessMat[,1]) , "guessingGroup"]
    gues[which(is.na(gues))] <- 0
  }  else  { gues <- NULL }
  return(gues)}

### ----------------------------------------------------------------------------

prepFixSlopeMatTAM <- function (fsm, allNam, qma, slopeMatDomainCol, slopeMatItemCol, slopeMatValueCol, dat, irtmodel){
  if(!is.null(fsm))  {
    fsm  <- eatTools::facToChar(fsm)
    if(!is.null( slopeMatDomainCol ) ) {
      allV  <- list(slopeMatDomainCol=slopeMatDomainCol , slopeMatItemCol=slopeMatItemCol, slopeMatValueCol =slopeMatValueCol)
      allNam<- c(allNam, lapply(allV, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = fsm, variable=ii)}))
      if ( ncol(qma) != 2) { stop ( "Duplicated item identifiers in 'fixSlopeMat' are only allowed for unidimensional models.\n") }
      mtch  <- eatTools::whereAre( colnames(qma)[2], fsm[, allNam[["slopeMatDomainCol"]]], verbose=FALSE)
      if ( length( mtch) < 2 ) { stop(cat(paste ( "Cannot found dimension '",colnames(qma)[2],"' in 'fixSlopeMat'. Found following values in '",allNam[["slopeMatDomainCol"]],"' column of 'fixSlopeMat': \n    '", paste( sort(unique(fsm[, allNam[["slopeMatDomainCol"]] ])), collapse="', '"),"'.\n",sep="")))}
      fsm   <- fsm[mtch, c(allNam[["slopeMatItemCol"]],allNam[["slopeMatValueCol"]])]
    }
    warning( "To date, fixing slopes only works for dichotomous unidimensional or between-item multidimensional models.")
    estVar <- TRUE
    weg2   <- setdiff(fsm[,1], allNam[["variablen"]])
    if(length(weg2)>0) {
      message(paste0("Following ",length(weg2), " items in matrix for items with fixed slopes ('fixSlopeMat') which are not in dataset:\n   ", paste(weg2, collapse=", "), "\nRemove these item(s) from 'fixSlopeMat' matrix."))
      fsm <- fsm[-match(weg2,fsm[,1]),]
      if(nrow(fsm)==0) {
        message("No items left in 'fixSlopeMat' if items which do not occur in data are removed. Set 'fixSlopeMat' to NULL and 'irtmodel' to '1PL'")
        return(list(allNam=allNam, estVar=FALSE, slopMat = NULL, irtmodel = "1PL"))
      }
    }
    weg3   <- setdiff(allNam[["variablen"]], fsm[,1])
    if(length(weg3)>0) {
      cat(paste("Following ",length(weg3), " items in dataset without fixed slopes in 'fixSlopeMat'. Slope(s) will be estimated freely.\n",sep=""))
      cat("   "); cat(paste(weg3, collapse=", ")); cat("\n")
    }
    if ( nrow(fsm) != length(unique(fsm[,1])) ) { stop ( "Item identifiers in 'fixSlopeMat' are not unique.\n")}
    fsm[,"reihenfolge"] <- 1:nrow(fsm)
    dims  <- (1:ncol(qma))[-1]
    slopMa<- do.call("rbind", by ( data = fsm, INDICES = fsm[,"reihenfolge"], FUN = function (zeile ) {
      zeile <- zeile[,-ncol(zeile)]
      stopifnot ( nrow(zeile) == 1 )
      qSel  <- qma[which( qma[,1] == zeile[[1]]),]
      anzKat<- length(unique(na.omit(dat[,as.character(zeile[[1]])])))
      zeilen<- anzKat * length(dims)
      block <- cbind ( rep ( match(zeile[[1]], allNam[["variablen"]]), times = zeilen), rep ( 1:anzKat, each = length(dims) ), rep ( 1:length(dims), times = anzKat), rep(0, zeilen))
      matchD<- which ( qSel[,-1] != 0 )
      stopifnot ( length( matchD ) == 1)
      match2<- intersect(which(block[,2] == max(block[,2])), which(block[,3] == (matchD)))
      stopifnot ( length( na.omit(match2 )) == 1)
      block[match2,4] <- as.numeric(zeile[[2]])
      return(block) }))
  }  else  {
    estVar   <- FALSE
    slopMa   <- NULL
  }
  return(list(allNam=allNam, estVar=estVar, slopMat = slopMa, irtmodel=irtmodel))}

### ----------------------------------------------------------------------------

prepEstSlopegroupsTAM <- function(esg, allNam){
  if(!is.null(esg))  {
    weg1 <- setdiff(allNam[["variablen"]], esg[,1])
    if(length(weg1)>0) {stop("Items in dataset which are not defined in design matrix for item groups with common slopes ('est.slopegroups').\n")}
    weg2 <- setdiff(esg[,1], allNam[["variablen"]])
    if(length(weg2)>0) {
      cat(paste("Following ",length(weg2), " Items in design matrix for item groups with common slopes ('est.slopegroups') which are not in dataset:\n",sep=""))
      cat("   "); cat(paste(weg2, collapse=", ")); cat("\n")
      cat("Remove these item(s) from design matrix.\n")
      esg <- esg[-match(weg2,esg[,1]),]
    }
    weg3 <- c(which(is.na(esg[,2])), which(esg[,2] ==""))
    if(length(weg3)>0) {stop("Items in 'est.slopegroups' with missing or empty values.\n")}
    esg  <- as.numeric(as.factor(as.character(esg[match(allNam[["variablen"]], esg[,1]),2])))
  }
  return(esg)}

### ----------------------------------------------------------------------------

prepAnchorTAM <- function (ank, allNam) {
  if(!is.null(ank)) {
    stopifnot(ncol(ank) == 2 )
    notInData   <- setdiff(ank[,1], allNam[["variablen"]])
    if(length(notInData)>0)  {ank <- ank[-match(notInData, ank[,1]),]}
    ank[,1]    <- match(as.character(ank[,1]), allNam[["variablen"]])
  }
  return(ank)}

### ----------------------------------------------------------------------------

personWithoutValidValues <- function (dat, allNam, remove.no.answers){
  if(inherits(try(datL  <- reshape2::melt(data = dat, id.vars = unique(unlist(allNam[-match("variablen", names(allNam))])), measure.vars = allNam[["variablen"]], na.rm=TRUE)  ),"try-error"))  {
    cat("W A R N I N G ! ! !   Error in melting for unknown reasons. Try workaround.\n"); flush.console()
    allHG <- setdiff(unique(unlist(allNam[-match("variablen", names(allNam))])), allNam[["ID"]] )
    stopifnot(length(allHG)>0)
    datL  <- reshape2::melt(data = dat, id.vars = allNam[["ID"]], measure.vars = allNam[["variablen"]], na.rm=TRUE)
    datL  <- merge(datL, dat[,unique(unlist(allNam[-match("variablen", names(allNam))]))], by = allNam[["ID"]], all=TRUE)
  }
  wegNV <- setdiff(dat[,allNam[["ID"]]], unique(datL[,allNam[["ID"]]]))
  perNA <- NULL
  if(length(wegNV)>0)   {
    cat(paste("Found ",length(wegNV)," cases with missings on all items.\n",sep=""))
    perNA<- dat[match(wegNV,dat[,allNam[["ID"]]] ), allNam[["ID"]]]
    if( remove.no.answers == TRUE)  {
      cat("Cases with missings on all items will be deleted.\n")
      dat  <- dat[-match(wegNV,dat[,allNam[["ID"]]] ) ,]
    }
    if( remove.no.answers == FALSE) {
      cat("Cases with missings on all items will be kept.\n")
    }
  }
  return(list(dat=dat, perNA=perNA, datL=datL))}

### ----------------------------------------------------------------------------

adaptMethod <- function(method, software,nodes){
  snodes <- NULL; QMC <- NULL
  if(method == "quasiMontecarlo" && software == "conquest") {
    cat("Method 'quasiMontecarlo' is not available for software 'conquest'. Set method to 'montecarlo'.\n")
    method <- "montecarlo"
  }
  if(method %in% c("montecarlo", "quasiMontecarlo"))  {
    if(is.null(nodes) )   {
      cat(paste("'",method,"' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 1000.\n",sep=""))
      nodes <- 1000
    } else {
      if (nodes < 500 ) {
        warning(paste0("Due to user specification, only ",nodes," nodes are used for '",method,"' estimation. Please note or re-specify your analysis."))
      }
    }
    if ( software == "tam" )   {snodes <- nodes; nodes <- NULL; QMC <- as.logical(car::recode ( method, "'montecarlo'=FALSE; 'quasiMontecarlo'=TRUE"))}
  }  else {
    if(is.null(nodes) )   {
      cat(paste("'",method,"' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 20.\n",sep=""))
      nodes <- 20
    }
    if ( software == "tam" )   { nodes <- seq(-6,6,len=nodes); snodes <- 0; QMC <- FALSE}
  }
  return(list(method=method, nodes=nodes, snodes=snodes, QMC=QMC))}

### ----------------------------------------------------------------------------

.substituteSigns <- function(dat, variable, all.Names = NULL ) {
  if(!is.null(variable)) {
    variableNew <- tolower(gsub("_|\\.|-", "", variable))
    if ( !is.null(all.Names)) {
      if (!is.null(all.Names[["DIF.var"]])) {
        variableNew <- unlist(isLetter(variableNew))
      }
    }
    cols        <- match(variable, colnames(dat))
    return(data.frame(cols=cols, old=variable,new=variableNew, stringsAsFactors = FALSE))
  }
  if(is.null(variable)) {return(data.frame(old=TRUE,new=TRUE))}
}

### ----------------------------------------------------------------------------

normalize.path <- function(string)
{string <- gsub("//","/",string)
string <- gsub("/","//",string)
string <- gsub("//","\\\\",string)
return(string)}

### ----------------------------------------------------------------------------

anker <- function(lab, prm, qMatrix, domainCol, itemCol, valueCol, multicore )  {
  stopifnot(ncol(lab)==2)
  if ( !ncol(prm) == 2 )   {
    if ( is.null(itemCol))  { stop("If anchor parameter frame has more than two columns, 'itemCol' must be specified.\n")}
    if ( is.null(valueCol)) { stop("If anchor parameter frame has more than two columns, 'valueCol' must be specified.\n")}
    allVars <- list(domainCol = domainCol, itemCol=itemCol, valueCol=valueCol)
    allNams <- lapply(allVars, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = prm, variable=ii)})
    notIncl <- setdiff ( colnames(qMatrix)[-1], prm[,allNams[["domainCol"]]])
    if ( length( notIncl ) > 0 ) { stop(paste ( "Q matrix contains domain(s) ",paste("'",paste(notIncl, collapse="', '"),"'",sep="")," which are not included in the '",allNams[["domainCol"]],"' column of the anchor parameter frame.\n",sep="")) }
    weg     <- setdiff ( unique(prm[,allNams[["domainCol"]]]), colnames(qMatrix)[-1])
    if ( length ( weg ) > 0 ) {
      ind <- eatTools::whereAre ( weg, prm[,allNams[["domainCol"]]], verbose=FALSE)
      cat(paste("Remove ",length(ind)," rows from the anchor parameter frame which do not belong to any of the specified domains in the Q matrix.\n",sep=""))
      prm <- prm[-ind,]
    }
    prm     <- prm[,c(allNams[["itemCol"]], allNams[["valueCol"]])]
  }
  colnames(prm) <- c("item","parameter")
  dopp<- which(duplicated(prm[,"item"]))
  if(length(dopp)>0) { cat(paste("W A R N I N G !!   Found ",length(dopp)," duplicate item identifiers in anchor list. Duplicated entries will be deleted.\n",sep="")) ; prm <- prm[which(!duplicated(prm[,"item"])), ] }
  ind <- intersect(lab[,"item"],prm[,"item"])
  if(length(ind) == 0) {stop("No common items found in 'anchor' list and data frame.\n")}
  if(length(ind) > 0)  {cat(paste(length(ind), " common items found in 'anchor' list and data frame.\n",sep="")) }
  if(!is.null(multicore) && multicore == TRUE) {
    txt <- capture.output(resT<- eatTools::mergeAttr(lab, prm, by = "item", sort = FALSE, all = FALSE, setAttr = FALSE, unitName = "item", xName = "item response data", yName = "anchor list", verbose = c("match", "unique")),type="message")
    if(length(txt)>0) { cat(txt, sep="\n")}
  }  else  {
    resT<- eatTools::mergeAttr(lab, prm, by = "item", sort = FALSE, all = FALSE, setAttr = FALSE, unitName = "item", xName = "item response data", yName = "anchor list", verbose = c("match", "unique"))
  }
  res <- data.frame(resT[sort(resT[,2],decreasing=FALSE,index.return=TRUE)$ix,], stringsAsFactors = FALSE)[,-1]
  stopifnot(nrow(res) == length(ind))
  return(list ( resConquest = res, resTam = resT[,-2]))}

### ----------------------------------------------------------------------------

desk.irt <- function(daten, itemspalten, na=NA,percent=FALSE,reduce=TRUE,codebook=list(datei=NULL,item=NULL,value=NULL,lab=NULL, komp=NULL), quiet = FALSE ) {
  daten <- eatTools::makeDataFrame(daten)
  if(!missing(itemspalten)) {daten <- daten[,itemspalten,drop=FALSE]}
  if (is.na(na[1])==FALSE) {
    recode.statement <- paste(na,"= NA",collapse="; ")
    daten            <- data.frame(sapply(daten,FUN=function(ii) {car::recode(ii,recode.statement)}),stringsAsFactors=FALSE)
  }
  specific.codes <- lapply(daten,function(ii){NULL})
  if(!is.null(codebook$datei) & !is.null(codebook$value))  {
    specific.codes <- lapply(as.list(colnames(daten)), FUN=function(ii) {
      codebook$datei[codebook$datei[,codebook$item] == ii,c(codebook$item,codebook$value)] } )
    kein.eintrag   <- which(sapply(specific.codes,FUN=function(ii) {nrow(ii)==0}))
    if(length(kein.eintrag)>0)  {cat(paste(length(kein.eintrag)," item(s) missing in codebook:\n",sep=""))
      cat(paste(colnames(daten)[kein.eintrag],collapse=", ")); cat("\n")}
  }
  results        <- lapply(1:ncol(daten), FUN=function(ii) {
    res.i <- eatTools::tablePattern(x=daten[,ii], pattern=specific.codes[[ii]]$value)
    namen.res.i <- names(res.i)
    if(length(res.i)==0) {
      if(quiet == FALSE ) { cat(paste("Item '",colnames(daten)[ii],"' without any valid values.\n",sep=""))}
      res.i <- 0
      namen.res.i <- NA}
    Label <- NA
    KB <- NA
    if(!is.null(codebook$lab))  {Label <- codebook$datei[codebook$datei[,codebook$item] == colnames(daten)[ii],codebook$lab]}
    if(!is.null(codebook$komp)) {KB    <- codebook$datei[codebook$datei[,codebook$item] == colnames(daten)[ii],codebook$komp]}
    res.i <- data.frame(item.nr = ii, item.name = colnames(daten)[ii], Label = Label, KB = KB, cases = length(daten[,ii]),Missing=sum(is.na(daten[,ii])),valid=sum(!is.na(daten[,ii])),Codes=namen.res.i,Abs.Freq=as.numeric(res.i),Rel.Freq=as.numeric(res.i)/sum(!is.na(daten[,ii])), item.p=mean(na.omit(daten[,ii])), stringsAsFactors=FALSE)
  })
  results        <- do.call("rbind",results)
  if(reduce == TRUE)  {
    weg   <- which ( results[,"Codes"] == min(results[,"Codes"]) )
    drin  <- setdiff ( 1:nrow(results), weg )
    zusatz<- setdiff ( results[weg,"item.name"], results[drin,"item.name"])
    if ( length(zusatz)>0) {
      zusatz <- match ( zusatz, results[,"item.name"])
      drin   <- unique ( c ( zusatz, drin ))
      weg    <- setdiff ( 1:nrow(results), drin )
    }
    results <- results[drin,]
  }
  if(percent == TRUE) {results$Rel.Freq <- 100 * results$Rel.Freq}
  return(results)}

### ----------------------------------------------------------------------------

item.diskrim <- function(daten, itemspalten, streng = TRUE) {
  if(!missing(itemspalten))  {daten <- daten[,itemspalten]}
  trenn <- suppressWarnings(eatTools::pwc(daten))
  if(streng) {return(data.frame(item.name=trenn[,"item"],item.diskrim = trenn[,"partWholeCorr"],stringsAsFactors = FALSE))} else {return(data.frame(item.name=trenn[,"item"],item.diskrim = trenn[,"corr"],stringsAsFactors = FALSE))}}
