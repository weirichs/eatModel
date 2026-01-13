### diverse functions called by defineModel()

doAufb <- function ( x, argL ) {
          if(!is.null(x[["qMatrix"]])) {
     ### check: wenn superSplitter BERUHEND AUF ITEM GROUPING genutzt wird, wird 'items'-Argument von 'defineModel' ignoriert bzw. ueberschrieben
             itemMis<- setdiff ( x[["qMatrix"]][,1], colnames(argL[["dat"]]))
             if(length(itemMis) > 0) {cli::cli_warn(c(paste0("Model split preparation for model No. ",x[["model.no"]], ", model name ",x[["model.name"]],": ", length(itemMis), " items from ",nrow(x[["qMatrix"]])," items listed the Q matrix not found in data:"), "i"=paste0("'", paste(itemMis, collapse= "', '"),"'")))}
             itemSel<- intersect ( x[["qMatrix"]][,1], colnames(argL[["dat"]]))
             qMatrix<- x[["qMatrix"]]
          }  else  {
             if ( is.null(argL[["items"]]) )  { stop(paste0("Model no. ",x[["model.no"]]," ('",x[["model.name"]],"'): no items defined.\n"))}
             itemSel<- argL[["items"]]                                          ### itemSel = "items selected"
          }
     ### Personen im Datensatz selektieren: Achtung: wenn keine Personen in "person.grouping", nimm alle!
          if(!is.null(x[["person.grouping"]])) {
             persMis<- setdiff ( x[["person.grouping"]][,1], argL[["dat"]][,argL[["id"]]])
             if( length ( persMis ) > 0) {warning(paste0("Model split preparation for model No. ",x[["model.no"]], ", model name ",x[["model.name"]],": ", length(persMis) ," from ",nrow(x[["person.grouping"]])," persons not found in data.\n"))}
             persons<- intersect ( x[["person.grouping"]][,1], argL[["dat"]][,argL[["id"]]])
             datSel <- argL[["dat"]][match(persons, argL[["dat"]][,argL[["id"]]]),]
          }  else  {
             datSel <- argL[["dat"]]
          }
     ### Unterverzeichnisse definieren
          if(is.null(argL[["dir"]])) { dirI <- NULL }  else  { dirI   <- file.path(argL[["dir"]], substring(x[["model.subpath"]],3)) }
          nameI  <- x[["model.name"]]
          if ( !is.null(argL[["qMatrix"]]) ) {qMatrix <- argL[["qMatrix"]]}
          if(!exists("qMatrix") && is.null(argL[["qMatrix"]]) ) {
             nDim    <- 1
             qMatrix <- NULL
          }   else   {
             if(!exists("qMatrix") ) { qMatrix <- argL[["qMatrix"]]}
             nDim <- ncol(qMatrix)-1
          }
     ### alles was hier definiert wurde, in ArgL ersetzen
          if(exists("qMatrix")) {argL[["qMatrix"]]<- qMatrix}
          if(exists("itemSel")) {argL[["items"]]  <- itemSel}
          if(exists("datSel"))  {argL[["dat"]]    <- datSel}
          if(!is.null(dirI))    {argL[["dir"]]    <- dirI}
          argL[["analysis.name"]] <- nameI
          argL[["nDim"]] <- nDim
     ### Zusatzargumente (definiert ueber add- und cross-Argumente aus splitModels() ) ggf. ergaenzen
          add  <- setdiff(intersect(names(x), names(formals(defineModel))), "qMatrix")
          if(length(add)>0) {
             for(a in add) {argL[[a]] <- x[[a]]}
          }
          return(argL)}

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
       if(!is.null(fsm))  {                                                     ### Achtung: wenn Items identifiers NICHT unique sind (z.B., Item gibt es global und domaenenspezifisch, dann wird jetzt 'fixSlopeMat'
           fsm  <- eatTools::facToChar(fsm)                                     ### auf die Dimension in der Q Matrix angepasst ... das ist nur erlaubt, wenn es ein eindimensionales Modell ist!!
           if(!is.null( slopeMatDomainCol ) ) {
                allV  <- list(slopeMatDomainCol=slopeMatDomainCol , slopeMatItemCol=slopeMatItemCol, slopeMatValueCol =slopeMatValueCol)
                allNam<- c(allNam, lapply(allV, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = fsm, variable=ii)}))
                if ( ncol(qma) != 2) { stop ( "Duplicated item identifiers in 'fixSlopeMat' are only allowed for unidimensional models.\n") }
                mtch  <- eatTools::whereAre( colnames(qma)[2], fsm[, allNam[["slopeMatDomainCol"]]], verbose=FALSE)
                if ( length( mtch) < 2 ) { stop(cat(paste ( "Cannot found dimension '",colnames(qma)[2],"' in 'fixSlopeMat'. Found following values in '",allNam[["slopeMatDomainCol"]],"' column of 'fixSlopeMat': \n    '", paste( sort(unique(fsm[, allNam[["slopeMatDomainCol"]] ])), collapse="', '"),"'.\n",sep="")))}
                fsm   <- fsm[mtch, c(allNam[["slopeMatItemCol"]],allNam[["slopeMatValueCol"]])]
           }
           if(is.null(slopeMatItemCol)) {
              if(ncol(fsm)!=2) {stop("If 'slopeMatItemCol' is not defined, 'fixSlopeMat' must have two colums.")}
              slopeMatItemCol <- allNam[["slopeMatItemCol"]] <- colnames(fsm)[1]
           }
           if(is.null(slopeMatValueCol)) {
              if(ncol(fsm)!=2) {stop("If 'slopeMatValueCol' is not defined, 'fixSlopeMat' must have two colums.")}
              slopeMatValueCol <-  allNam[["slopeMatValueCol"]] <- colnames(fsm)[2]
           }
           cat( "Warning: To date, fixing slopes only works for dichotomous unidimensional or between-item multidimensional models.\n")
           estVar <- TRUE
           weg2   <- setdiff(fsm[,slopeMatItemCol], allNam[["variablen"]])
           if(length(weg2)>0) {
              message(paste0("Following ",length(weg2), " items in matrix for items with fixed slopes ('fixSlopeMat') which are not in dataset:\n   ", paste(weg2, collapse=", "), "\nRemove these item(s) from 'fixSlopeMat' matrix."))
              fsm <- fsm[-match(weg2,fsm[,slopeMatItemCol]),]
     ### Achtung: wenn nach dem Entfernen der nicht im Datensatz enthaltenen Items keine Items in fsm uebrig bleiben, muss slopeMa auf NULL und irtmodel auf 1pl gesetzt werden
              if(nrow(fsm)==0) {
                 message("No items left in 'fixSlopeMat' if items which do not occur in data are removed. Set 'fixSlopeMat' to NULL and 'irtmodel' to '1PL'")
                 return(list(allNam=allNam, estVar=FALSE, slopMat = NULL, irtmodel = "1PL"))
              }
           }
           weg3   <- setdiff(allNam[["variablen"]], fsm[,slopeMatItemCol])
           if(length(weg3)>0) {
              cat(paste("Following ",length(weg3), " items in dataset without fixed slopes in 'fixSlopeMat'. Slope(s) will be estimated freely.\n",sep=""))
              cat("   "); cat(paste(weg3, collapse=", ")); cat("\n")
           }
     ### Achtung, grosser Scheiss: wenn man nicht (wie oben) eine Reihenfolgespalte angibt, aendert die untere 'by'-Schleife die Sortierung!
           if ( nrow(fsm) != length(unique(fsm[,slopeMatItemCol])) ) { stop ( "Item identifiers in 'fixSlopeMat' are not unique.\n")}
           fsm[,"reihenfolge"] <- 1:nrow(fsm)
           dims  <- (1:ncol(qma))[-1]                                           ### Slopematrix muss itemweise zusammengebaut werden
     ### Achtung: fuer partial credit brauche ich jetzt die Anzahl der Kategorien des Items mit den meisten Kategorien!
           ncatMx<- max(as.numeric(names(eatTools::tableUnlist(dat[,allNam[["variablen"]]])))) + 1
           slopMa<- do.call("rbind", by ( data = fsm, INDICES = fsm[,"reihenfolge"], FUN = function (zeile ) {
                    zeile <- zeile[,-ncol(zeile)]
                    stopifnot ( nrow(zeile) == 1 )
                    qSel  <- qma[which( qma[,1] == zeile[[1]]),]
                    zeilen<- ncatMx * length(dims)                              ### fuer jedes Items gibt es [max. Anzahl Kategorien] * [Anzahl Dimensionen] Zeilen in der TAM matrix
                    block <- cbind ( rep ( match(zeile[[1]], allNam[["variablen"]]), times = zeilen), rep ( 1:ncatMx, each = length(dims) ), rep ( 1:length(dims), times = ncatMx), rep(0, zeilen))
                    matchD<- which ( qSel[,-1] != 0 )
                    stopifnot ( length( matchD ) == 1)
                    match2<- intersect(which(block[,2] == min(block[,2])+1), which(block[,3] == (matchD)))
                    stopifnot ( length( na.omit(match2 )) == 1)
                    block[match2,4] <- as.numeric(zeile[[slopeMatValueCol]])
     ### unteres jetzt nur fuer partial credit
                    if(ncatMx > 2) {
                       block[which(block[,4] == 0)[-1],4] <- 1
                       block[,4] <- block[,4] * (block[,2] - 1)
                    }
                    return(block) }))
       }  else  {
           estVar   <- FALSE
           slopMa   <- NULL
       }
       return(list(allNam=allNam, estVar=estVar, slopMat = slopMa, irtmodel=irtmodel))}


### ----------------------------------------------------------------------------

prepEstSlopegroupsTAM <- function(esg, allNam){
       esgNm <- NULL                                                            ### initialisieren
       if(!is.null(esg))  {
           weg1 <- setdiff(allNam[["variablen"]], esg[,1])
           if(length(weg1)>0) {stop("Items in dataset which are not defined in design matrix for item groups with common slopes ('est.slopegroups').\n")}
           weg2 <- setdiff(esg[,1], allNam[["variablen"]])
           if(length(weg2)>0) {
              cat(paste("Following ",length(weg2), " Items in design matrix for item groups with common slopes ('est.slopegroups') which are not in dataset:\n",sep=""))
              cat("   "); cat(paste(weg2, collapse=", ")); cat("\n")
              cat("Remove these item(s) from design matrix.\n")
              esg <- esg[-match(weg2,esg[,1]),]
           }                                                                    ### untere zeile: pruefen, ob keine fehlenden oder leeren Eintraege in der Liste sind
           weg3 <- c(which(is.na(esg[,2])), which(esg[,2] ==""))
           if(length(weg3)>0) {stop("Items in 'est.slopegroups' with missing or empty values.\n")}
           esgNm<- as.numeric(as.factor(as.character(esg[match(allNam[["variablen"]], esg[,1]),2])))
       }
       return(list(esg=esg, esgNm=esgNm))}

### ----------------------------------------------------------------------------

prepAnchorTAM <- function (dfm, skeleton = NULL) {                              ### dfm = defineModelObject
        ank <- NULL                                                             ### initialisieren
        if(!is.null(dfm[["anchor"]][["ank"]])) {
            ank         <- dfm[["anchor"]][["ank"]]
            allNam      <- dfm[["anchor"]][["allNam"]]
            if(ncol(ank) != 2 && dfm[["irtmodel"]] %nin% c("PCM", "GPCM", "GPCM.groups")) {stop("Anchor parameter frame must have two columns for non-PCM models without specifying item and/or domain column.")}
            notInData   <- setdiff(ank[,1], allNam[["variablen"]])              ### Untere Zeile: Wichtig! Sicherstellen, dass Reihenfolge der Items in Anker-Statement der Reihenfolge im datensatz entspricht
            if(length(notInData)>0)  {ank <- ank[-match(notInData, ank[,1]),]}  ### messages entfernt, denn die werden ja schon in anker() durch mergeAttr() ausgegeben
            if(dfm[["irtmodel"]] %in% c("PCM", "GPCM", "GPCM.groups") && !is.null(skeleton)) {
               ankLong  <- ank |> dplyr::mutate(name = paste(item,category, sep="_"))
               weg      <- which(rownames(skeleton) %nin% ankLong[,"name"])     ### partial credit anchoring using skeleton
               if(length(weg)>0) {skeleton <- skeleton[-weg,]}
               stopifnot(length(ankLong[,"name"]) == length(unique(ankLong[,"name"])))
               stopifnot(length(rownames(skeleton)) == length(unique(rownames(skeleton))))
               ank        <- skeleton
               ank[,"xsi"]<- ankLong[match(rownames(ank), ankLong[,"name"]),"parameter"]
            } else {                                                            ### conventional anchoring
               ank[,1]    <- match(as.character(ank[,1]), allNam[["variablen"]])
            }
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
        snodes <- NULL; QMC <- NULL                                             ### initialisieren
        if(software == "mirt") {
           cat("Specifying 'method' and 'nodes' not yet implemented for 'mirt'.\n")
           return(list(method=method, nodes=nodes, snodes=snodes, QMC=QMC))
        }
        if(method == "quasiMontecarlo" && software == "conquest") {
           cat("Method 'quasiMontecarlo' is not available for software 'conquest'. Set method to 'montecarlo'.\n")
           method <- "montecarlo"
        }
        if(method %in% c("montecarlo", "quasiMontecarlo"))  {                   ### stochastische Integration
           if(is.null(nodes) )   {
              cat(paste("'",method,"' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 1000.\n",sep=""))
              nodes <- 1000
           } else {
              if (nodes < 500 ) {
                  cat(paste0("Warning: Due to user specification, only ",nodes," nodes are used for '",method,"' estimation. Please note or re-specify your analysis.\n"))
              }
           }                                                                    ### untere Zeile: Reihenfolge ist wichtig! erst snodes auf nodes setzen, dann nodes auf NULL
           if ( software == "tam" )   {snodes <- nodes; nodes <- NULL; QMC <- as.logical(car::recode ( method, "'montecarlo'=FALSE; 'quasiMontecarlo'=TRUE"))}
        }  else {                                                               ### numerische Integration
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

anker <- function(lab, prm, qMatrix, domainCol, itemCol, valueCol, catCol)  {
                  stopifnot(ncol(lab)==2)
                  if(!ncol(prm) == 2 ) {                                        ### wenn itemliste nicht unique ... 'domain'-Spalte kann ausgelassen werden
                     if(is.null(itemCol) || is.null(valueCol))  { stop("If anchor parameter frame has more than two columns, 'itemCol' and 'valueCol' must be specified.\n")}
                     allVars <- list(domainCol = domainCol, itemCol=itemCol, valueCol=valueCol, catCol=catCol)
                     allNams <- lapply(allVars, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = prm, variable=ii)})
                     notIncl <- setdiff ( colnames(qMatrix)[-1], prm[,allNams[["domainCol"]]])
                     if(length(notIncl) > 0 && !is.null(allNams[["domainCol"]]) ) {stop(paste ( "Q matrix contains domain(s) ",paste("'",paste(notIncl, collapse="', '"),"'",sep="")," which are not included in the '",allNams[["domainCol"]],"' column of the anchor parameter frame.\n",sep=""))}
                     weg     <- setdiff ( unique(prm[,allNams[["domainCol"]]]), colnames(qMatrix)[-1])
                     if(length(weg) > 0 ) {
                        ind <- eatTools::whereAre ( weg, prm[,allNams[["domainCol"]]], verbose=FALSE)
                        cat(paste("Remove ",length(ind)," rows from the anchor parameter frame which do not belong to any of the specified domains in the Q matrix.\n",sep=""))
                        prm <- prm[-ind,]
                     }
                     prm     <- prm[,c(allNams[["itemCol"]], allNams[["valueCol"]], allNams[["catCol"]], allNams[["domainCol"]])]
                  }
                  colnames(prm) <- c("item", "parameter", "category", "domain")[1:ncol(prm)]
                  if("category" %in% colnames(prm)) {                           ### wenn es die Spalte "category" gibt, darf sie nur Werte enthalten, die "Cat1", "Cat2", ... lauten
                     if(!all(unique(prm[,"category"]) %in% paste0("Cat", 1:9))) {stop(paste0("Invalid values in 'category' column: '",paste( setdiff(unique(prm[,"category"]), paste0("Cat", 1:9)), collapse="', '"), "'"))}
                  }
                  nru <- nrow(unique(prm))                                      ### nru = nrow unique
                  if(nru != nrow(prm)) {
                     dopp<- which(duplicated(prm[,"item"]))
                     warning(paste0("Found ", length(dopp)," duplicate item identifiers (per domain and/or category) in anchor list. Duplicated entries will be deleted."))
                     prm <- unique(prm)
                  }
                  ind <- intersect(lab[,"item"],prm[,"item"])
                  if(length(ind) == 0) {stop("No common items found in 'anchor' list and data frame.\n")}
                  if(length(ind) > 0)  {cat(paste(length(ind), " common items found in 'anchor' list and data frame.\n",sep="")) }
                  resT<- eatTools::mergeAttr(lab, prm, by = "item", sort = FALSE, all = FALSE, setAttr = FALSE, unitName = "item", xName = "item response data", yName = "anchor list", verbose = c("match", "unique"))
                  res <- resT |> dplyr::arrange(dplyr::across(tidyselect::all_of(intersect(c("itemNr", "category"), names(resT))))) |> dplyr::select(-tidyselect::any_of("item"))
                  return(list ( resConquest = res, resTam = resT[,-2]))}

### ----------------------------------------------------------------------------

desk.irt <- function(daten, itemspalten, reduce=TRUE) {
             daten <- eatTools::makeDataFrame(daten, verbose=FALSE)
             if(!missing(itemspalten)) {daten <- daten[,itemspalten,drop=FALSE]}
             res   <- do.call("rbind", lapply(1:ncol(daten), FUN=function(ii) {
                      pVals<- unlist(lapply(1:max(daten[,ii], na.rm=TRUE), FUN = function(j) {
                              if(j %in% daten[,ii]) {                           ### dieses damit mittendrin fehlende Kategorien nicht mitgezaehlt werden
                                 y <- eatTools::num.to.cat(daten[,ii], j-0.0001, c(0,1))
                                 return(mean(y, na.rm=TRUE))}}))
                      modm <- eatTools::makeDataFrame(model.matrix( as.formula(paste0("~factor(",colnames(daten)[ii], ")-1")), data = daten), verbose=FALSE)
                      res.i<- lapply(modm, table)
                      codes<- eatTools::halveString(names(res.i), "\\.", first=FALSE)[,2]
                      abs  <- unlist(lapply(res.i, FUN = function(j) {j[2]}))
                      res.i <- data.frame(item.nr = ii, item.name = colnames(daten)[ii], cases = length(daten[,ii]),Missing=sum(is.na(daten[,ii])),valid=sum(!is.na(daten[,ii])),category = paste0("Cat",codes), Codes=codes,Abs.Freq=abs,item.p=c(NA, pVals), stringsAsFactors=FALSE)
                      return(res.i)}))
             if(reduce == TRUE)  {
                weg   <- which ( res[,"Codes"] == min(res[,"Codes"]) )
                drin  <- setdiff ( 1:nrow(res), weg )
                zusatz<- setdiff ( res[weg,"item.name"], res[drin,"item.name"])
                if ( length(zusatz)>0) {
                     zusatz <- match ( zusatz, res[,"item.name"])
                     drin   <- unique ( c ( zusatz, drin ))
                     weg    <- setdiff ( 1:nrow(res), drin )
                }
                res <- res[drin,]
             }
             return(res)}

### ----------------------------------------------------------------------------

item.diskrim <- function(daten, itemspalten, streng = TRUE) {
  if(!missing(itemspalten))  {daten <- daten[,itemspalten]}
  trenn <- suppressWarnings(eatTools::pwc(daten))
  if(streng) {return(data.frame(item.name=trenn[,"item"],item.diskrim = trenn[,"partWholeCorr"],stringsAsFactors = FALSE))} else {return(data.frame(item.name=trenn[,"item"],item.diskrim = trenn[,"corr"],stringsAsFactors = FALSE))}}

cleanifySplittedModels <- function (lst, argL) {
  if(length(lst) == 4L & !is.null(lst[["models"]]) &  length(nrow( lst[["models"]]) > 0)>0 ) {
     if(!is.symbol ( argL[["analysis.name"]] ) ) {
        cat(paste("Analysis name is already specified by the 'splitted models' object. User-defined analysis name '",argL[["analysis.name"]],"' will be used as prefix.\n",sep=""))
        lst[["models"]][,"model.name"] <- paste(argL[["analysis.name"]], "_", lst[["models"]][,"model.name"],sep="")
        for ( u in 1:length(lst[["models.splitted"]]) ) { lst[["models.splitted"]][[u]][["model.name"]] <- paste(argL[["analysis.name"]], "_", lst[["models.splitted"]][[u]][["model.name"]],sep="") }
     }
     if(nrow(lst[[1L]])>0) {
        mods   <- intersect(lst[["models"]][,"model.no"], unlist(lapply(lst[["models.splitted"]], FUN = function ( l ) {l[["model.no"]]})))
     }  else  {
        mods <- unlist(lapply(lst[["models.splitted"]], FUN = function ( l ) {l[["model.no"]]}))
     }
  }  else  {
     mods <- unlist(lapply(lst[["models.splitted"]], FUN = function ( l ) {l[["model.no"]]}))
  }
  if(length(mods) == 0) { stop("Inconsistent model specification in 'splitted models'.\n") } else { if(argL[["verbose"]] == TRUE) { cat(paste("\nSpecification of 'qMatrix' and 'person.groups' results in ",length(mods)," model(s).\n",sep="")) } }
  if(!is.null(lst[["nCores"]] ) ) {
     if( lst[["nCores"]] > 1 ) {
        cat(paste ( "Use multicore processing. Models are allocated to ",lst[["nCores"]]," cores.\n",sep=""))
        flush.console()
     }
  }
  return(lst)}

generateConsoleInfo <- function (argL, x) {
    items   <- eatTools::existsBackgroundVariables(dat = argL[["dat"]], variable=argL[["items"]])
    listing <- data.frame ( name = c("Model name", "Number of items", "Number of persons", "Number of dimensions"), wert = as.character(c(x[["model.name"]],length(items), nrow(argL[["dat"]]) , argL[["nDim"]])), stringsAsFactors = FALSE)
    zusatz  <- setdiff(names(x),  c("model.no", "model.name", "model.subpath", "dim", "Ndim", "group", "Ngroup", "qMatrix", "person.grouping"))
    if(length(zusatz)>0) {listing <- rbind(listing, data.frame ( name = zusatz, wert = unlist(x[zusatz]), stringsAsFactors=FALSE))}
    listing[,"leerz"] <- max (nchar(listing[,"name"])) - nchar(listing[,"name"]) + 1
    txt     <- apply(listing, MARGIN = 1, FUN = function ( j ) { paste("\n    ", j[["name"]], ":", paste(rep(" ", times = j[["leerz"]]), sep="", collapse=""), j[["wert"]], sep="")})
    nDots   <- max(nchar(listing[,"name"])) + max(na.omit(nchar(listing[,"wert"]))) + 6
    txtRet  <- capture.output(cat(paste("\n\n",paste(rep("=",times = nDots), sep="", collapse=""),"\nModel No. ",x[["model.no"]], paste(txt,sep="", collapse=""), "\n",paste(rep("=",times = nDots), sep="", collapse=""),"\n\n", sep="")))
    return(txtRet)}

