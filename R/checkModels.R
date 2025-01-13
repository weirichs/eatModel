### called by splitModels() and/or defineModel().
### Contains all functions that check certain aspects of models.

### called by splitModels() ----------------------------------------------------

checkPersonGroupsConsistency <- function(d){
  d <- eatTools::makeDataFrame(d, verbose=FALSE)
  ### Eintraege in erster Spalte muessen unique sein und duerfen keine missings enthalten
  if(any(is.na(d[,1]))){
    stop("Person identifier in first column of 'person.groups' has missing values.")}
  ### die naechsten checks erfolgen jeweils fuer alle weiteren spalten
  chk1 <- lapply(colnames(d)[-1], FUN = function (x){
    ### gruppierungsvariablen duerfen nicht konstant sein und keine fehlenden Werte haben
    if(length(unique(d[,x])) == 1){
      stop(paste0("Column '",x,"' of 'person.groups' is constant."))}
    if(any(is.na(d[,x]))){
      stop(paste0("Column '",x,"' of 'person.groups' has missing values."))}
  })
  # data frame needs at least 2 columns: Person ID, data, etc.
  checkmate::assert_data_frame(d, min.cols = 2)

  return(d)
}

### called by defineModel() and splitModels() ----------------------------------

checkQmatrixConsistency <- function(qmat, errorWhenNot01 = FALSE){
  if(is.null(qmat)) {return(qmat)}

  qmat <- eatTools::makeDataFrame(qmat, name = "Q matrix", onlyWarn=FALSE)
  if(!inherits(qmat[,1], "character")){
    qmat[,1] <- as.character(qmat[,1])}
  nClass<- sapply(qmat[,-1,drop=FALSE], inherits, what=c("numeric", "integer"))

  # all columns - except the first - must be numeric or integer
  if(!all(nClass)){
    warning(paste0("Found non-numeric indicator column(s) in the Q matrix. Transform column(s) '",
                   paste(colnames(qmat)[ which(nClass==FALSE)+1], collapse = "', '") ,"' into numeric format."))
    qmat <- data.frame(qmat[,1,drop=FALSE], eatTools::asNumericIfPossible(qmat[,-1,drop=FALSE]), stringsAsFactors = FALSE)}

  #' There should only be values 0 and 1 (no missings).
  #' In rare cases (conquest 2pl with fixed Itemladungen) other values than 0/1 are ok,
  #' that's why there's a warning here, instead of an error.
  #' Exception: function is called by splitModels() -> HAS to throw an error with values other than 0/1.
  werte <- eatTools::tableUnlist(qmat[, -1, drop=FALSE], useNA="always")
  if(length(setdiff(names(werte), c("0","1", "NA"))) < 0){
    eval(parse(text=paste0("cli::cli_",ifelse(errorWhenNot01, "abort", "warn"),
                           "(c(\"Expected values for Q matrix are 0 and 1.\", \"",
                           ifelse(errorWhenNot01, "x", "i"),
                           "\"=paste0(\"Found unexpected values: '\", paste(names(werte), collapse= \"', '\"),\"'\")))")))}
  if(werte[match("NA", names(werte))] > 0){
    stop("Missing values in Q matrix.\n")}

  # Indikatorspalten duerfen nicht konstant 0 sein (konstant 1 ginge, das waere dann within item multidimensionality)
  wertes <- lapply(qmat[, -1, drop=FALSE], FUN = function(col) {all(col == 0)})
  konst <- which(wertes == TRUE)
  if(length(konst) > 0){              # sind alle Indikatorspalten ausschliesslich 0 -> Fehler
    if(length(konst) == length(wertes)){
      stop("All indicator columns in Q matrix have 0 values.")
    }
    cat(paste0("Column(s) '",paste(names(konst), collapse = "', '"),
               "' in Q matrix are constant with value 0. Delete column(s).\n"))
    qmat <- qmat[,-match(names(konst), colnames(qmat)), drop=FALSE]
  }

  # no doubled input in item columns.
  doppel <- which(duplicated(qmat[,1]))
  if(length(doppel) > 0){
    cat("Found duplicated elements in the item id column of the q matrix. Duplicated elements will be removed.\n")
    chk  <- table(qmat[,1])                           # es wird hier vorher gecheckt, ob - wenn ein Item zweimal in der Q Matrix
    chk  <- chk[which(chk > 1)]                       # auftritt - es beidemale auf dieselben latenten Dimensionen laedt.
    chkL <- lapply(names(chk), FUN = function(ch){
      qChk <- qmat[which(qmat[,1] == ch),]
      pste <- apply(qChk, 1, FUN = function(x) {paste(x[-1], collapse="")})
      if(!all(pste == pste[1])){
        stop(paste0("Inconsistent Q matrix. Item '", ch, "' occurs ", nrow(qChk),
                    " times with incoherent loading structure: \n",
                    eatTools::print_and_capture(qChk, spaces = 3)))
        }
      })
    qmat <- qmat[!duplicated(qmat[,1]),]
  }

  # delete items, that don't load on any dimension.
  zeilen <- apply(qmat, 1, FUN = function(y) {all(names(table(y[-1])) == "0")} )
  weg    <- which(zeilen == TRUE)
  if(length(weg) > 0){
    cat(paste("Note: Following ",length(weg)," item(s) in Q matrix do not belong to any dimension. Delete these item(s) from Q matrix.\n",
              sep=""))
    cat("    "); cat(paste(qmat[weg,1],collapse=", ")); cat("\n")
    qmat <- qmat[-weg,]
  }
  return(qmat)
}

### called by defineModel() ----------------------------------------------------

checkContextVars <- function(x, varname, type = c("weight", "DIF", "group", "HG"), itemdata, suppressAbort = FALSE, internal = FALSE){
  type <- match.arg(arg = type, choices = c("weight", "DIF", "group", "HG"))
  stopifnot(length(x) == nrow(itemdata))
  if(missing(varname)){
    varname <- "ohne Namen"}
  if(!inherits(x, c("numeric", "integer")) && isTRUE(internal)){
    if (type == "weight"){
      stop(paste(type, " variable has to be 'numeric' necessarily. Automatic transformation is not recommended. Please transform by yourself.\n",sep=""))
    }
    cat(paste(type, " variable has to be 'numeric'. Variable '",varname,"' of class '",class(x),"' will be transformed to 'numeric'.\n",sep=""))
    x <- suppressWarnings(unlist(eatTools::asNumericIfPossible(x = data.frame(x, stringsAsFactors = FALSE), transform.factors = TRUE, maintain.factor.scores = FALSE, force.string = FALSE)))
    if(!inherits(x, "numeric")){
      x <- as.numeric(as.factor(x))}
    cat(paste("    '", varname, "' was converted into numeric variable of ",length(table(x))," categories. Please check whether this was intended.\n",sep=""))
    if(length(table(x)) < 12){
      cat(paste("    Values of '", varname, "' are: ",paste(names(table(x)), collapse = ", "),"\n",sep=""))}
  }

  toRemove <- NULL
  mis      <- length(unique(x))
  if(mis == 0){
    if(suppressAbort == FALSE){
      stop(paste("Error: ",type," Variable '",varname,"' without any values.",sep=""))
    }  else  {
      cat(paste0("Warning: ", type," Variable '",varname,"' without any values. '",varname,"' will be removed.\n"))
      toRemove <- varname
    }
  }
  if(mis == 1 ){
    if(suppressAbort == FALSE){
      stop(paste("Error: ",type," Variable '",varname,"' is a constant.",sep=""))
    }  else  {
      cat(paste0(type," Variable '",varname,"' is a constant. '",varname,"' will be removed.\n"))
      toRemove <- varname
    }
  }
  if(type == "DIF" | type == "group"){
    if(mis > 10 && isTRUE(internal))   {warning(paste0(type," Variable '",varname,"' with more than 10 categories. Recommend recoding."))}
  }

  wegDifMis <- NULL; wegDifConst <- NULL; char <- 1; weg <- which(is.na(1:12)); info <- NULL
  if(is.null(toRemove)){
    char <- max(nchar(as.character(na.omit(x))))
    weg  <- which(is.na(x))
    if(length(weg) > 0){
      warning(paste0("Found ",length(weg)," cases with missing on ",type," variable '",varname,"'. Conquest probably will collapse unless cases are not deleted.\n"))}
    if(type == "DIF"){
      if(mis > 2 && isTRUE(internal)){
        cat(paste(type, " Variable '",varname,"' does not seem to be dichotomous.\n",sep=""))
      }
      y       <- paste0("V", x)
      n.werte <- lapply(itemdata, FUN = function(iii) {by(iii, INDICES = list(y), FUN=table, simplify=FALSE)} )
      completeMissingGroupwise <- data.frame(t(sapply(n.werte, function(ll){
        lapply(ll, FUN = function(uu){
          length(uu[uu > 0])
        })
      })), stringsAsFactors = FALSE)

      for(iii in seq(along=completeMissingGroupwise)){
        missingCat.i <- which(completeMissingGroupwise[,iii] == 0)
        if(length(missingCat.i) > 0){
          cat(paste("Warning: Following ", length(missingCat.i), " items with no values in ", type, " variable '",
                    varname, "', group ", substring(colnames(completeMissingGroupwise)[iii], 2), ": \n", sep=""))
          wegDifMis <- c(wegDifMis, rownames(completeMissingGroupwise)[missingCat.i])
          cat(paste0("   ", paste(rownames(completeMissingGroupwise)[missingCat.i],collapse=", "), "\n"))
          info <- plyr::rbind.fill(info,
                                   data.frame(varname = varname, varlevel = substring(colnames(completeMissingGroupwise)[iii], 2),
                                              nCases = table(y)[colnames(completeMissingGroupwise)[iii]], type = "missing",
                                              vars = rownames(completeMissingGroupwise)[missingCat.i], stringsAsFactors = FALSE))
        }
        constantCat.i <- which(completeMissingGroupwise[,iii] == 1)
        if(length(constantCat.i) > 0){
          cat(paste("Warning: Following ", length(constantCat.i), " items are constants in ", type, " variable '",
                    varname, "', group ", substring(colnames(completeMissingGroupwise)[iii], 2), ":\n",sep=""))
          wegDifConst <- c(wegDifConst, rownames(completeMissingGroupwise)[constantCat.i])
          values <- n.werte[rownames(completeMissingGroupwise)[constantCat.i]]
          values <- lapply(values, FUN = function(v) {v[[colnames(completeMissingGroupwise)[iii]]]} )
          cat(paste0("   ", paste(rownames(completeMissingGroupwise)[constantCat.i],collapse=", "), "\n"))
          info <- plyr::rbind.fill(info,
                                   data.frame(varname = varname, varlevel = substring(colnames(completeMissingGroupwise)[iii], 2),
                                              nCases = table(y)[colnames(completeMissingGroupwise)[iii]], type = "constant",
                                              vars =names(values), value =  sapply(values, names), nValue = unlist(values), stringsAsFactors = FALSE))
        }
      }
    }
  }
  return(list(x = x, char = char, weg = weg, varname=varname, wegDifMis = wegDifMis, wegDifConst = wegDifConst, toRemove = toRemove, info=info))
}

### called by defineModel() ----------------------------------------------------

checkBGV <- function(allNam, dat, software, remove.no.answersHG, remove.vars.DIF.missing, namen.items.weg, remove.vars.DIF.constant){
  weg.dif <- NULL; weg.hg <- NULL; weg.weight <- NULL; weg.group <- NULL
  if(length(allNam[["HG.var"]])>0 || length(allNam[["group.var"]])>0 || length(allNam[["DIF.var"]])>0 || length(allNam[["weight.var"]]) >0 || length(allNam[["add.vars"]]) >0 ) {
    varClass<- sapply(c(allNam[["HG.var"]],allNam[["group.var"]],allNam[["DIF.var"]], allNam[["weight.var"]], allNam[["add.vars"]]),FUN = function(ii) {class(dat[,ii])})
    if ( isFALSE(all(sapply(varClass, length) == 1)) ) {
      fehler <- which(sapply(varClass, length) != 1)
      stop("Following ",length(fehler), " variables with more that one class: \n", eatTools::print_and_capture(varClass[names(fehler)], spaces = 5))
    }
  }
  if(length(allNam[["add.vars"]])>0)  { stopifnot(all(sapply(allNam[["add.vars"]], FUN = function(ii) { inherits(dat[,ii], c("integer", "numeric"))})))}
  if(length(allNam[["HG.var"]])>0)    {
    varClass<- sapply(allNam[["HG.var"]], FUN = function(ii) { inherits(dat[,ii], c("integer", "numeric"))})
    if(!all(varClass)) {
      vnam<- names(varClass)[which(varClass == FALSE)]
      cat(paste("Background variable(s) '",paste(vnam, collapse="', '"),"' of class \n    '",paste(sapply(dat[,vnam, drop=FALSE], class),collapse="', '"),"' will be converted to indicator variables.\n",sep=""))
      ind <- do.call("cbind", lapply ( vnam, FUN = function ( yy ) {
        if ( length(which(is.na(dat[,yy])))>0) { stop(paste0("Found ",length(which(is.na(dat[,yy]))), " missings on background variable '",yy,"'."))}
        dat[,yy] <- eatTools::cleanifyString(dat[,yy])
        newFr <- model.matrix( as.formula (paste("~",yy,sep="")), data = dat)[,-1,drop=FALSE]
        cat(paste("    Variable '",yy,"' was converted to ",ncol(newFr)," indicator(s) with name(s) '",paste(colnames(newFr), collapse= "', '"), "'.\n",sep=""))
        return(newFr) }))
      if(software == "conquest") {
        subNm <- .substituteSigns(dat=ind, variable=colnames(ind))
        if(!all(subNm$old == subNm$new)) {
          sn  <- subNm[which( subNm$old != subNm$new),]
          colnames(ind) <- eatTools::recodeLookup(colnames(ind), sn[,c("old", "new")])
        }
      }
      allNam[["HG.var"]] <- c(setdiff(allNam[["HG.var"]],vnam), colnames(ind))
      if ( length(allNam[["HG.var"]]) > 99 && software == "conquest" ) {
        warning(paste0(length(allNam[["HG.var"]]), " background variables might be problematic in 'Conquest'. Recommend to use 'TAM' instead."))
      }
      dat <- data.frame ( dat, ind, stringsAsFactors = FALSE )
    }
    hg.info <- lapply(allNam[["HG.var"]], FUN = function(ii) {checkContextVars(x = dat[,ii], varname=ii, type="HG", itemdata=dat[,allNam[["variablen"]], drop = FALSE], suppressAbort = TRUE, internal=TRUE )})
    for ( i in 1:length(hg.info)) { dat[, hg.info[[i]][["varname"]] ] <- hg.info[[i]]$x }
    wegVar  <- unlist(lapply(hg.info, FUN = function ( uu ) { uu[["toRemove"]] }))
    if(length(wegVar)>0) { allNam[["HG.var"]] <- setdiff ( allNam[["HG.var"]], wegVar) }
    weg.hg  <- unique(unlist(lapply(hg.info, FUN = function ( y ) {y$weg})))
    if(length(weg.hg)>0) {
      if ( remove.no.answersHG == TRUE ) {
        cat(paste("Remove ",length(weg.hg)," cases with missings on at least one HG variable.\n",sep=""))
      }  else  {
        cat(paste(length(weg.hg)," cases with missings on at least one HG variable will be kept according to 'remove.no.answersHG = FALSE'.\n",sep=""))
        weg.hg <- NULL
      }
    }
  }
  if(length(allNam[["group.var"]])>0)  {
    group.info <- lapply(allNam$group.var, FUN = function(ii) {checkContextVars(x = dat[,ii], varname=ii, type="group", itemdata=dat[,allNam[["variablen"]], drop = FALSE], internal=TRUE)})
    for ( i in 1:length(group.info)) { dat[, group.info[[i]]$varname ] <- group.info[[i]]$x }
    weg.group  <- unique(unlist(lapply(group.info, FUN = function ( y ) {y$weg})))
    if(length(weg.group)>0)  {
      cat(paste("Remove ",length(weg.group)," cases with missings on group variable.\n",sep=""))
    }
  }
  if(length(allNam[["DIF.var"]])>0)  {
    dif.info <- lapply(allNam[["DIF.var"]], FUN = function(ii) {checkContextVars(x = dat[,ii], varname=ii, type="DIF", itemdata=dat[,allNam[["variablen"]], drop = FALSE], internal = TRUE)})
    if ( remove.vars.DIF.missing == TRUE ) {
      for ( uu in 1:length(dif.info)) { if (length(dif.info[[uu]]$wegDifMis) >0) {
        cat(paste("Remove item(s) which only have missing values in at least one group of DIF variable '",dif.info[[uu]]$varname,"'.\n", sep=""))
        namen.items.weg <- c(namen.items.weg,dif.info[[uu]]$wegDifMis) }
      }
    }
    if ( remove.vars.DIF.constant == TRUE ) {
      for ( uu in 1:length(dif.info)) { if (length(dif.info[[uu]]$wegDifConst) >0) {
        cat(paste("Remove item(s) which are constant in at least one group of DIF variable '",dif.info[[uu]]$varname,"'.\n",sep=""))
        namen.items.weg <- c(namen.items.weg,dif.info[[uu]]$wegDifConst) }
      }
    }
    for ( i in 1:length(dif.info)) { dat[, dif.info[[i]]$varname ] <- dif.info[[i]]$x }
    weg.dif  <- unique(unlist(lapply(dif.info, FUN = function ( y ) {y$weg})))
    if(length(weg.dif)>0)  {
      cat(paste("Remove ",length(weg.dif)," cases with missings on DIF variable.\n",sep=""))
    }
  }
  if(length(allNam[["weight.var"]])>0)  {
    if(length(allNam[["weight.var"]])!=1) {stop("Use only one weight variable.")}
    weight.info <- lapply(allNam[["weight.var"]], FUN = function(ii) {checkContextVars(x = dat[,ii], varname=ii, type="weight", itemdata=dat[,allNam[["variablen"]], drop = FALSE], internal = TRUE)})
    for ( i in 1:length(weight.info)) { dat[, weight.info[[i]]$varname ] <- weight.info[[i]]$x }
    weg.weight  <- unique(unlist(lapply(weight.info, FUN = function ( y ) {y$weg})))
    if(length(weg.weight)>0) {
      cat(paste("Remove ",length(weg.weight)," cases with missings on weight variable.\n",sep=""))
    }

  }
  namen.all.hg <- unique(c(allNam[["HG.var"]],allNam[["group.var"]],allNam[["DIF.var"]],allNam[["weight.var"]], allNam[["add.vars"]]))
  weg.all <- unique(c(weg.dif, weg.hg, weg.weight, weg.group))
  perExHG <- NULL
  if(length(weg.all)>0) {
    cat(paste("Remove",length(weg.all),"case(s) overall due to missings on at least one explicit variable.\n"))
    perExHG<- dat[weg.all, allNam[["ID"]] ]
    dat    <- dat[-weg.all,]
  }
  return(list(dat=dat, allNam=allNam, namen.items.weg=namen.items.weg,perExHG=perExHG, namen.all.hg=namen.all.hg))}

### called by defineModel() ----------------------------------------------------

### Hilfsfunktion fuer checkItemConsistency
createNamenItemsWeg <- function (crit, remove) {
  if(remove == TRUE) {
    niw <- names(crit)
    mess<- "Remove these items from the data set: "
  } else {
    niw  <- NULL
    mess <- "These items are nevertheless kept in the data set: "
  }
  return(list(niw=niw, mess=mess))}

### Hilfsfunktion fuer defineModel
checkItemConsistency <- function(dat, allNam, remove.missing.items, verbose, removeMinNperItem, minNperItem, remove.constant.items, model.statement){
  options(warn=1)                                                       ### alle Warnungen in dieser Funktion sollen immer angezeigt werden
  namen.items.weg <- NULL                                               ### initialisieren
### Wandle NaN in NA, falls es welche gibt
  is.NaN <- do.call("cbind", lapply(dat[,allNam[["variablen"]], drop = FALSE], FUN = function (uu) { is.nan(uu) } ) )
  if(sum(is.NaN) > 0 ) {
    cat(paste("Found ",sum(is.NaN)," 'NaN' values in the data. Convert 'NaN' to 'NA'.\n",sep=""))
    for ( j in allNam[["variablen"]]) {
      weg <- which ( is.nan(dat[,j] ))
      if(length(weg)>0) {  dat[weg,j] <- NA }
    }
  }
### sind die responses numerisch bzw. stehen da Ziffern drin? (notfalls sowas wie as.character(1) )
  datL    <- dplyr::mutate_at(reshape2::melt(dat, measure.vars = allNam[["variablen"]], id.vars = allNam[["ID"]], na.rm=TRUE), .vars = "value", .funs = as.character)
  zahl    <- grep("[[:digit:]]", datL[,"value"])                        ### sind das alles Ziffern? (auch wenn die Spalten als "character" klassifiziert sind)
  noZahl  <- setdiff(1:nrow(datL), zahl)
  if (length( noZahl ) > 0 ) {
    itemNoZ <- unique(datL[noZahl,"variable"])
    cli::cli_warn(c(paste0("Found {length(noZahl)} non-numeric value{?s} in ",length(itemNoZ)," of ",length(allNam[["variablen"]])," items:"),"i" = paste0("Items: '", paste( itemNoZ, collapse= "', '"), "'"),"i" = paste0("Non-numeric values: '", paste( unique(datL[noZahl,"value"]), collapse= "', '"), "'")))
  }
  klasse  <- unlist( lapply(dat[,allNam[["variablen"]], drop = FALSE], class) )
  if(any(unlist(lapply(dat[,allNam[["variablen"]], drop = FALSE], inherits, what=c("integer", "numeric"))) == FALSE)) {
    warn <- c(unlist(lapply(setdiff(unique(klasse),c("integer", "numeric")), FUN = function (kls) {paste0(length(names(klasse)[which(klasse == kls)])," item columns of class '",kls, "': '",paste(names(klasse)[which(klasse == kls)], collapse="', '"), "'")})),"All item columns will be transformed to be 'numeric'. Recommend to edit your data manually prior to analysis")
    names(warn) <- rep("i", length(warn))
    cli::cli_warn(c("Found unexpected class type(s) in item response columns:", warn))
    dat  <- dplyr::mutate_at(dat, .vars = allNam[["variablen"]], .funs = eatTools::asNumericIfPossible, force.string = TRUE)
  }
  values  <- lapply(dat[,allNam[["variablen"]], drop = FALSE], FUN = function ( ii ) { table(ii)})
  isDichot<- unlist(lapply(values, FUN = function ( vv ) { identical(c("0","1"), names(vv)) }))
  n.werte <- sapply(values, FUN=function(ii) {length(ii)})
  n.mis   <- which(n.werte == 0)
### identifiziere Items ohne jegliche gueltige Werte
  if(length(n.mis) >0) {
    weg  <- createNamenItemsWeg(n.mis, remove = remove.missing.items)
    namen.items.weg <- c(namen.items.weg, weg[["niw"]])
    cli::cli_warn(c(paste0("{length(n.mis)} testitem{?s} without any values:"), "i"=paste0(weg[["mess"]], "'", paste(names(n.mis), collapse="', '"),"'")))
  }
### identifiziere Items mit Anzahl gueltiger Werte < minNperItem
  nValid <- unlist(lapply(dat[,allNam[["variablen"]], drop = FALSE], FUN = function ( ii ) { length(na.omit ( ii )) }))
  below  <- which ( nValid < minNperItem )                              ### identifiziere Items mit weniger gueltigen Werte als in 'minNperItem' angegeben (nur wenn 'removeMinNperItem' = TRUE)
  if(length(below) > 0 ) {
    weg  <- createNamenItemsWeg(below, remove = removeMinNperItem)
    namen.items.weg <- c(namen.items.weg, weg[["niw"]])
    cli::cli_warn(c(paste0("{length(below)} testitem{?s} with less than ", minNperItem, " valid responses."), "i"=paste0(weg[["mess"]], "'", paste(names(below), collapse="', '"),"'")))
  }
### identifiziere konstante Items (Items ohne Varianz)
  constant <- which(n.werte == 1)
  if(length(constant) >0) {
    weg             <- createNamenItemsWeg(constant, remove = remove.constant.items)
    namen.items.weg <- c(namen.items.weg, weg[["niw"]])
    uniqueVal       <- sapply(names(constant), FUN = function (ii) {unique(na.omit(dat[,ii]))})
    nVal            <- sapply(names(constant), FUN = function (ii) {length(which(!is.na(dat[,ii])))})
    cli::cli_warn(c(paste0("{length(constant)} testitem{?s} {?is/are} constants. ", weg[["mess"]]), "i"=paste(paste0("Item '", names(constant), "', only value '", uniqueVal, "' occurs: ", nVal, " valid responses."), sep="\n")))
  }
### identifiziere alle Items, die nicht dichotom (="ND") sind
  n.rasch  <- which( !isDichot )                                        ### (aber nicht die, die bereits wegen konstanter Werte aussortiert wurden!)
  if(length(n.rasch) >0 )   {                                           ### also polytome Items oder Items, die mit 1/2 anstatt 0/1 kodiert sind
    valND <- values[ which(names(values) %in% names(n.rasch)) ]
    valND <- valND[which(sapply(valND, length) > 1)]
    if(length(valND)>0) {
      cat(paste("Warning: ",length(valND)," variable(s) are not strictly dichotomous with 0/1.\n",sep=""))
      for (ii in 1:length(valND))  {
        max.nchar <-  max(nchar(names(table(dat[,names(valND)[ii]]))))
        if(max.nchar>1) {
          cat(paste("Arity of variable",names(valND)[ii],"exceeds 1.\n"))
        }
        if(verbose == TRUE) {
          cat(paste(names(valND)[ii],": ", paste( names(table(dat[,names(valND)[ii]])),collapse=", "),"\n",sep=""))
        }
      }
      cat("Expect a rating scale model or partial credit model.\n")
      if(model.statement == "item") { warning("Sure you want to use 'model statement = item' even when items are not dichotomous?")}
    }
  }
  options(warn=0)
  return(list(dat=dat,allNam=allNam, namen.items.weg=unique(namen.items.weg)))}


### called by defineModel() ----------------------------------------------------

checkID_consistency <- function(dat, allNam, software){
  dat[,allNam[["ID"]] ] <- as.character(dat[,allNam[["ID"]] ])
  doppelt     <- which(duplicated(dat[,allNam[["ID"]]]))
  if(length(doppelt)>0)  {stop(paste( length(doppelt) , " duplicate IDs found!",sep=""))}
  if(software == "conquest") {
    notAllowed  <- grep("-|\\.", dat[,allNam[["ID"]] ])
    if ( length(notAllowed)>0) {
      cat("Conquest neither allows '.' nor '-' in ID variable. Delete signs from ID variable.\n")
      dat[,allNam[["ID"]] ] <- eatTools::removePattern(string = eatTools::removePattern(string=dat[,allNam[["ID"]] ], pattern="\\."), pattern = "-")
      if ( length ( which(duplicated(dat[,allNam[["ID"]]])))>0) {
        dat[,allNam[["ID"]] ] <- paste0(1:nrow(dat),dat[,allNam[["ID"]] ])
      }
    }
  }
  return(dat)}

### called by defineModel() ----------------------------------------------------

checkDir <- function(dir, software) {
  if(!is.null(dir)) {
    dir <- eatTools::crop(dir,"/")
    if(dir.exists(dir) == FALSE) {
      cat(paste("Warning: Specified folder '",dir,"' does not exist. Create folder ... \n",sep=""))
      dir.create(dir, recursive = TRUE)
    }
  }  else  {
    if (software == "conquest") {stop("Argument 'dir' must be specified if software = 'conquest'.\n")}
  }
  return(dir)}

### called by defineModel() ----------------------------------------------------

checkBoundary <- function(dat, allNam, boundary, remove.boundary) {
  datL.valid  <- reshape2::melt(dat, id.vars = allNam[["ID"]], measure.vars = allNam[["variablen"]], na.rm=TRUE)
  if(nrow(datL.valid) == 0) {warning("No valid item values. Skip data preparation."); return(NULL)}
  nValid      <- table(datL.valid[,allNam[["ID"]]])
  inval       <- nValid[which(nValid<boundary)]
  if(length(inval)>0) {
    if ( length( inval > 5)) {auswahl  <- sort ( inval)[c(1, round(length(inval)/2)  ,length(inval))] }  else { auswahl <- sort (inval)[c(1, 3 , length(inval))] }
    cat(paste( length(inval), " subject(s) with less than ",boundary," valid item responses: ", paste(names(auswahl),auswahl,sep=": ", collapse="; ")," ... \n",sep=""))
    if(remove.boundary==TRUE) {
      cat(paste("subjects with less than ",boundary," valid responses will be removed.\n    Caution! This can result in loosing some items likewise.\n",sep="") )
      weg <- match(names(inval), dat[,allNam[["ID"]]])
      stopifnot(length(which(is.na(weg))) == 0 ) ; flush.console()
      dat <- dat[-weg,]
    }
  }
  return(dat)}

### called by defineModel() ----------------------------------------------------

checkPersonSumScores <- function(datL, allNam, dat, remove.failures){
  minMax<- do.call("rbind", by ( data = datL, INDICES = datL[,"variable"], FUN = function ( v ) {
    v[,"valueMin"] <- min(v[,"value"])
    v[,"valueMax"] <- max(v[,"value"])
    return(v)}))
  datW  <- reshape2::dcast(minMax, as.formula(paste(allNam[["ID"]], "~variable",sep="")), value.var = "value")
  datMin<- reshape2::dcast(minMax, as.formula(paste(allNam[["ID"]], "~variable",sep="")), value.var = "valueMin")
  datMax<- reshape2::dcast(minMax, as.formula(paste(allNam[["ID"]], "~variable",sep="")), value.var = "valueMax")
  allFal<- datW[ which ( rowSums ( datW[,-1], na.rm = TRUE ) == rowSums ( datMin[,-1], na.rm = TRUE ) ), allNam[["ID"]] ]
  allTru<- datW[ which ( rowSums ( datW[,-1], na.rm = TRUE ) == rowSums ( datMax[,-1], na.rm = TRUE ) ), allNam[["ID"]] ]
  per0  <- NULL; perA <- NULL
  if(length(allFal)>0) {
    num <- rowSums(datMax[ which ( datMax[,1] %in% allFal), -1], na.rm = TRUE)
    numF<- data.frame ( id = allFal, itemsVisited = num)
    numF<- data.frame(numF[sort(numF[,"itemsVisited"],decreasing=FALSE,index.return=TRUE)$ix,])
    if ( nrow( numF) > 5) { auswahl  <- numF[c(1, round(nrow(numF)/2), nrow(numF)),] }  else { auswahl <- na.omit(numF[c(1, 2, nrow(numF)),]) }
    cat(paste( length(allFal), " subject(s) do not solve any item:\n   ", paste(auswahl[,"id"], " (",auswahl[,"itemsVisited"]," false)",sep="",collapse=", ")," ... \n",sep=""))
    weg0<- na.omit(match(allFal, dat[,allNam[["ID"]]]))
    per0<- data.frame ( numF, itemsSolved = 0, stringsAsFactors = FALSE)
    if (isTRUE(remove.failures))  {
      cat("   Remove subjects without any correct response.\n"); flush.console()
      dat <- dat[-weg0,]
    }
  }
  if(length(allTru)>0) {
    num <- rowSums(datMax[ which ( datMax[,1] %in% allTru), -1], na.rm = TRUE)
    numT<- data.frame ( id = allTru, itemsVisited = num, itemsSolved = num)
    numT<- data.frame(numT[sort(numT[,"itemsSolved"],decreasing=FALSE,index.return=TRUE)$ix,])
    if ( nrow( numT) > 5) { auswahl  <- numT[c(1, round(nrow(numT)/2), nrow(numT)),] }  else { auswahl <- na.omit(numT[c(1, 2, nrow(numT)),]) }
    cat(paste( length(allTru), " subject(s) solved each item: ", paste(auswahl[,"id"], " (",auswahl[,"itemsSolved"] ," correct)",sep="", collapse=", ")," ... \n",sep=""))
    perA<- numT
  }
  return(list(dat=dat, per0=per0, perA=perA))}

