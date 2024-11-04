defineModel <- function(dat, items, id, splittedModels = NULL,
                        irtmodel = c("1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", "2PL.groups", "GPCM.design", "3PL"),
                        qMatrix=NULL, DIF.var=NULL, HG.var=NULL, group.var=NULL,
                        weight.var=NULL, anchor = NULL, domainCol=NULL, itemCol=NULL,
                        valueCol=NULL,check.for.linking = TRUE, minNperItem = 50,
                        removeMinNperItem = FALSE, boundary = 6, remove.boundary = FALSE,
                        remove.no.answers = TRUE, remove.no.answersHG = TRUE,
                        remove.missing.items = TRUE, remove.constant.items = TRUE,
                        remove.failures = FALSE, remove.vars.DIF.missing = TRUE,
                        remove.vars.DIF.constant = TRUE, verbose=TRUE, software = c("conquest","tam"),
                        dir = NULL, analysis.name, schooltype.var = NULL, model.statement = "item",
                        compute.fit = TRUE, pvMethod = c("regular", "bayesian"),
                        fitTamMmlForBayesian = TRUE, n.plausible=5, seed = NULL,
                        conquest.folder=NULL, constraints=c("cases","none","items"),
                        std.err=c("quick","full","none"), distribution=c("normal","discrete"),
                        method=c("gauss", "quadrature", "montecarlo", "quasiMontecarlo"),
                        n.iterations=2000, nodes=NULL, p.nodes=2000, f.nodes=2000,
                        converge=0.001, deviancechange=0.0001, equivalence.table=c("wle","mle","NULL"),
                        use.letters=FALSE, allowAllScoresEverywhere = TRUE, guessMat = NULL,
                        est.slopegroups = NULL, fixSlopeMat = NULL, slopeMatDomainCol=NULL,
                        slopeMatItemCol=NULL, slopeMatValueCol=NULL, progress = FALSE,
                        Msteps = NULL, increment.factor=1 , fac.oldxsi=0,
                        export = list(logfile = TRUE, systemfile = FALSE, history = TRUE,
                                      covariance = TRUE, reg_coefficients = TRUE, designmatrix = FALSE) ){

### checking/asserting the arguments -------------------------------------------
  # dat
  ismc <- attr(dat, "multicore")
  dat  <- eatTools::makeDataFrame(dat, name = "dat")

  # list: splitted Models
  checkmate::assert_list(splittedModels, null.ok = TRUE)
  # subset: irtmodel
  checkmate::assert_subset(irtmodel, choices = c("1PL", "2PL", "PCM", "PCM2","RSM",
                                                 "GPCM", "2PL.groups", "GPCM.design", "3PL"))

  # data frame: qMatrix
  checkQmatrixConsistency(qMatrix)

  #' assert vector:
  #' HG.var, group.var, weight.var, schooltype.var
  lapply(c(HG.var, group.var, weight.var), checkmate::assert_vector, null.ok = TRUE)
  checkmate::assert_vector(schooltype.var, len = 1, null.ok = TRUE)

  #' assert numeric:
  #' minNperItem, boundary, n.plausible, seed, n.iterations, nodes, converge,
  #' deviancechange, Msteps
  lapply(c(minNperItem, boundary, n.plausible, n.iterations, nodes, Msteps),
         checkmate::assert_numeric, len = 1, lower = 0)
  lapply(c(converge, deviancechange), checkmate::assert_numeric, len = 1)
  checkmate::assert_numeric(seed, len = 1, null.ok = TRUE)

  #' assert character via match.arg:
  # software
  software <- match.arg(arg = software, choices = c("conquest","tam"))
  if(software == "conquest" && is.null(conquest.folder)) {
    conquest.folder <- identifyConquestFolder()
  }
  # method
  method <- match.arg(arg = method, choices = c("gauss", "quadrature", "montecarlo", "quasiMontecarlo"))

  #' assert logical:
  lapply(c(check.for.linking, removeMinNperItem, remove.boundary, remove.no.answers,
           remove.no.answersHG, remove.missing.items, remove.constant.items,
           remove.failures, remove.vars.DIF.missing, remove.vars.DIF.constant,
           verbose),
         checkmate::assert_logical, len = 1)

  # arguments specific to `conquest`:
  if(software == "conquest"){
    # character: dir, analysis.name, conquest.folder, model.statement, export
    lapply(c(dir, analysis.name),
           checkmate::assert_character, len = 1)
    checkmate::assert_character(model.statement, len = 1, null.ok = TRUE)
    checkmate::assert_character(export)
    checkmate::assert_directory(conquest.folder)
    # logical: compute.fit, use.letters, allowAllScoresEverywhere
    lapply(c(compute.fit, use.letters, allowAllScoresEverywhere),
           checkmate::assert_logical, len = 1)
    # character via match.arg: constraints, std.err, distribution, equivalence.table
    constraints <- match.arg(arg = constraints, choices = c("cases","none","items"))
    std.err <- match.arg(arg = std.err, choices = c("quick","full","none"))
    distribution <- match.arg(arg = distribution, choices = c("normal","discrete"))
    equivalence.table <- match.arg(arg = equivalence.table, choices = c("wle","mle","NULL"))
    # numeric: p.nodes, f.nodes
    lapply(c(p.nodes, f.nodes), checkmate::assert_numeric, len = 1)
  }

  # arguments specific to `tam`:
  if(software == "tam"){
    # character: dir, analysis.name
    lapply(c(dir, analysis.name), checkmate::assert_character, len = 1, null.ok = TRUE)
    # character via match.arg: pvMethod, constraints
    pvMethod <- match.arg(arg = pvMethod, choices = c("regular", "bayesian"))

    constraints <- match.arg(arg = constraints, choices = c("cases", "none", "items"))
    if(constraints == "none"){
      stop("tbd: You can't use constraints = 'none' when using 'tam'.")
    }
    # logical: fitTamMmlForBayesian, progress
    lapply(c(fitTamMmlForBayesian, progress), checkmate::assert_logical, len = 1)
    # named data frames: guessMat, est.slopegroups, fixSlopeMat
    lapply(c(guessMat, est.slopegroups, fixSlopeMat), checkmate::assert_data_frame, col.names = "named", ncols = 2)
    warning(paste0("The first column of the data frame `guessMat` should be called `item`, but is called ",
                   colnames(guessMat)[1], "."))
    warning(paste0("The first column of the data frame `est.slopegroups` should be called `item`, but is called ",
                   colnames(est.slopegroups)[1], "."))

    lapply(c(guessMat[,1], est.slopegroups[,1], fixSlopeMat[,1]), checkmate::assert_character)
    lapply(c(guessMat[,2], est.slopegroups[,2], fixSlopeMat[,2]), checkmate::assert_numeric)
  }

### function -------------------------------------------------------------------

  if(!is.null(splittedModels)) {
    if(length(splittedModels) == 4L & !is.null(splittedModels[["models"]]) &  length(nrow( splittedModels[["models"]]) > 0)>0 ) {
      if ( !missing ( analysis.name ) ) {
        cat(paste("Analysis name is already specified by the 'splittedModels' object. User-defined analysis name '",analysis.name,"' will be used as prefix.\n",sep=""))
        splittedModels[["models"]][,"model.name"] <- paste(analysis.name, "_", splittedModels[["models"]][,"model.name"],sep="")
        for ( u in 1:length(splittedModels[["models.splitted"]]) ) { splittedModels[["models.splitted"]][[u]][["model.name"]] <- paste(analysis.name, "_", splittedModels[["models.splitted"]][[u]][["model.name"]],sep="") }
      }
      if(nrow(splittedModels[[1L]])>0) {
        mods   <- intersect(splittedModels[["models"]][,"model.no"], unlist(lapply(splittedModels[["models.splitted"]], FUN = function ( l ) {l[["model.no"]]})))
      }  else  {
        mods <- unlist(lapply(splittedModels[["models.splitted"]], FUN = function ( l ) {l[["model.no"]]}))
      }
    }  else  {
      mods <- unlist(lapply(splittedModels[["models.splitted"]], FUN = function ( l ) {l[["model.no"]]}))
    }
    if(length(mods) == 0) { stop("Inconsistent model specification in 'splittedModels'.\n") } else { if(verbose == TRUE) { cat(paste("\nSpecification of 'qMatrix' and 'person.groups' results in ",length(mods)," model(s).\n",sep="")) } }
    if(!is.null(splittedModels[["nCores"]] ) ) {
      if( splittedModels[["nCores"]] > 1 ) {
        cat(paste ( "Use multicore processing. Models are allocated to ",splittedModels[["nCores"]]," cores.\n",sep=""))
        flush.console()
      }
    }
    cl1 <- as.list(match.call(definition = defineModel))
    cl1[["analysis.name"]] <- NULL
    cl1[["dat"]] <- as.name("dat")
    cl1[["splittedModels"]] <- splittedModels
    cll <- list()
    for ( u in 2:length(cl1)) {cll[[u-1]] <- eval(cl1[[u]])}
    names(cll) <- names(cl1)[-1]
    anf <- mods[1]
    if(is.null ( splittedModels[["nCores"]] ) | splittedModels[["nCores"]] == 1 ) {
      models <- lapply ( mods, FUN = doAufb, matchCall = cll, anf=anf, verbose=verbose, dir=dir, multicore=FALSE)
    }  else  {
      txt <- capture.output ( models <- lapply ( mods, FUN = doAufb, matchCall = cll, anf=anf, verbose=verbose, dir=dir, multicore=TRUE) )
      doIt<- function (laufnummer,  ... ) {
        if(!exists("getResults"))  { library(eatModel) }
        txt <- capture.output ( res <- do.call("defineModel", args = models[[laufnummer]] ) )
        return(list ( res=res, txt=txt)) }
      beg <- Sys.time()
      if(splittedModels[["mcPackage"]] == "parallel") {
        cl  <- makeCluster(splittedModels[["nCores"]], type = "SOCK")
      }  else  {
        cl  <- future::makeClusterPSOCK(splittedModels[["nCores"]], verbose=FALSE)
      }
      mods<- clusterApply(cl = cl, x = 1:length(models), fun = doIt)
      stopCluster(cl)
      cat(paste ( length(models), " models were prepared for estimation: ", sep="")); print( Sys.time() - beg, digits = 3)
      models <- lapply(mods, FUN = function ( m ) { m[["res"]] } )
      txts<- lapply(mods, FUN = function ( m ) { m[["txt"]] } )
      luec<- which(txt == "")
      pos <- luec[which ( diff(luec) == 1L )]
      dif2<- which(diff(pos) == 1L)
      if(length(dif2)>0) { pos <- pos [ -dif2 ] }
      pos <- c(pos, length(txt)+1)
      txtP<- lapply ( 1:(length(pos)-1), FUN = function ( u ) { txt[ pos[u] : (pos[u+1]-1) ] })
      txtG<- NULL
      stopifnot(length(txtP) == length(txts))
      for ( j in 1:length(txtP) ) {
        txtG <- c(txtG, txtP[[j]], txts[[j]])
      }
      fl  <- min(grep( pattern = "====", x = txt))
      if ( fl > 1) {
        if (!all(txt[1:(fl-1)] == "" )) {
          txtG <- c("", txt[1:(fl-3)], txtG)
        }
      }
      cat(txtG, sep="\n")
    }
    attr(models, "split") <- splittedModels
    class(models)    <- c("defineMultiple", "list")
    return(models)
  }  else  {
    if ( is.null(Msteps) ) {
      if ( irtmodel == "3PL" ) { Msteps <- 10 } else { Msteps <- 4 }
    }
    method   <- match.arg(method)
    pvMethod <- match.arg(pvMethod)
    if(software == "conquest") {
      original.options <- options("scipen")
      options(scipen = 20)
      if(missing(analysis.name)) {stop("Please specify 'analysis.name' or use 'software = \"tam\"'\n")}
    }  else  {
      if(missing(analysis.name)) {analysis.name <- "not_specified"}
    }
    checkmate::assert_character(model.statement, len = 1)
    if(missing(dat))   {stop("No dataset specified.\n") }
    if(is.null(items)) {stop("Argument 'items' must not be NULL.\n",sep="")}
    if(length(items) == 0 ) {stop("Argument 'items' has no elements.\n",sep="")}
    if ( length(items) != length(unique(items)) ) {
      warning(paste0("Warning: Item identifier is not unique. Only ",length(unique(items))," unique item identifiers will be used."))
      items <- unique(items)
    }
    if(length(id) != 1 ) {stop("Argument 'id' must be of length 1.\n",sep="")}
    if(model.statement != "item") {
      vars <- setdiff(eatTools::crop(unlist(strsplit(model.statement, "\\+|-|\\*"))), "item")
      mis  <- which(!vars %in% colnames(dat))
      if ( length(mis)>0) {stop(paste0("Variables '",paste(vars[mis], collapse="', '"), "' from 'model.statement' not found in data."))}
    }  else  {
      vars <- NULL
    }
    allVars     <- list(ID = id, variablen=items, DIF.var=DIF.var, HG.var=HG.var, group.var=group.var, weight.var=weight.var, schooltype.var = schooltype.var, add.vars = vars)
    all.Names   <- lapply(allVars, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = dat, variable=ii)})
    if(software == "conquest") {
      if(max(nchar(all.Names[["variablen"]]))>11) {stop("In Conquest, maximum length of variable names must not exceed 11 characters. Please shorten variables names.\n")}
    }
    dat <- checkID_consistency(dat=dat, allNam=all.Names, software=software)
    dir <- checkDir(dir=dir, software=software)
    dat <- checkBoundary(dat=dat, allNam=all.Names, boundary=boundary, remove.boundary=remove.boundary)
    subsNam <- .substituteSigns(dat=dat, variable=unlist(all.Names[-c(1:2)]), all.Names = all.Names)
    if(software == "conquest" || !is.null(all.Names[["DIF.var"]])) {
      if(!all(subsNam$old == subsNam$new)) {
        sn     <- subsNam[which( subsNam$old != subsNam$new),]
        message("Conquest neither allows '.', '-', and '_' nor upper case letters in explicit variable names and numbers in DIF variable name. Delete signs from variables names for explicit and DIF variables: \n\n", eatTools::print_and_capture (sn, spaces = 5), "\n")
        colnames(dat) <- eatTools::recodeLookup(colnames(dat), sn[,c("old", "new")])
        all.Names     <- lapply(all.Names, FUN = function ( y ) {eatTools::recodeLookup(y, sn[,c("old", "new")]) })
        if(model.statement != "item") {
          cat("    Remove deleted signs from variables names for explicit variables also in the model statement. Please check afterwards for consistency!\n")
          for ( uu in 1:nrow(sn))  {model.statement <- gsub(sn[uu,"old"], sn[uu,"new"], model.statement)}
        }
      }
      if("item" %in% unlist(all.Names[-c(1:2)])) { stop("Conquest does not allow labelling explicit variable(s) with 'Item' or 'item'.\n") }
    }
    if(length(intersect(all.Names$DIF.var, all.Names$variablen))>0)    {stop("Test items and DIF variable have to be mutually exclusive.\n")}
    if(length(intersect(all.Names$weight.var, all.Names$variablen))>0) {stop("Test items and weighting variable have to be mutually exclusive.\n")}
    if(length(intersect(all.Names$HG.var, all.Names$variablen))>0)     {stop("Test items and HG variable have to be mutually exclusive.\n")}
    if(length(intersect(all.Names$group.var, all.Names$variablen))>0)  {stop("Test items and group variable have to be mutually exclusive.\n")}
    if(is.null(qMatrix)) { qMatrix <- data.frame ( item = all.Names$variablen, Dim1 = 1, stringsAsFactors = FALSE) } else {
      qMatrix <- checkQmatrixConsistency(qMatrix)
      notInDat<- setdiff(qMatrix[,1], all.Names$variablen)
      notInQ  <- setdiff( all.Names$variablen , qMatrix[,1])
      if(length(notInDat)>0) {
        cat(paste("Following ", length(notInDat)," item(s) missed in data frame will be removed from Q matrix: \n    ",paste(notInDat,collapse=", "),"\n",sep=""))
        qMatrix <- checkQmatrixConsistency(qMatrix[-match(notInDat, qMatrix[,1]),])
      }
      if(length(notInQ)>0) {
        cat(paste("Following ", length(notInQ)," item(s) missed in Q matrix will be removed from data: \n    ",paste(notInQ,collapse=", "),"\n",sep=""))
      }
      all.Names[["variablen"]] <- qMatrix[,1]  } ;   flush.console()
    cic <- checkItemConsistency(dat=dat, allNam = all.Names, remove.missing.items=remove.missing.items, verbose=verbose, removeMinNperItem=removeMinNperItem, minNperItem=minNperItem, remove.constant.items=remove.constant.items, model.statement=model.statement)
    cbc <- checkBGV(allNam = cic[["allNam"]], dat=cic[["dat"]], software=software, remove.no.answersHG=remove.no.answersHG, remove.vars.DIF.missing=remove.vars.DIF.missing, namen.items.weg=cic[["namen.items.weg"]], remove.vars.DIF.constant=remove.vars.DIF.constant)
    if(length(cbc[["namen.items.weg"]])>0)  {
      cat(paste("Remove ",length(unique(cbc[["namen.items.weg"]]))," test item(s) overall.\n",sep=""))
      cbc[["allNam"]]$variablen <- setdiff(cbc[["allNam"]]$variablen, unique(cbc[["namen.items.weg"]]) )
      qMatrix             <- qMatrix[match(cbc[["allNam"]]$variablen, qMatrix[,1]),]
    }
    pwvv<- personWithoutValidValues(dat=cbc[["dat"]], allNam=cbc[["allNam"]], remove.no.answers=remove.no.answers)
    cpsc<- checkPersonSumScores(datL = pwvv[["datL"]], allNam = cbc[["allNam"]], dat=pwvv[["dat"]], remove.failures=remove.failures)
    if(check.for.linking == TRUE) {
      linkNaKeep <- checkLink(dataFrame = cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], remove.non.responser = FALSE, verbose = FALSE )
      linkNaOmit <- checkLink(dataFrame = cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], remove.non.responser = TRUE, verbose = FALSE )
      if(linkNaKeep == FALSE & linkNaOmit == TRUE )  {cat("Note: Dataset is not completely linked. This is probably only due to missings on all cases.\n")}
      if(linkNaKeep == TRUE )                        {cat("Dataset is completely linked.\n")}
    }
    met <- adaptMethod(method=method, software=software, nodes=nodes)
    if(length(cbc[["namen.all.hg"]])>0) {all.hg.char <- sapply(cbc[["namen.all.hg"]], FUN=function(ii) {max(nchar(as.character(na.omit(cpsc[["dat"]][,ii]))))})} else {all.hg.char <- NULL}
    if ( software == "conquest" )   {
      var.char  <- sapply(cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], FUN=function(ii) {max(nchar(as.character(na.omit(ii))))})
      no.number <- setdiff(1:length(var.char), grep("[[:digit:]]",var.char))
      if(length(no.number)>0) {var.char[no.number] <- 1}
      if(use.letters == TRUE)   {
        rec.statement <- paste(0:25,"='",LETTERS,"'",sep="",collapse="; ")
        daten.temp <- cpsc[["dat"]]
        for (i in cbc[["allNam"]][["variablen"]])  {
          cpsc[["dat"]][,i] <- car::recode(cpsc[["dat"]][,i], rec.statement)}
        var.char <- rep(1,length(cbc[["allNam"]][["variablen"]]))}
    }
    daten   <- data.frame(ID=as.character(cpsc[["dat"]][,cbc[["allNam"]][["ID"]]]), cpsc[["dat"]][,cbc[["namen.all.hg"]], drop = FALSE], cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], stringsAsFactors = FALSE)
    deskRes <- desk.irt(daten = daten, itemspalten = match(cbc[["allNam"]][["variablen"]], colnames(daten)), percent = TRUE)
    crit    <- which (deskRes[,"valid"] < minNperItem)
    if ( length(crit)>0) {
      cat ( paste ( "Following ",length(crit), " items with less than ",minNperItem," item responses:\n",sep=""))
      options(width=1000)
      print(deskRes[crit,-match(c("item.nr", "Label", "KB", "Codes", "Abs.Freq", "Rel.Freq"), colnames(deskRes))], digits = 3)
    }
    if(inherits(try(discrim <- item.diskrim(daten,match(cbc[["allNam"]][["variablen"]], colnames(daten)))  ),"try-error"))  {
      discrim <- item.diskrim(daten.temp,match(cbc[["allNam"]][["variablen"]], colnames(daten.temp)))
      rm(daten.temp)
    }
    if ( length ( cbc[["allNam"]][["schooltype.var"]] ) > 0 ) {
      deskS <- by ( data = cpsc[["dat"]], INDICES = cpsc[["dat"]][, cbc[["allNam"]][["schooltype.var"]] ], FUN = function ( st ) {
        drst <- desk.irt(daten = st, itemspalten = match(cbc[["allNam"]][["variablen"]], colnames(st)), percent = TRUE)
        colnames(drst) <- car::recode (colnames(drst) , paste0("'item.p'='item.p.",st[1,cbc[["allNam"]][["schooltype.var"]]],"'") )
        return(drst)})
      for ( uu in 1:length( deskS) ) {
        matchU <- match(c("item.nr","Label", "KB", "cases", "Missing", "valid", "Codes" , "Abs.Freq", "Rel.Freq"), colnames(deskS[[uu]]))
        stopifnot ( length (which(is.na(matchU))) == 0 , ncol(deskS[[uu]]) - length ( matchU) == 2)
        deskRes <- merge ( deskRes, deskS[[uu]][,-matchU], by = "item.name", all = TRUE)
      }
    }
    lab <- data.frame(itemNr = 1:length(cbc[["allNam"]][["variablen"]]), item = cbc[["allNam"]][["variablen"]], stringsAsFactors = FALSE)
    if(!is.null(anchor))  {
      ankFrame <- anker (lab = lab, prm = anchor, qMatrix = qMatrix, domainCol=domainCol, itemCol=itemCol, valueCol=valueCol, multicore = ismc)
    } else {
      ankFrame <- NULL
      if ( fitTamMmlForBayesian == FALSE ) {
        cat("   Note: 'anchor' is necessary if 'fitTamMmlForBayesian' is FALSE. Because 'anchor' is NULL, 'fitTamMmlForBayesian' is set to be TRUE now.\n")
        fitTamMmlForBayesian <- TRUE
      }
    }
    if ( software == "conquest" )   {
      daten$ID <- gsub ( " ", "0", formatC(daten$ID, width=max(as.numeric(names(table(nchar(daten$ID)))))) )
      fixed.width <- c(as.numeric(names(table(nchar(daten[,"ID"])))), all.hg.char, rep(max(var.char),length(var.char)))
      txt  <-  capture.output ( gdata::write.fwf(daten , colnames = FALSE,rownames = FALSE, sep="",quote = FALSE,na=".", width=fixed.width))
      stopifnot(length(table(nchar(txt)))==1)
      rm(txt)
      gdata::write.fwf(daten , file.path(dir,paste(analysis.name,".dat",sep="")), colnames = FALSE,rownames = FALSE, sep="",quote = FALSE,na=".", width=fixed.width)
      colnames(lab) <- c("===>","item")
      write.table(lab,file.path(dir,paste(analysis.name,".lab",sep="")),col.names = TRUE,row.names = FALSE, dec = ",", sep = " ", quote = FALSE)
      batch <- paste( normalize.path(conquest.folder),paste(analysis.name,".cqc",sep=""), sep=" ")
      write(batch, file.path(dir,paste(analysis.name,".bat",sep="")))
      foo <- gen.syntax(Name=analysis.name, daten=daten, all.Names = cbc[["allNam"]], namen.all.hg = cbc[["namen.all.hg"]], all.hg.char = all.hg.char, var.char= max(var.char), model=qMatrix, anchored=anchor, pfad=dir, n.plausible=n.plausible, compute.fit = compute.fit,
                        constraints=constraints, std.err=std.err, distribution=distribution, method=met[["method"]], n.iterations=n.iterations, nodes=met[["nodes"]], p.nodes=p.nodes, f.nodes=f.nodes, converge=converge,deviancechange=deviancechange, equivalence.table=equivalence.table, use.letters=use.letters, model.statement=model.statement, conquest.folder = conquest.folder, allowAllScoresEverywhere = allowAllScoresEverywhere, seed = seed, export = export)
      if(!is.null(anchor))  {
        write.table(ankFrame[["resConquest"]], file.path(dir,paste(analysis.name,".ank",sep="")) ,sep=" ", col.names = FALSE, row.names = FALSE, quote = FALSE)
      }
      if(file.exists( file.path ( dir,  paste(analysis.name,".log",sep=""))) )  {
        cat(paste("Found existing log file '",paste(analysis.name,".log",sep=""), "' in folder '",dir,"'\nConquest analysis will overwrite log file. Original log file will be saved as '",paste(analysis.name,"_old.log'\n",sep=""),sep=""))
        do <- file.rename(from = file.path(dir, paste(analysis.name,".log",sep="")), to = file.path(dir, paste(analysis.name,"_old.log",sep="")))
      }
      options(scipen = original.options); flush.console()
      ret <- list ( software = software, input = paste("\"", file.path(dir, paste(analysis.name,"cqc",sep=".")), "\"", sep=""), conquest.folder = paste("\"", conquest.folder, "\"", sep=""), dir=dir, analysis.name=analysis.name, model.name = analysis.name, qMatrix=qMatrix, all.Names=cbc[["allNam"]], deskRes = deskRes, discrim = discrim, perNA=pwvv[["perNA"]], per0=cpsc[["per0"]], perA = cpsc[["perA"]], perExHG = cbc[["perExHG"]], itemsExcluded = cbc[["namen.items.weg"]], daten=daten, method=met[["method"]], nodes=met[["nodes"]], p.nodes=p.nodes, f.nodes=f.nodes)
      class(ret) <-  c("defineConquest", "list")
      return ( ret )  }
    if ( software == "tam" )   {
      cat(paste("Q matrix specifies ",ncol(qMatrix)-1," dimension(s).\n",sep=""))
      anchor          <- prepAnchorTAM(ank = ankFrame[["resTam"]], allNam = cbc[["allNam"]])
      est.slopegroups <- prepEstSlopegroupsTAM(esg = est.slopegroups, allNam = cbc[["allNam"]])
      fixSlopeMat     <- prepFixSlopeMatTAM(fsm = fixSlopeMat, allNam = cbc[["allNam"]], qma =  qMatrix, slopeMatDomainCol=slopeMatDomainCol, slopeMatItemCol=slopeMatItemCol, slopeMatValueCol=slopeMatValueCol, dat=daten, irtmodel=irtmodel)
      guessMat        <- prepGuessMat(guessMat, allNam = fixSlopeMat[["allNam"]])
      control <- list ( snodes = met[["snodes"]] , QMC=met[["QMC"]], convD = deviancechange ,conv = converge , convM = .0001 , Msteps = Msteps , maxiter = n.iterations, max.increment = 1 ,
                        min.variance = .001 , progress = progress , ridge=0 , seed = seed , xsi.start0=FALSE,  increment.factor=increment.factor , fac.oldxsi= fac.oldxsi)
      if ( !is.null(met[["nodes"]])) { control$nodes <- met[["nodes"]] }
      ret     <- list ( software = software, constraint = match.arg(constraints) , qMatrix=qMatrix, anchor=anchor,  all.Names=fixSlopeMat[["allNam"]], daten=daten, irtmodel=fixSlopeMat[["irtmodel"]], est.slopegroups = est.slopegroups, guessMat=guessMat, control = control, n.plausible=n.plausible, dir = dir, analysis.name=analysis.name, deskRes = deskRes, discrim = discrim, perNA=pwvv[["perNA"]], per0=cpsc[["per0"]], perA = cpsc[["perA"]], perExHG = cbc[["perExHG"]], itemsExcluded = cbc[["namen.items.weg"]], fixSlopeMat = fixSlopeMat[["slopMat"]], estVar = fixSlopeMat[["estVar"]], pvMethod = pvMethod,  fitTamMmlForBayesian=fitTamMmlForBayesian)
      class(ret) <-  c("defineTam", "list")
      return ( ret )    }   }
}
