### runModel is usually used after defineModel().

runModel <- function(defineModelObj, show.output.on.console = FALSE,
                     show.dos.console = TRUE, wait = TRUE) {

### checks ---------------------------------------------------------------------

  checkmate::assert_list(defineModelObj)
  lapply(c(show.output.on.console, show.dos.console, wait), checkmate::assert_logical, len = 1)

### function -------------------------------------------------------------------

  if (inherits(defineModelObj, "defineMultiple") ) {
    if(is.null ( attr(defineModelObj, "split")[["nCores"]] ) || attr(defineModelObj, "split")[["nCores"]] == 1 ) {
      res <- lapply(defineModelObj, FUN = function ( r ) {
        ret <- runModel ( defineModelObj = r, show.output.on.console = show.output.on.console, show.dos.console = show.dos.console, wait = wait)
        return(ret)})
    }  else  {
      doIt<- function (laufnummer,  ... ) {
        if(!exists("runModel"))  { library(eatModel) }
        ret <- runModel ( defineModelObj = defineModelObj[[laufnummer]], show.output.on.console = show.output.on.console, show.dos.console = show.dos.console, wait = TRUE)
        return(ret) }
      beg <- Sys.time()
      if ( attr(defineModelObj, "split")[["mcPackage"]] == "parallel") {
        cl  <- makeCluster(attr(defineModelObj, "split")[["nCores"]], type = "SOCK")
      }  else  {
        cl  <- future::makeClusterPSOCK(attr(defineModelObj, "split")[["nCores"]], verbose=FALSE)
      }
      res <- clusterApply(cl = cl, x = 1:length(defineModelObj), fun = doIt , show.output.on.console = show.output.on.console, show.dos.console = show.dos.console, wait = wait)
      stopCluster(cl)
      cat(paste ( length(defineModelObj), " analyses finished: ", sep="")); print( Sys.time() - beg, digits=3)
    }
    class(res) <- c("runMultiple", "list")
    attr(res, "split") <- attr(defineModelObj, "split")
    return(res)
  } else {
    if(inherits(defineModelObj, "defineConquest")) {
      oldPfad <- getwd()
      setwd(defineModelObj$dir)
      suppressWarnings(system(paste(defineModelObj$conquest.folder," ",defineModelObj$input,sep=""),invisible=!show.dos.console,show.output.on.console=show.output.on.console, wait=wait) )
      if(wait == FALSE) { Sys.sleep(0.2) }
      setwd(oldPfad)
      class(defineModelObj) <- c("runConquest", "list")
      return ( defineModelObj )
    }
    if(inherits(defineModelObj, "defineTam")) {
      if ( show.output.on.console == TRUE ) { defineModelObj[["control"]][["progress"]] <- TRUE }
      if(length( defineModelObj[["all.Names"]][["HG.var"]])>0)     { Y <- defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["HG.var"]], drop=FALSE] } else { Y <- NULL }
      if(length( defineModelObj[["all.Names"]][["weight.var"]])>0) { wgt <- as.vector(defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["weight.var"]]])} else {wgt <- NULL}
      if(length( defineModelObj[["all.Names"]][["group.var"]])>0)  { group <- as.vector(defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["group.var"]]])} else {group <- NULL}
      stopifnot(all(defineModelObj[["qMatrix"]][,1] == defineModelObj[["all.Names"]][["variablen"]]))
      if(length(defineModelObj[["all.Names"]][["DIF.var"]]) == 0 ) {
        if( defineModelObj[["irtmodel"]] %in% c("1PL", "PCM", "PCM2", "RSM")) {
          if ( isTRUE(defineModelObj[["fitTamMmlForBayesian"]]) ) {
            mod  <- tam.mml(resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], constraint = defineModelObj[["constraint"]], pid = defineModelObj[["daten"]][,"ID"], Y = Y, Q = defineModelObj[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = defineModelObj[["anchor"]], irtmodel = defineModelObj[["irtmodel"]], pweights = wgt, control = defineModelObj[["control"]], group=group)
          }  else  {
            mod  <- tamObjForBayesianPV (anchor = defineModelObj[["anchor"]], qMatrix = defineModelObj[["qMatrix"]], resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y=Y)
          }
        }
        if( defineModelObj[["irtmodel"]] %in% c("2PL", "GPCM", "2PL.groups", "GPCM.design", "3PL") )  {
          if( defineModelObj[["irtmodel"]] == "3PL") {
            if ( isTRUE(defineModelObj[["fitTamMmlForBayesian"]]) ) {
              mod  <- tam.mml.3pl(resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y = Y, Q = defineModelObj[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = defineModelObj[["anchor"]], pweights = wgt, est.guess =defineModelObj[["guessMat"]],  est.variance = defineModelObj[["estVar"]], control = defineModelObj[["control"]], group=group)
            }  else  {
              mod  <- tamObjForBayesianPV (anchor = defineModelObj[["anchor"]], qMatrix = defineModelObj[["qMatrix"]], resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y=Y, slopeMatrix = defineModelObj[["fixSlopeMat"]])
            }
          }  else {
            if ( defineModelObj[["fitTamMmlForBayesian"]] == TRUE ) {
              mod  <- tam.mml.2pl(resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y = Y, Q = defineModelObj[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = defineModelObj[["anchor"]], irtmodel = defineModelObj[["irtmodel"]], est.slopegroups=defineModelObj[["est.slopegroups"]],pweights = wgt, B.fixed = defineModelObj[["fixSlopeMat"]], est.variance = defineModelObj[["estVar"]], control = defineModelObj[["control"]], group=group)
            }  else  {
              mod  <- tamObjForBayesianPV (anchor = defineModelObj[["anchor"]], qMatrix = defineModelObj[["qMatrix"]], resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y=Y, slopeMatrix = defineModelObj[["fixSlopeMat"]])
            }
          }
        }
      } else {
        assign(paste("DIF_",defineModelObj[["all.Names"]][["DIF.var"]],sep="") , as.data.frame (defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["DIF.var"]]]) )
        formel   <- as.formula(paste("~item - ",paste("DIF_",defineModelObj[["all.Names"]][["DIF.var"]],sep="")," + item * ",paste("DIF_",defineModelObj[["all.Names"]][["DIF.var"]],sep=""),sep=""))
        facetten <- as.data.frame (defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["DIF.var"]]])
        colnames(facetten) <- paste("DIF_",defineModelObj[["all.Names"]][["DIF.var"]],sep="")
        if ( isTRUE(defineModelObj[["fitTamMmlForBayesian"]]) ) {
          mod  <- tam.mml.mfr(resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], facets = facetten, constraint = defineModelObj[["constraint"]], formulaA = formel, pid = defineModelObj[["daten"]][,"ID"], Y = Y, Q = defineModelObj[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = defineModelObj[["anchor"]], irtmodel = defineModelObj[["irtmodel"]], pweights = wgt, control = defineModelObj[["control"]], group=group)
        }  else  {
          mod  <- tamObjForBayesianPV (anchor = defineModelObj[["anchor"]], qMatrix = defineModelObj[["qMatrix"]], resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y=Y, slopeMatrix = defineModelObj[["fixSlopeMat"]])
        }
      }
      attr(mod, "defineModelObj") <- defineModelObj[-match("daten", names(defineModelObj))]
      attr(mod, "Y")              <- Y
      return(mod)  }  }
}

### runModel() specific help functions -----------------------------------------

qMatToB <- function(qma, slp) {
  zei <- match( qma[,"item"], slp[,1])
  for ( i in 1:length(zei) ) {
    ind <- which(qma[i,] ==1 )
    stopifnot(length(ind)==1, qma[i,"item"] == slp[zei[i],1])
    qma[i,ind] <- slp[zei[i],2] }
  return(qma)}

tamObjForBayesianPV <- function(anchor, qMatrix, slopeMatrix = NULL, resp, pid, Y) {
  warning("To date, bayesian plausible values imputation only works for binary between-item dimensionality models.")
  if ( !is.null(slopeMatrix)) {
    qMatrix <- qMatToB ( qma = qMatrix, slp = slopeMatrix)
  }
  xsi.obj<- as.matrix(data.frame ( V1 = 0, V2 = anchor[,"parameter"] * (-1)))
  B.obj  <- array(unlist(lapply(2:ncol(qMatrix),
                                FUN = function (col) {data.frame(Cat0 = 0, Cat1 = qMatrix[,col])})),
                  dim = c(nrow(qMatrix), 2, ncol(qMatrix)-1),
                  dimnames = list(qMatrix[,"item"], c("Cat0", "Cat1"), paste0("Dim0", 1:(ncol(qMatrix)-1)) ))
  tamObj <- list(AXsi = xsi.obj, B = B.obj, resp = resp, Y=Y, pid = pid)
  class(tamObj) <- c("list", "tamBayes")
  return(tamObj)}





