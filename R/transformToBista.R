transformToBista <- function(equatingList, refPop, cuts, weights = NULL,
                             defaultM = 500, defaultSD = 100, q3bound = .20,
                             roman = FALSE, vera = TRUE, idVarName = NULL,
                             years = NULL){
### checks ---------------------------------------------------------------------

  weights <- eatTools::makeDataFrame(weights)
  lapply(c(defaultM, defaultSD, q3bound), checkmate::assert_numeric, len = 1, lower = 0)
  checkmate::assert_numeric(years, len = 2, null.ok = TRUE)
  checkmate::assert_character(idVarName, len = 1, null.ok = TRUE)

### function -------------------------------------------------------------------

  mr  <- FALSE
  if(missing(refPop)) {
    mr  <- TRUE
    cat("'refPop' was not defined. Treat current sample as drawn from the reference population.\n")
    flush.console()
  }  else  {
    refPop <- eatTools::makeDataFrame(refPop)
    for ( i in 2:ncol(refPop)) {
      if(!inherits(refPop[,i], c("integer", "numeric"))) {stop("All columns of 'refPop' except for the first one must be numeric.")}
    }
  }
  if( missing(cuts)) { cutsMis <- TRUE }  else  { cutsMis <- FALSE }
  nam1<- names(equatingList[["items"]])                                    ### hier stehen in der regel die beiden Modellnamen, also quasi names(table(equatingList[["results"]][,"model"]))
  isRunM<- all(c("model" , "source" , "var1" , "var2" , "type" , "indicator.group",
                 "group", "par", "derived.par", "value") %in% names(equatingList[["results"]]))
  if (isRunM) {
    it     <- itemFromRes(equatingList[["results"]])
  }  else  {
    it     <- equatingList[["results"]]
  }
  dims   <- unique(it[,"dimension"])
  if (is.null(dims)) {
    dims <- unique(na.omit(equatingList[["results"]][,"group"]))
    warning(paste0("Cannot extract dimensions from 'results' object. This should only occur for bayesian plausible values imputation. Assume following dimensions: \n    '",
                   paste(dims, collapse = "', '"),"'."))
    if(vera==TRUE) {
      warning("'vera' must be FALSE for bayesian plausible values imputation. Set 'vera' to FALSE.")
      vera <- FALSE
    }
  }
  if (isRunM) {
    id     <- unique(equatingList[["results"]][intersect(which(equatingList[["results"]][,"type"] == "tech"),
                                                         which(equatingList[["results"]][,"par"] == "ID")),"derived.par"])
    if(length(id)!=1) {
      id   <- getIdVarName(id=NULL, idVarName, verbose=TRUE)
    }
    refList<- lapply ( dims, FUN = function (dimname) {
      rex  <- pvFromRes(equatingList[["results"]][unique(c(which(equatingList[["results"]][,"group"] == dimname),
                                                           which(equatingList[["results"]][,"type"] == "tech"))), ],
                        toWideFormat = FALSE, idVarName=idVarName, verbose=FALSE)
      if (is.null(rex)) {return(NULL)}
      if ( is.null(weights) ) {
        txt <- capture.output ( msd <- eatRep::repMean ( datL = rex, ID = id, imp = "imp",
                                                         dependent = "value", na.rm = TRUE))
        msd <- adaptEatRepVersion(msd)
      }  else  {
        rex <- eatTools::mergeAttr ( rex, weights , by.x = id, by.y = colnames(weights)[1],
                                     all.x = TRUE, all.y = FALSE,  setAttr = FALSE,
                                     unitName = "cases", xName = paste0("plausible values for dimension ",dimname),
                                     yName = "weights", verbose = c("match", "dataframe"))
        mis <- which(is.na(rex[,colnames(weights)[2]]))
        if ( length(mis) > 0 ) {
          if(length(mis) == nrow(rex)) {stop(paste("Mergin of weights and plausible values for '",
                                                   dimname, "' failed. No common units."))}
          cat(paste ( "Found ",length(mis)," missing values in the 'weights' frame.\n    Cases with missing values on weighting variable will be ignored for transformation.\n",
                      sep=""))
          rex <- rex[-mis,]
        }
        txt <- capture.output(msd <- eatRep::repMean(datL = rex, ID = id, imp = "imp",
                                                     wgt = colnames(weights)[2],
                                                     dependent = "value", na.rm = TRUE))
        msd <- adaptEatRepVersion(msd)
      }
      rp <- data.frame(domain = dimname, m = msd[intersect(which(msd[,"parameter"] == "mean"),
                                                           which(msd[,"coefficient"] == "est")),"value"],
                       sd = msd[intersect(which(msd[,"parameter"] == "sd"), which(msd[,"coefficient"] == "est")),"value"])
      return(list (msd = msd , rp=rp))})
    names(refList) <- dims
    ref    <- do.call("rbind", lapply(refList, FUN = function ( u ) { u[["rp"]] }))
    if ( isTRUE(mr) ) {
      refPop <- ref
    }   else  {
      mis <- which(is.na(refPop))
      if ( length(mis) >0) {
        stopifnot ( nrow(ref ) == nrow(refPop))
        mat <- merge( 1:nrow(refPop), 1:ncol(refPop), by = NULL)
        refPop[unique(mat[mis,"x"]), unique(mat[mis,"y"])] <- ref[unique(mat[mis,"x"]),
                                                                  unique(mat[mis,"y"])]
      }
    }
  }
  if(ncol ( refPop ) == 3) {
    cat ( paste("The 'refPop' data.frame does not include information about reference population mean/SD on Bista metric. Values will be defaulted to ",
                defaultM,"/",defaultSD,".\n",sep=""))
    refPop[,4] <- defaultM; refPop[,5] <- defaultSD
  }  else  {
    if ( ncol ( refPop) != 5 ) { stop ( "Invalid 'refPop'.\n") }
  }
  if(!all(dims %in% refPop[,1]))  {warning(paste0("Dimension names in the 'results' object '",
                                                  paste(dims, collapse="', '"),
                                                  "' do not match to names in the first columns of 'refPop': '",
                                                  paste(refPop[,1], collapse="', '"),"."))}
  if(!all(dims %in% names(cuts))) {warning(paste0("Dimension names in the 'results' object '",
                                                  paste(dims, collapse="', '"),
                                                  "' do not match to names in the cut scores object: '",
                                                  paste(names(cuts), collapse="', '"),"."))}
  modN<- lapply(nam1, FUN = function ( mod ) {
    nam2 <- names(equatingList[["items"]][[mod]])
    dimN <- lapply(nam2, FUN = function ( dims ) {
      if (isRunM) {
        resMD<- equatingList[["results"]][unique(c(intersect(which(equatingList[["results"]][,"model"] == mod),
                                                             which(equatingList[["results"]][,"group"] == dims)),
                                                   which(equatingList[["results"]][,"type"] == "tech"))),]
        rex  <- pvFromRes(resMD, toWideFormat = TRUE, idVarName = idVarName, verbose=FALSE)
        if (!is.null(rex)) {
          if ( length ( rex[,id]) != unique(length ( rex[,id])) ) {
            stop(paste( "Model '",mod,"', Dimension '",dims,"': cases according to '",
                        id,"' variable are not unique.\n",sep=""))
          }
        }
        offSet  <- grep("offset", as.character(resMD[,"par"]))
        if(length(offSet)>0) {  resMD[,"par"] <- car::recode ( resMD[,"par"], "'offset'='est'") }
        itFrame <- itemFromRes(resMD)
      }  else  {
        itFrame <- it[intersect(which(it[,"model"] == mod), which(it[,"dimension"] == dims)),]
      }
      if ( !is.null(itFrame) && !itFrame[1,"dimension"] %in% refPop[,1] ) {
        cat(paste("Cannot found dimension '",itFrame[1,"dimension"],"' in the first column of the 'refPop' argument. Skip transformation ... \n",
                  sep=""))
        return ( list ( itempars = NULL, personpars = NULL, rp = NULL))
      }  else  {
        if ( is.null ( itFrame ) ) {
          cat(paste0("Model '",mod,"', dimension '",dims,"': No item parameters found. This should only occur for bayesian plausible values imputation. Transformation of item parameters will be skipped.\n"))
        }  else  {
          if ( isFALSE(cutsMis) ) {
            if ( !itFrame[1,"dimension"] %in% names(cuts) ) {
              stop(paste("Cannot found dimension '",itFrame[1,"dimension"],
                         "' in the 'cuts' list.",sep=""))}
            mat1<- match( itFrame[1,"dimension"], names(cuts))
            if ( !"values" %in% names(cuts[[mat1]]) ) {
              stop(paste("'cuts' must be a named list. Cannot found 'values' element for dimension '",
                         itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
            if ( length(cuts[[mat1]])>1) {
              if ( !"labels" %in% names(cuts[[mat1]]) ) {
                stop(paste("'cuts' must be a named list. Cannot found 'labels' element for dimension '",
                           itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
            }
          }
          mat <- match( itFrame[1,"dimension"], refPop[,1])
          stopifnot(length(mat)==1)
          itFrame[,"estTransf"] <- itFrame[,"est"] + equatingList[["items"]][[mod]][[dims]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dims]][["method"]] ]]
          if ( isTRUE(mr) ) {
            if ( equatingList[["items"]][[mod]][[dims]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dims]][["method"]] ]] != 0) {
              cat("W A R N I N G: Preceding Equating without 'refPop' definition. Sure you want to use current sample as drawn from the reference population?\n")
              refPop[mat,2] <- refPop[mat,2]+ equatingList[["items"]][[mod]][[dims]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dims]][["method"]] ]]
            }
          }
          itFrame[,"estTransf625"]   <- itFrame[,"estTransf"] + log(0.625/(1-0.625))
          itFrame[,"estTransfBista"] <- (itFrame[,"estTransf625"] - refPop[mat,2]) / refPop[mat,3] * refPop[mat,5] + refPop[mat,4]
          if ( isFALSE(cutsMis) ) {
            traitLevel            <- eatTools::num.to.cat(x = itFrame[,"estTransfBista"], cut.points = cuts[[mat1]][["values"]], cat.values = cuts[[mat1]][["labels"]])
            itFrame[,"traitLevel"]<- traitLevel
          }
          itFrame[,"linkingConstant"]<- equatingList[["items"]][[mod]][[dims]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dims]][["method"]] ]]
          itFrame[,"linkingMethod"]  <- equatingList[["items"]][[mod]][[dims]][["method"]]
          itFrame[,"nLinkitems"]     <- equatingList[["items"]][[mod]][[dims]][["eq"]][["descriptives"]][["N.Items"]]
          itFrame[,"linkingError"]   <- equatingList[["items"]][[mod]][[dims]][["eq"]][["descriptives"]][["linkerror"]]
          itFrame[,"linkingErrorTransfBista"] <- ( (itFrame[,"linkingError"]^2) * (refPop[mat,5]^2) / (refPop[mat,3]^2) )^0.5
        }
        if (isRunM) {
          pv  <- pvFromRes(resMD, toWideFormat = FALSE, idVarName=idVarName, verbose=FALSE)
          equ <- equatingList[["items"]][[mod]][[dims]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dims]][["method"]] ]]
          if (!exists("mat")) { mat <- match(dims,  refPop[,1]) }
          pv[,"valueTransfBista"] <- (pv[,"value"] + equ - refPop[mat,2]) / refPop[mat,3] * refPop[mat,5] + refPop[mat,4]
          if (!is.null(pv)) {
            if ( is.null(weights) ) {
              txt <- capture.output(msdF <- eatRep::repMean(datL = pv, ID = id,
                                                            imp = "imp",
                                                            dependent = "valueTransfBista",
                                                            na.rm = TRUE))
              msdF<- adaptEatRepVersion(msdF)
            }  else  {
              pvF <- eatTools::mergeAttr(pv, weights , by.x = id, by.y = colnames(weights)[1],
                                         all.x = TRUE, all.y = FALSE,  setAttr = FALSE,
                                         unitName = "cases", xName = paste0("plausible values for dimension ",dims),
                                         yName = "weights", verbose = c("match", "dataframe"))
              mis <- which(is.na(pvF[,colnames(weights)[2]]))
              if ( length(mis) > 0 ) {
                cat(paste ( "Found ",length(mis)," missing values in the 'weights' frame.\n    Cases with missing values on weighting variable will be ignored for transformation.\n",
                            sep=""))
                pvF <- pvF[-mis,]
              }
              txt <- capture.output(msdF <- eatRep::repMean(datL = pvF, ID = id,
                                                            imp = "imp", wgt = colnames(weights)[2],
                                                            dependent = "valueTransfBista", na.rm = TRUE))
              msdF<- adaptEatRepVersion(msdF)
            }
            msdFok <- c(msdF[intersect(which(msdF[,"parameter"] == "mean"),
                                       which(msdF[,"coefficient"] == "est")),"value"],
                        msdF[intersect(which(msdF[,"parameter"] == "sd"),
                                       which(msdF[,"coefficient"] == "est")),"value"])
          }  else  {
            message("Results object does not contain any plausible values. Skip transformation of linking error for competence levels.")
          }
          if ( !is.null(pv) && !is.null ( itFrame )) {
            if ( cutsMis == FALSE & !is.null ( equatingList[["items"]] )) {
              cts <- c( -10^6, cuts[[mat1]][["values"]], 10^6)
              le  <- do.call("rbind", lapply ( (length(cts)-1):1 , FUN = function ( l ) {
                kmp<- c(cts[l], cts[l+1])
                a1 <- sum ( dnorm ( ( kmp - refPop[mat,4]) / refPop[mat,5] ) * c(-1,1) / refPop[mat,5] )
                a2 <- sum ( dnorm ( ( kmp - msdFok[1]) / msdFok[2] ) * c(-1,1) / msdFok[2] )
                if(a2 == 0 ) {cat("mutmasslicher fehler.\n")}
                del<- ( (  a1^2 + a2^2 ) * (unique(itFrame[,"linkingErrorTransfBista"])^2) / 2  )^0.5
                del<- data.frame ( traitLevel = attr(traitLevel, "cat.values")[l],
                                   linkingErrorTraitLevel = del )
                return(del)}))
              ori <- colnames(itFrame)
              chk <- unique(le[,"traitLevel"]) %in% unique(itFrame[,"traitLevel"])
              if ( length( which(chk == FALSE)) > 0) {
                warning(paste("Model '",unique(itFrame[,"model"]),"', dimension '",
                              unique(itFrame[,"dimension"]),"': No items on trait level(s) '",
                              paste(unique(le[,"traitLevel"])[which(chk == FALSE)],
                                    collapse = "', '"), "'.", sep=""))
              }
              itFrame <- eatTools::mergeAttr(itFrame, le, by = "traitLevel", sort = FALSE,
                                             all.x = TRUE, all.y = FALSE,  setAttr = FALSE,
                                             unitName = "trait levels", xName = "item parameter list",
                                             yName = "linking error list", verbose = c("match"))
              itFrame <- itFrame[,c(ori, "linkingErrorTraitLevel")]
            }
            itFrame[,"refMean"]        <- refPop[mat,2]
            itFrame["refSD"]           <- refPop[mat,3]
            itFrame[,"refTransfMean"]  <- refPop[mat,4]
            itFrame[,"refTransfSD"]    <- refPop[mat,5]
          }
          if(!exists("mat1") ) {mat1 <- match(dims, names(cuts)); stopifnot(length(mat1)==1)}
          if (!is.null(pv)) {
            if ( isFALSE(cutsMis) ) { pv[,"traitLevel"] <- eatTools::num.to.cat(x = pv[,"valueTransfBista"],
                                                                                cut.points = cuts[[mat1]][["values"]],
                                                                                cat.values = cuts[[mat1]][["labels"]])}
            pv[,"dimension"]  <- pv[,"group"]
            if(!exists("le")) {
              warning("Skip check whether all competence levels are occupied (due to bayesian plausible values imputation).")
            }  else  {
              chk <- unique(le[,"traitLevel"]) %in% unique(pv[,"traitLevel"])
              if ( length( which(chk == FALSE)) > 0) {
                warning(paste("Model '",unique(itFrame[,"model"]),"', dimension '",
                              unique(itFrame[,"dimension"]),"': No plausible values on trait level(s) '",
                              paste( unique(le[,"traitLevel"])[which(chk == FALSE)],
                                     collapse = "', '"), "'.", sep=""))
              }
              stopifnot ( length( unique ( na.omit(itFrame[,"linkingErrorTransfBista"]))) %in% 0:1)
              pv[,"linkingError"] <- equatingList[["items"]][[mod]][[dims]][["eq"]][["descriptives"]][["linkerror"]]
              pv[,"linkingErrorTransfBista"] <- unique ( itFrame[,"linkingErrorTransfBista"])
            }
            ori <- colnames(pv)
            if ( cutsMis == FALSE && !is.null ( equatingList[["items"]]) && exists("le") ) {
              pv  <- eatTools::mergeAttr(pv, le, by = "traitLevel", sort = FALSE,
                                         all.x = TRUE, all.y = FALSE, setAttr = FALSE,
                                         unitName = "trait levels", xName = "plausible values",
                                         yName = "linking error list", verbose = c("match"))
              pv  <- pv[,c(ori, "linkingErrorTraitLevel")]
            }
            if (!is.null(weights)) {
              pv  <- merge ( pv, weights , by.x = id, by.y = colnames(weights)[1],
                             all.x = TRUE, all.y = FALSE)
            }
          }  else  {
            msdFok <- c(NA, NA)
          }
        }
        rp  <- refPop[mat,]
        colnames(rp) <- c("domain", "refMean", "refSD", "bistaMean", "bistaSD")
        if (isRunM) {
          rp  <- cbind ( model = mod, rp, focusMean = msdFok[1], focusSD = msdFok[2])
        }  else  {
          pv <- NULL
        }
        return(list ( itempars = itFrame, personpars = pv, rp = rp))
      }  })
    itempars<- do.call(plyr::rbind.fill, lapply ( dimN, FUN = function ( x ) { x[["itempars"]]}))
    perspar <- do.call("rbind", lapply ( dimN, FUN = function ( x ) { x[["personpars"]]}))
    rp      <- do.call("rbind", lapply ( dimN, FUN = function ( x ) { x[["rp"]]}))
    return( list ( itempars = itempars, personpars = perspar, rp=rp)) } )
  personpars <- do.call("rbind", lapply ( modN, FUN = function ( x ) { x[["personpars"]]}))
  itempars   <- do.call(plyr::rbind.fill, lapply ( modN, FUN = function ( x ) { x[["itempars"]]}))
  rp         <- do.call("rbind", lapply ( modN, FUN = function ( x ) { x[["rp"]]}))
  if ( isFALSE(vera) || isFALSE(isRunM) ) {
    itemVera <- NULL
  }  else  {
    itemVera <- createItemVeraObj(itempars=itempars, roman=roman,
                                  results = equatingList[["results"]], q3bound=q3bound)
  }
  if (!is.null(years)) {
    stopifnot(length(years) == 2 && length(unique(years)) == 2)
    leo  <- createLinkingErrorObject(itempars=itempars, years=years)
  }  else  {
    leo  <- NULL
  }
  if (isRunM) {
    context    <- equatingList[["results"]][which(equatingList[["results"]][,"type"]=="tech"),]
  }  else  {
    context    <- NULL
  }
  ret        <- list(itempars = itempars, personpars = personpars, refPop = refPop,
                     means = rp, all.Names = context, itemparsVera = itemVera, linkingErrors = leo)
  class(ret) <- c("list", "transfBista")
  return( ret ) }
