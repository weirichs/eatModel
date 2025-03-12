transformToBista <- function(equatingList, refPop, cuts, weights = NULL,
                             defaultM = 500, defaultSD = 100, q3bound = .20,
                             roman = FALSE, vera = TRUE, idVarName = NULL,
                             years = NULL){
### checks ---------------------------------------------------------------------

  checkmate::assert_data_frame(weights, ncols = 2, null.ok = TRUE)
  lapply(c(defaultM, defaultSD, q3bound), checkmate::assert_numeric, len = 1, lower = 0)
  checkmate::assert_numeric(years, len = 2, null.ok = TRUE)
  checkmate::assert_character(idVarName, len = 1, null.ok = TRUE)

### function -------------------------------------------------------------------
       mr    <- FALSE
       if(missing(refPop)) {mr <- TRUE}
       checkmate::assert_character(idVarName, null.ok = TRUE, len = 1)
    ### Wenn die Funktion ein Objekt kriegt, wo nur zwei Itemparameterlisten equated wurden, muss hier anders verfahren werden, dann gibt es bspw. keine Linkingfehler auf Kompetenzstufenmetrik. Das kriegt die Funktion hier erstmal raus. Wenn 'isRunM' gleich TRUE, dann default
       isRunM<- all(c("model" , "source" , "var1" , "var2" , "type" , "indicator.group", "group", "par", "derived.par", "value") %in% names(equatingList[["results"]]))
       if (isRunM) {
           it     <- itemFromRes(equatingList[["results"]])
           mods   <- names(equatingList[["items"]])
           id     <- unique(equatingList[["results"]][intersect(which(equatingList[["results"]][,"type"] == "tech"), which(equatingList[["results"]][,"par"] == "ID")),"derived.par"])
           if(length(id)!=1) { id   <- getIdVarName(id=NULL, idVarName, verbose=TRUE)}
       }  else  {
           it     <- equatingList[["results"]]
           mods   <- "model1"
           id     <- NULL
       }
       dims   <- unique(it[,"dimension"])
       if (is.null(dims)) {
           dims <- unique(na.omit(equatingList[["results"]][,"group"]))
           warning(paste0("Cannot extract dimensions from 'results' object. This should only occur for bayesian plausible values imputation. Assume following dimensions: \n    '",paste(dims, collapse = "', '"),"'."))
           if(vera==TRUE) {
              warning("'vera' must be FALSE for bayesian plausible values imputation. Set 'vera' to FALSE.")
              vera <- FALSE
           }
       }
       if( missing(cuts)) { cutsMis <- TRUE }  else  { cutsMis <- FALSE }
    ### wenn equatet wurde, sollte auch 'refPop' definiert sein (es sei denn, es wurde verankert skaliert) ... in dem Fall wird 'refPop' auf Konsistenz gecheckt
    ### wenn 'refPop' fehlt, wird es fuer alle gegebenen Dimensionen anhand der Gesamtstichprobe berechnet ... beides macht die folgende Funktion
       refPop <- generateOrCheckRefPop(equatingList = equatingList, refPop = refPop, mods=mods, dims=dims, isRunM = isRunM, id=id, weights=weights, defaultM=defaultM,defaultSD=defaultSD)
    ### fuer die Transformation selbst wird nach Modellen und Dimensionen getrennt
       modN<- lapply(mods, FUN = function ( mod ) {                             ### aeussere Schleife: geht ueber modelle
              dimN <- lapply(dims, FUN = function ( dimname ) {                 ### innere Schleife: geht ueber Dimensionen (innerhalb von modellen)
    ### check: sind Personen innerhalb jeder Dimension (und jeder Imputation) unique? ... 'redMD' ist ein reduziertes Results-Objekt: nur die interessierende Dimension des interessierenden Modells + saemtliche "tech"-Variablen
                      if (isRunM) {
                          wahl <- intersect(which(equatingList[["results"]][,"model"] == mod), which(equatingList[["results"]][,"group"] == dimname))
                          if(length(wahl)==0) {return(NULL)}
                          resMD<- equatingList[["results"]][unique(c(wahl,  which(equatingList[["results"]][,"type"] == "tech"))),]
                          rex  <- pvFromRes(resMD, toWideFormat = TRUE, idVarName = idVarName, verbose=FALSE)
                          if (!is.null(rex)) {
                              if ( length ( rex[,id]) != unique(length ( rex[,id])) ) {
                                   stop(paste( "Model '",mod,"', Dimension '",dimname,"': cases according to '", id,"' variable are not unique.\n",sep=""))
                              }
                          }
    ### check: keine verankerten parameter?
                          offSet  <- grep("offset", as.character(resMD[,"par"]))
                          if(length(offSet)>0) {  resMD[,"par"] <- car::recode ( resMD[,"par"], "'offset'='est'") }
                          itFrame <- itemFromRes(resMD)
                      }  else  {
                          itFrame <- it[intersect(which(it[,"model"] == mod), which(it[,"dimension"] == dimname)),]
                      }
                      if ( !is.null(itFrame) && !itFrame[1,"dimension"] %in% refPop[,"domain"] ) {
                            cat(paste("Cannot found dimension '",itFrame[1,"dimension"],"' in the first column of the 'refPop' argument. Skip transformation ... \n",sep=""))
                            return ( list ( itempars = NULL, personpars = NULL, rp = NULL))
                      }  else  {
                            if ( is.null ( itFrame ) ) {
                                cat(paste0("Model '",mod,"', dimension '",dimname,"': No item parameters found. This should only occur for bayesian plausible values imputation. Transformation of item parameters will be skipped.\n"))
                            }  else  {
    ### wenn cuts vom user definiert, wird hier auf plausibilitaet geprueft (kuenftig ggf. auslagern in separate funktion)
                                if ( isFALSE(cutsMis) ) {
                                     if ( !itFrame[1,"dimension"] %in% names(cuts) ) { stop(paste("Cannot found dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.",sep=""))}
                                     mat1<- match( unique(itFrame[,"dimension"]), names(cuts))
                                     if ( !"values" %in% names(cuts[[mat1]]) ) { stop(paste("'cuts' must be a named list. Cannot found 'values' element for dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
                                     if ( length(cuts[[mat1]])>1) {
                                         if ( !"labels" %in% names(cuts[[mat1]]) ) { stop(paste("'cuts' must be a named list. Cannot found 'labels' element for dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
                                     }
                                }
                                mat <- merge(unique(itFrame[,c("dimension", "model")]), unique(refPop), by.x = c("dimension", "model"), by.y = c("domain", "model"),all=FALSE)
                                stopifnot(nrow(mat)==1)
    ### 1. Transformation fuer Itemparameter
                                itFrame[,"estTransf"] <- itFrame[,"est"] + equatingList[["items"]][[mod]][[dimname]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dimname]][["method"]] ]]
    ### Achtung, heikel: wenn equatet wurde, aber der Datensatz aus der Normpopulation kommt, werden hier die empirischen Mittelwerte,
    ### die oben (mit oder ohne Gewichte) berechnet wurden, nochmal transformiert ... sollte praktisch nie der Fall sein.
                                if ( isTRUE(mr) ) {
                                     if ( equatingList[["items"]][[mod]][[dimname]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dimname]][["method"]] ]] != 0) {
                                          cat("W A R N I N G: Preceding Equating without 'refPop' definition. Sure you want to use current sample as drawn from the reference population?\n")
                                          mat[,3] <- mat[,3] + equatingList[["items"]][[mod]][[dimname]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dimname]][["method"]] ]]
                                     }
                                }
                                itFrame[,"estTransf625"]   <- itFrame[,"estTransf"] + log(0.625/(1-0.625))
                                itFrame[,"estTransfBista"] <- (itFrame[,"estTransf625"] - mat[,3]) / mat[,4] * mat[,6] + mat[,5]
                                if ( isFALSE(cutsMis) ) {
    ### Achtung: dieser Umweg ist notwendig, weil 'num.to.cat' Attribute ausgibt die unten wieder gebraucht werden!
                                     traitLevel            <- eatTools::num.to.cat(x = itFrame[,"estTransfBista"], cut.points = cuts[[mat1]][["values"]], cat.values = cuts[[mat1]][["labels"]])
                                     itFrame[,"traitLevel"]<- traitLevel
                                }
    ### Achtung!! Linkingfehler sollte eigentlich nur ausgegeben werden, wenn das vorherige equating NICHT durchgeschleift wurde!
                                itFrame[,"linkingConstant"]<- equatingList[["items"]][[mod]][[dimname]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dimname]][["method"]] ]]
                                itFrame[,"linkingMethod"]  <- equatingList[["items"]][[mod]][[dimname]][["method"]]
                                itFrame[,"nLinkitems"]     <- equatingList[["items"]][[mod]][[dimname]][["eq"]][["descriptives"]][["N.Items"]]
                                itFrame[,"linkingError"]   <- equatingList[["items"]][[mod]][[dimname]][["eq"]][["descriptives"]][["linkerror"]]
    ### Transformation des Linkingfehlers entsprechend der Rechenregeln fuer Varianzen. ist geprueft, dass dasselbe rauskommt, wie wenn man Parameter transformiert und dann Linkingfehler bestimmt
                                itFrame[,"linkingErrorTransfBista"] <- ( (itFrame[,"linkingError"]^2) * (mat[,6]^2) / (mat[,4]^2) )^0.5
    ### Deltamethode, wie in eatTrend (Funktion 'seKompstuf'). Dazu wird MW und SD der Fokuspopulation benoetigt! (wurde oben als 'msd' berechnet)
    ### das ganze findet nur statt, wenn sowohl cut scores bereits definiert sind und wenn equatet wurde (denn nur dann gibt es einen Linkingfehler, den man transformieren kann)
                            }
                            if (isRunM) {
                                pv  <- pvFromRes(resMD, toWideFormat = FALSE, idVarName=idVarName, verbose=FALSE)
                                equ <- equatingList[["items"]][[mod]][[dimname]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dimname]][["method"]] ]]
    ### Hotfix fuer bayesianisch
                                if (!exists("mat")) { mat <- refPop[match(dimname,  refPop[,"domain"]),] }
                                pv[,"valueTransfBista"] <- (pv[,"value"] + equ - mat[,3]) / mat[,4] * mat[,6] + mat[,5]
    ### Dazu muss zuerst Mittelwert und SD der Fokuspopulation bestimmt werden.
                                if (!is.null(pv)) {
                                    if ( is.null(weights) ) {
                                         msdF<- eatRep::repMean ( datL = pv, ID = id, imp = "imp", dependent = "valueTransfBista", na.rm = TRUE, verbose = FALSE, progress = FALSE)
                                         msdF<- adaptEatRepVersion(msdF)
                                    }  else  {
                                         pvF <- eatTools::mergeAttr ( pv, weights , by.x = id, by.y = colnames(weights)[1], all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "cases", xName = paste0("plausible values for dimension ",dimname), yName = "weights", verbose = c("match", "dataframe"))
                                         mis <- which(is.na(pvF[,colnames(weights)[2]]))
                                         if ( length(mis) > 0 ) {               ### missings in the weights frame are not allowed
                                              cat(paste ( "Found ",length(mis)," missing values in the 'weights' frame.\n    Cases with missing values on weighting variable will be ignored for transformation.\n",sep=""))
                                              pvF <- pvF[-mis,]
                                         }
                                         msdF<- eatRep::repMean ( datL = pvF, ID = id, imp = "imp", wgt = colnames(weights)[2], dependent = "valueTransfBista", na.rm = TRUE, verbose = FALSE, progress = FALSE)
                                         msdF<- adaptEatRepVersion(msdF)
                                    }
                                    msdFok <- c(msdF[intersect(which(msdF[,"parameter"] == "mean"), which(msdF[,"coefficient"] == "est")),"value"], msdF[intersect(which(msdF[,"parameter"] == "sd"), which(msdF[,"coefficient"] == "est")),"value"])
                                }  else  {
                                    message("Results object does not contain any plausible values. Skip transformation of linking error for competence levels.")
                                }
                                if ( !is.null(pv) && !is.null ( itFrame )) {    ### Cuts mit Schwelle nach unten und nach oben offen
                                    if ( cutsMis == FALSE & !is.null ( equatingList[["items"]] )) {
                                         cts <- c( -10^6, cuts[[mat1]][["values"]], 10^6)
                                         le  <- do.call("rbind", lapply ( (length(cts)-1):1 , FUN = function ( l ) {
                                                kmp<- c(cts[l], cts[l+1])       ### Linkingfehler fuer einzelnen Kompetenzintervalle; absteigend wie bei karoline
                                                a1 <- sum ( dnorm ( ( kmp - mat[,5]) / mat[,6] ) * c(-1,1) / mat[,6] )
                                                a2 <- sum ( dnorm ( ( kmp - msdFok[1]) / msdFok[2] ) * c(-1,1) / msdFok[2] )
    ### Achtung! der 'mutmassliche Fehler' kann auch auftreten, wenn das Equating zuvor durchgeschleift wurde und deshalb gar keine Linkingfehler berechnet werden koennen
                                                if(a2 == 0 ) {cat("mutmasslicher fehler.\n")}
                                                del<- ( (  a1^2 + a2^2 ) * (unique(itFrame[,"linkingErrorTransfBista"])^2) / 2  )^0.5
                       			                    del<- data.frame ( traitLevel = attr(traitLevel, "cat.values")[l], linkingErrorTraitLevel = del )
                       			                    return(del)}))
    ### ggf. weg! 'linkingErrorTraitLevel' ergibt ja fuer Items keinen Sinn, nur fuer Personenparameter
                                         ori <- colnames(itFrame)
                                         chk <- unique(le[,"traitLevel"]) %in% unique(itFrame[,"traitLevel"])
                                         if ( length( which(chk == FALSE)) > 0) {
                                             warning(paste("Model '",unique(itFrame[,"model"]),"', dimension '",unique(itFrame[,"dimension"]),"': No items on trait level(s) '",paste( unique(le[,"traitLevel"])[which(chk == FALSE)], collapse = "', '"), "'.", sep=""))
                                         }
                                         itFrame <- eatTools::mergeAttr ( itFrame, le, by = "traitLevel", sort = FALSE, all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "trait levels", xName = "item parameter list", yName = "linking error list", verbose = "match")
                                         itFrame <- itFrame[,c(ori, "linkingErrorTraitLevel")]
                                    }
                                    itFrame <- itFrame |> dplyr::mutate(refMean= mat[,3], refSD = mat[,4], refTransfMean=mat[,5], refTransfSD= mat[,6])
                                }
    ### 2. Transformation der Personenparameter: kann auch dann stattfinden, wenn PVs bayesianisch gezogen wurden
                                if(!exists("mat1") ) {mat1 <- match(dimname, names(cuts)); stopifnot(length(mat1)==1)}
                                if (!is.null(pv)) {
                                    if ( isFALSE(cutsMis) ) { pv[,"traitLevel"]   <- eatTools::num.to.cat(x = pv[,"valueTransfBista"], cut.points = cuts[[mat1]][["values"]], cat.values = cuts[[mat1]][["labels"]])}
                                    pv[,"dimension"]  <- pv[,"group"]
                                    if(!exists("le")) {
                                        warning("Skip check whether all competence levels are occupied (due to bayesian plausible values imputation).")
                                    }  else  {
                                        chk <- unique(le[,"traitLevel"]) %in% unique(pv[,"traitLevel"])
                                        if ( length( which(chk == FALSE)) > 0) {
                                             warning(paste("Model '",unique(itFrame[,"model"]),"', dimension '",unique(itFrame[,"dimension"]),"': No plausible values on trait level(s) '",paste( unique(le[,"traitLevel"])[which(chk == FALSE)], collapse = "', '"), "'.", sep=""))
                                        }
                                        stopifnot ( length( unique ( na.omit(itFrame[,"linkingErrorTransfBista"]))) %in% 0:1)
                                        pv[,"linkingError"] <- equatingList[["items"]][[mod]][[dimname]][["eq"]][["descriptives"]][["linkerror"]]
                                        pv[,"linkingErrorTransfBista"] <- unique ( itFrame[,"linkingErrorTransfBista"])
                                    }
                                    ori <- colnames(pv)                             ### nur wenn untere Bedingung == TRUE, gibt es das Objekt 'le', das gemergt werden soll
                                    if ( cutsMis == FALSE && !is.null ( equatingList[["items"]]) && exists("le") ) {
                                         pv  <- eatTools::mergeAttr ( pv, le, by = "traitLevel", sort = FALSE, all.x = TRUE, all.y = FALSE, setAttr = FALSE, unitName = "trait levels", xName = "plausible values", yName = "linking error list", verbose = c("match"))
                                         pv  <- pv[,c(ori, "linkingErrorTraitLevel")]
                                    }
    ### ggf. Gewichte an Personenframe mit dranhaengen
                                    if (!is.null(weights)) {
                                         pv  <- merge ( pv, weights , by.x = id, by.y = colnames(weights)[1], all.x = TRUE, all.y = FALSE)
                                    }
                                }  else  {
                                    msdFok <- c(NA, NA)
                                }
                            }
    ### 'refPop' Informationstabelle bauen
                            colnames(mat) <- c("domain", "model", "refMean", "refSD", "bistaMean", "bistaSD")
                            if (!isRunM) {pv <- NULL}
                            mat  <- cbind(mat, focusMean = msdFok[1], focusSD = msdFok[2])
                            return(list( itempars = itFrame, personpars = pv, rp = mat))
                      }  })                                                     ### untere Zeile: hier muss man rbind.fill nehmen, das gibt sonst im LV2021 einen Fehler, wenn bei der PV-Ziehung manche Items verankert sind, andere nicht
              itempars<- do.call(plyr::rbind.fill, lapply ( dimN, FUN = function ( x ) { x[["itempars"]]}))
              perspar <- do.call("rbind", lapply ( dimN, FUN = function ( x ) { x[["personpars"]]}))
              rp      <- do.call("rbind", lapply ( dimN, FUN = function ( x ) { x[["rp"]]}))
              return( list ( itempars = itempars, personpars = perspar, rp=rp)) } )
       personpars <- do.call("rbind", lapply ( modN, FUN = function ( x ) { x[["personpars"]]}))
       itempars   <- do.call(plyr::rbind.fill, lapply ( modN, FUN = function ( x ) { x[["itempars"]]}))
       rp         <- do.call("rbind", lapply ( modN, FUN = function ( x ) { x[["rp"]]}))
    ### jetzt noch die Itemparameterliste fuer die Vergleichsarbeiten reduzieren und aufbereiten
       if ( isFALSE(vera) || isFALSE(isRunM) ) {
           itemVera <- NULL
       }  else  {
           itemVera <- createItemVeraObj(itempars=itempars, roman=roman, results = equatingList[["results"]], q3bound=q3bound)
       }
    ### optional: separate linkingfehlerobjekte erzeugen
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
       ret        <- list ( itempars = itempars, personpars = personpars, refPop = refPop, means = rp, all.Names = context, itemparsVera = itemVera, linkingErrors = leo)
       class(ret) <- c("list", "transfBista")
       return( ret ) }

### Hilfsfunktion fuer 'transformToBista'
generateOrCheckRefPop <- function (equatingList, refPop, dims, mods, isRunM, id, weights, defaultM,defaultSD) {
       if(missing(refPop)) {                                                    ### 'refPop' erzeugen
           if(isRunM) {
               cat("'refPop' was not defined. Treat current sample as drawn from the reference population.\n")
               refPop <- do.call("rbind", lapply(mods, FUN = function ( mod ) { ### aeussere Schleife: geht ueber modelle
                         nam2 <- names(equatingList[["items"]][[mod]])          ### innere Schleife: geht ueber Dimensionen (innerhalb von modellen)
                         refL2<- do.call("rbind", lapply(nam2, FUN = function ( dimname ) {
                             rex  <- pvFromRes(equatingList[["results"]][unique(c(intersect(which(equatingList[["results"]][,"model"] == mod), which(equatingList[["results"]][,"group"] == dimname)),  which(equatingList[["results"]][,"type"] == "tech"))),], toWideFormat = FALSE, idVarName=id, verbose=FALSE)
                             if (is.null(rex)) {return(NULL)}                   ### NULL wird zurueckgegeben, wenn keine PVs in der Ergebnisstrauktur vorhanden waren
                             if ( is.null(weights) ) {
                                  msd <- eatRep::repMean ( datL = rex, ID = id, imp = "imp", dependent = "value", na.rm = TRUE, verbose = FALSE, progress = FALSE) |> adaptEatRepVersion()
                             }  else  {                                         ### in transformToBista() muessen die messages von mergeAttr nicht uber capture.output abgefangen werden wie in anker(), da transformToBista() immer nur single core aufgerufen werden kann
                                rex <- eatTools::mergeAttr ( rex, weights , by.x = id, by.y = colnames(weights)[1], all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "cases", xName = paste0("plausible values for dimension ",dimname), yName = "weights", verbose = c("match", "dataframe"))
                                mis <- which(is.na(rex[,colnames(weights)[2]]))
                                if ( length(mis) > 0 ) {                        ### missings in the weights frame are not allowed
                                     if(length(mis) == nrow(rex)) {stop(paste("Mergin of weights and plausible values for '", dimname, "' failed. No common units."))}
                                     cat(paste ( "Found ",length(mis)," missing values in the 'weights' frame.\n    Cases with missing values on weighting variable will be ignored for transformation.\n",sep=""))
                                     rex <- rex[-mis,]
                                }
                                msd <- eatRep::repMean ( datL = rex, ID = id, imp = "imp", wgt = colnames(weights)[2], dependent = "value", na.rm = TRUE, verbose = FALSE, progress = FALSE) |> adaptEatRepVersion()
                             }
                             rp <- data.frame ( domain = dimname , model = mod, m = msd[intersect(which(msd[,"parameter"] == "mean"), which(msd[,"coefficient"] == "est")),"value"], sd = msd[intersect(which(msd[,"parameter"] == "sd"), which(msd[,"coefficient"] == "est")),"value"])
                             return(rp)}))
                         return(refL2)}))
           }  else  {
               return(NULL)
           }
       }  else  {                                                               ### 'refPop' checken
           refPop <- eatTools::makeDataFrame(refPop)                            ### ab Spalte 2 muss alles numerisch sein
           lapply(refPop[,-1], checkmate::assert_numeric)
           checkmate::assert_character(refPop[,1], any.missing = FALSE, unique=TRUE)
           if(!all(refPop[,1] %in% dims)) {
              notIncl <- setdiff(dims,refPop[,1])
              stop(paste0("Following ",length(notIncl), " dimension(s) not included in 'refPop': '", paste(notIncl, collapse = "', '"),"'."))
           }
           refPop <- merge(expand.grid(model = mods, domain = dims), refPop, by.x = "domain", by.y = colnames(refPop)[1], all=TRUE)
       }
       if(ncol ( refPop ) < 5) {
           cat ( paste("The 'refPop' data.frame does not include information about reference population mean/SD on Bista metric. Values will be defaulted to ",defaultM,"/",defaultSD,".\n",sep=""))
           refPop <- data.frame(refPop, defaultMean = defaultM, defaultSD = defaultSD, stringsAsFactors = FALSE)
       } else {
           if ( ncol ( refPop) != 6 ) { stop ( "Invalid 'refPop'.\n") }
       }
       return(refPop)}
