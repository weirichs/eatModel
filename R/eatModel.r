



simEquiTable <- function ( anchor, mRef, sdRef, addConst = 500, multConst = 100, cutScores) {
                anchor<- eatTools::makeDataFrame(anchor)
                if ( ncol(anchor) != 2) {
                     warning(paste0("'anchor' has ",ncol(anchor)," columns. First column is used as item ID, second column is used as item parameter."))
                }
                if(!inherits(anchor[,2], c("integer", "numeric"))) {stop("Item parameter column must be numeric.")}
                if(length(unique(anchor[,1])) != nrow(anchor)) {stop("Item ID column has duplicated entries.")}
                dtmp  <- data.frame(rbind(1*(lower.tri(matrix(1, nrow = nrow(anchor), ncol = nrow(anchor)))),1))
                dtmp  <- data.frame(dtmp, score = rowSums(dtmp) , irtoys::wle(dtmp, cbind(1, anchor[,2], 0)), stringsAsFactors = FALSE)
                dtmp[,"bista"] <- (dtmp[,"est"] - mRef) / sdRef * multConst + addConst
                dtmp[,"ks"]    <- eatTools::num.to.cat ( x = dtmp[,"bista"], cut.points = cutScores[["values"]], cat.values = cutScores[["labels"]])
                dig   <- 0
                while ( length(which(round(dtmp[,"bista"], digits = dig) %in% cutScores[["values"]])) > 0) {dig <- dig + 1}
                shrt  <- do.call("rbind", by ( data = dtmp, INDICES = dtmp[,"ks"], FUN = function ( sks ) { data.frame ( score = paste(c(min(sks[,"score"]), max(sks[,"score"])), collapse=" bis "), estimate = paste(round(c(min(sks[,"est"]), max(sks[,"est"])),digits=2), collapse=" bis "), bista = paste(round(c(min(sks[,"bista"]), max(sks[,"bista"])),digits=dig), collapse=" bis "), ks=unique(sks[,"ks"]), stringsAsFactors=FALSE)}))
                return(list ( complete = dtmp[,c("score", "est", "bista", "ks")], short = shrt))}


checkItemParLists <- function (prmNorm, item, domain, testlet, value, dims = NULL) {
           if ( ncol ( prmNorm ) == 2 ) {
                if ( is.null(item) && is.null(value) ) {
                     item <- colnames(prmNorm)[1]
                     value<- colnames(prmNorm)[2]
                }
                if ( is.null(item) && !is.null(value) || !is.null(item) && is.null(value)) {
                     stop("If 'prmNorm' has two columns, either both 'item' and 'value' or none of them should be specified.")
                }
           }  else  {
                if ( is.null(item) || is.null(value)) { stop("If 'prmNorm' has more than two columns, 'item' and 'value' columns must be specified explicitly.") }
           }
           allF <- list(item=item, domain = domain, testlet=testlet, value = value)
           allF <- lapply(allF, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = prmNorm, variable=ii)})
           nomis<- sapply(prmNorm[,unlist(allF)], FUN = function ( i ) { length(which(is.na(i)))})
           if ( any(nomis>0)) {
                warning("Found ", length(which(nomis>0)), " column(s) in 'prmNorm' with missing values: '", paste(names(nomis[which(nomis>0)]), collapse= "', '"), "'")
           }
           tab  <- table(prmNorm[,c(allF[["item"]], allF[["domain"]]), drop=FALSE])
           if (!all(tab %in% 0:1)) {stop("Items must be unique for each domain in reference parameter frame 'prmNorm'.")}
           if(!inherits(prmNorm[,allF[["value"]]], "numeric")) {stop("Parameter value column in 'prmNorm' must be numeric.")}
           if (!is.null ( allF[["domain"]]) && !is.null(dims) ) {
                mis <- setdiff ( dims,  names(table(prmNorm[, allF[["domain"]] ])) )
                if ( length( mis ) > 0 ) { stop ( paste ( "Domain '",mis,"' is missing in 'prmNorm'.\n",sep="")) }
                uni <- by ( data = prmNorm, INDICES = prmNorm[, allF[["domain"]] ], FUN = function ( g ) {
                       if (!length(g[,allF[["item"]]]) == length(unique(g[,allF[["item"]]]))) { stop(paste ( "Item identifiers are not unique in 'prmNorm' for domain '",g[1,allF[["domain"]]],"'.\n",sep=""))}
                       }, simplify = FALSE)
           }
           return(allF)}

transformItemParListIntoResults <- function(results, itemF, domainF, testletF, valueF){
           allF <- checkItemParLists(prmNorm =results, item = itemF, domain = domainF, testlet = testletF, value = valueF)
           if (!is.null(allF[["domain"]])) {
                dims <- names(table(results[,allF[["domain"]]]))
           }  else  {
                allF[["domain"]] <- "domaene"
                results[, allF[["domain"]]] <- dims <- "global"
           }
           results[,"model"] <- results[, allF[["domain"]]]
           weg  <- intersect ( colnames (results ) , setdiff  ( c("item", "dimension", "est"), unlist(allF) ))
           if ( length ( weg ) > 0 )  {
                results <- results[, -match(weg, colnames(results))]
           }
           allF <- allF[which(!sapply(allF, is.null))]
           toRec<- lapply(names(allF), FUN = function ( ff ) { paste ( "'",allF[[ff]],"'='",car::recode(ff, "'item'='item'; 'domain'='dimension'; 'value'='est'"),"'",sep="")})
           toRec<- paste(toRec, collapse = "; ")
           colnames(results) <- car::recode (colnames(results), toRec)
           return(list(results=results, dims=dims))}

buildEmptyResultsObject <- function (d, method, results ) {
           it  <- itemFromRes(d)
           if ( "estOffset" %in% colnames ( it ) ) {
                 d[,"par"] <- car::recode ( d[,"par"], "'offset'='est'")
                 it <- itemFromRes(d)
           }
           if ( !is.null(it)) {
                 dimN <- by ( data = it, INDICES = it[,"dimension"], FUN = function ( prmDim ) {
                         eq <- list(B.est = c(Mean.Mean=0 , Haebara =0, Stocking.Lord=0), descriptives = c(N.Items =0, SD=NA,  Var=NA, linkerror=NA))
                         return ( list ( eq = eq, items = prmDim, method = method ) ) }, simplify = FALSE)
           }  else  {
                 resX <- results[which(!is.na(results[,"group"])),]
                 dimN <- by ( data = results, INDICES = results[,"group"], FUN = function ( prmDim ) {
                         eq <- list(B.est = c(Mean.Mean=0 , Haebara =0, Stocking.Lord=0), descriptives = c(N.Items =0, SD=NA,  Var=NA, linkerror=NA))
                         return ( list ( eq = eq, items = prmDim, method = method ) ) }, simplify = FALSE)
           }
           return(dimN)}

printToConsole <- function(d, nMods, it, prmDim, eq, allN, method, estimation, eqh, eqr, mess1) {
           cat(paste("\n",paste(rep("=",100),collapse=""),"\n \nModel No. ",match(d[1,"model"], names(nMods)),"\n    Model name:                ",d[1,"model"],"\n    Number of dimension(s):    ",length(unique(it[,"dimension"])),"\n    Name(s) of dimension(s):   ", paste( names(table(as.character(it[,"dimension"]))), collapse = ", "),"\n",sep=""))
           if  ( length(names(table(as.character(it[,"dimension"])))) > 1) {  cat(paste("    Name of current dimension: ",names(table(prmDim[,"dimension"]))," \n",sep=""))}
           cat(paste("    Number of linking items:   " , eq[["descriptives"]][["N.Items"]],"\n",sep=""))
           if ( !is.null(allN[["testlet"]]) ) { cat(paste( "    Number of testlets:        ",  eq[["ntl"]],"\n",sep="")) }
           if(!is.null(mess1)) {add <- " (excluding testlets)"} else {add <- NULL}
           cat(paste("    Linking method:            " , method,add, "\n",sep=""))
           if (method == "robust") { cat(paste("    Optimal trimming param.:   " , eqr[["kopt"]],"\n",sep="")) }
           if (method == "Haberman") {
               cat(paste("    Estimation method:         " , car::recode(estimation,"'OLS'='ordinary least squares'; 'BSQ'='bisquare weighted regression'; 'HUB'='regression using Huber weights'; 'MED'='median regression'; 'LTS'='trimmed least squares'; 'L1'='median polish'; 'L0'='minimizing number of interactions'"), "\n",sep=""))
               tf <- capture.output(summary(eqh))
               i1 <- grep("Used trimming factor", tf)
               i2 <- grep("Estimation information item intercepts", tf)
               i3 <- min(i1[which(i1>i2)])
               i4 <- unlist(strsplit(tf[i3], "="))
               cat(paste("    Used trimming factor:      " , round(as.numeric(eatTools::crop(i4[length(i4)])), digits = 3), "\n",sep=""))   }}

handleLinkingDif <- function(prmDim,prbl, eq, difBound, dif, method, excludeLinkingDif, iterativ,prmM, allN) {
               cat(paste ( "\nDimension '", prmDim[1,"dimension"], "': ", length( prbl), " of ", nrow( eq[["anchor"]]), " items with linking DIF > ",difBound," identified.\n",sep=""))
               dskr <- data.frame ( item = eq[["anchor"]][prbl,"item"], dif = dif[prbl], linking.constant = eq[["B.est"]][[method]], linkerror = eq[["descriptives"]][["linkerror"]] )
               if ( !excludeLinkingDif) { info<- dskr }
               if ( excludeLinkingDif ) {
                    if ( iterativ == FALSE ) {
                          cat(paste("   Exclude ",length( prbl), " items.\n",sep=""))
                          qp1 <- prmM[-match ( dskr[,"item"], prmM[,allN[["item"]]]),]
                          eq1 <- equAux ( x=prmDim[ ,c("item", "est")], y = qp1[,c(allN[["item"]], allN[["value"]], allN[["testlet"]])] )
                          info<- data.frame ( method = "nonIterativ", rbind ( data.frame ( itemExcluded = "" , linking.constant = eq[["B.est"]][[method]], linkerror = eq[["descriptives"]][["linkerror"]] ), data.frame ( itemExcluded = paste ( prmM[match ( dskr[,"item"], prmM[,allN[["item"]]]),allN[["item"]]] , collapse = ", "), linking.constant = eq1[["B.est"]][[method]], linkerror = eq1[["descriptives"]][["linkerror"]] ) ))
                          eq  <- eq1
                    }  else  {
                          info<- data.frame ( method = "iterativ", iter = 0 , itemExcluded = "" , DIF.excluded="", linking.constant = eq[["B.est"]][[method]], linkerror = eq[["descriptives"]][["linkerror"]] )
                          qp1 <- prmM
                          iter<- 1
                          while  ( length ( prbl ) > 0 ) {
                                   maxV<- eq[["anchor"]][,"TransfItempar.Gr1"] - eq[["anchor"]][,"Itempar.Gr2"]
                                   maxV<- maxV[which(abs(maxV) == max(abs(maxV)))]
                                   maxD<- which ( abs ( eq[["anchor"]][,"TransfItempar.Gr1"] - eq[["anchor"]][,"Itempar.Gr2"] ) == max ( abs (eq[["anchor"]][,"TransfItempar.Gr1"] - eq[["anchor"]][,"Itempar.Gr2"])) )
                                   wegI<- eq[["anchor"]][maxD,"item"]
                                   cat ( paste ( "   Iteration ", iter,": Exclude item '",wegI,"'.\n",sep=""))
                                   qp1 <- qp1[-match ( wegI, qp1[,allN[["item"]]]),]
                                   eq  <- equAux ( x = prmDim[ ,c("item", "est")], y = qp1[,c(allN[["item"]], allN[["value"]], allN[["testlet"]])] )
                                   dif <- eq[["anchor"]][,"TransfItempar.Gr1"] - eq[["anchor"]][,"Itempar.Gr2"]
                                   prbl<- which ( abs ( dif ) > difBound )
                                   info<- rbind(info, data.frame ( method = "iterativ", iter = iter , itemExcluded = wegI, DIF.excluded=as.character(round(maxV,digits=3)), linking.constant = round ( eq[["B.est"]][[method]],digits = 3), linkerror = round ( eq[["descriptives"]][["linkerror"]], digits = 3) ))
                                   iter<- iter + 1
                          }
                    }
               }
               return(list(eq=eq, info=info, info2=dskr))}

noLinkingDif <- function (method, eq, eqr, eqh) {
              if (method %in% c("Mean.Mean", "Haebara", "Stocking.Lord")) {
                   info <- data.frame ( linking.constant = eq[["B.est"]][[method]], linkerror = eq[["descriptives"]][["linkerror"]] )
              }
              if ( method == "robust" ) {
                   werte<- list (c(eqr[["meanpars"]][["k0"]], eqr[["meanpars"]][[names(eqr[["ind.kopt"]])]]), c(eqr[["se"]][["k0"]], eqr[["se"]][[names(eqr[["ind.kopt"]])]]))
                   werte<- lapply( werte , FUN = function ( wert) { ifelse(is.null(wert), NA, wert)})
                   info <- data.frame ( linking.method = c("non robust", "robust"), linking.constant = werte[[1]], linkerror = werte[[2]], stringsAsFactors=FALSE)
              }
              if ( method == "Haberman" ) {
                   wert <- eqh[["transf.itempars"]][2,"B_b"]
                   wert <- ifelse(is.null(wert), NA, wert)
                   info <- data.frame ( linking.constant = wert, linkerror = NA, stringsAsFactors=FALSE)
              }
              return(list(eq=eq, info=info))}

createOutput <- function (method, eqr, prm, eqh, info){
              if (method == "robust") {
                  wert<- list(eqr[["meanpars"]][[names(eqr[["ind.kopt"]])]], eqr[["se"]][[names(eqr[["ind.kopt"]])]])
                  wert<- lapply(wert, FUN = function (w) {ifelse(is.null(w), NA, w)})
                  eq  <- list ( B.est = data.frame ( robust = wert[[1]]), descriptives = data.frame ( N.items = nrow(prm), linkerror = wert[[2]], stringsAsFactors =FALSE) )
              }
              if (method == "Haberman") {
                  wert <- list(eqh[["transf.itempars"]][2,"B_b"], nrow(eqh[["joint.itempars"]]))
                  wert <- lapply(wert, FUN = function (w) {ifelse(is.null(w), NA, w)})
                  eq <- list ( B.est = data.frame ( Haberman = wert[[1]]), descriptives = data.frame ( N.items = wert[[2]], linkerror = NA, stringsAsFactors =FALSE) )
              }
              return(list(eq=eq, info=info))}


equAux  <- function ( x, y ) {
           eq  <- sirt::equating.rasch(x = x, y = y[,1:2])
           if ( ncol(y)==3) {
                colnames(x)[1] <- colnames(y)[1] <- "item"
                dfr <- merge( x, y, by = "item", all = FALSE)
                stopifnot ( ncol ( dfr ) == 4 )
                if ( nrow ( dfr ) < 1 ) { stop ( "No common items for linking.\n")}
                txt <- capture.output ( eqJk<- sirt::equating.rasch.jackknife(dfr[ , c(4 , 2  , 3 , 1 ) ], display = FALSE ) )
                if(!all ( unlist(lapply(txt, nchar)) == 0  ) ) { cat(txt, sep="\n")}
                eq[["descriptives"]][["linkerror"]] <- eqJk[["descriptives"]][["linkerror.jackknife"]]
                eq[["ntl"]]  <- length(unique(dfr[,4]))
           }
           return(eq)}

adaptEatRepVersion <- function ( x ) {
     if ( inherits(x, "data.frame"))  {
           return ( x )
     }  else  {
           x <- x[[1]][[1]]
           stopifnot ( inherits(x, "data.frame") )
           return(x)
     } }

createLinkingErrorObject <- function (itempars, years) {
     res <- do.call("rbind", by(data = itempars, INDICES = itempars[,"dimension"], FUN = function (d) {
            r1 <- do.call("rbind", lapply(c("value", "valueTransfBista", "traitLevel"), FUN = function (av) {
                  if ( av %in% c("value", "valueTransfBista")) {
                      prm <- "mean"
                      le  <- unique(d[,car::recode(av, "'value'='linkingError'; 'valueTransfBista'='linkingErrorTransfBista'")])
                      stopifnot(length(le) == length(prm))
                  } else {
                      dat <- unique(d[,c("traitLevel", "linkingErrorTraitLevel")])
                      stopifnot(length(dat[,1]) == length(unique(dat[,1])))
                      prm <- dat[,"traitLevel"]
                      le  <- dat[,"linkingErrorTraitLevel"]
                  }
                  dfr <- data.frame ( trendLevel1 = years[1], trendLevel2 = years[2], depVar = av, domain = d[1,"dimension"], parameter = prm, linkingError = le, stringsAsFactors = FALSE)
                  return(dfr)}))
            return(r1)}))
     return(res)}

createItemVeraObj <- function(itempars, roman, results, q3bound){
       pCols      <- colnames(itempars)[grep("^itemP", colnames(itempars))]
       allCols    <- na.omit(match ( c("dimension","item", pCols, "itemDiscrim", "estTransf", "infit", "estTransfBista", "traitLevel"), colnames(itempars)))
       itemVera   <- itempars[,allCols]
       colnames(itemVera) <- car::recode ( colnames(itemVera), "'dimension'='domain'; 'item'='iqbitem_id'; 'itemDiscrim'='trennschaerfe'; 'estTransf'='logit'; 'estTransfBista'='bista'; 'traitLevel'='kstufe'")
       colnames(itemVera)[match(pCols, colnames(itemVera))] <- paste0("lh", eatTools::removePattern ( string = pCols, pattern = "^itemP"))
       if ( roman == TRUE ) {
            if (!all(itemVera[,"kstufe"] %in% c("1a", "1b", 1:5))) {stop(paste("Competence levels do not match allowed values. '1a', '1b', '1', '2', '3', '4', '5' is allowed. '",paste(names(table(itemVera[,"kstufe"])), collapse = "', '"),"' was found.\n",sep=""))}
            itemVera[,"kstufe"] <- car::recode (itemVera[,"kstufe"], "'1a'='Ia'; '1b'='Ib'; '1'='I'; '2'='II'; '3'='III'; '4'='IV'; '5'='V'")
       }
       if ( length ( unique ( itemVera[,"iqbitem_id"])) != length ( itemVera[,"iqbitem_id"]) ) {
            cat("Found duplicated entries in 'item-ID' column. This should only occur for subject 'math' in grade 3.\n")
            tab  <- table(itemVera[,c("domain", "iqbitem_id")])
            if ( !"GL" %in% rownames(tab)) {
                 cat("Cannot find 'global' entry in the 'domain' column. Cancel reshaping.\n")
            }  else  {
                 if ( !sum(tab[which(rownames(tab) == "GL"),]) == ncol(tab)) {
                     cat("Found items without values on the 'global' domain. Cancel reshaping.\n")
                 }  else  {
                     if ( !all(colSums(tab) == 2) ) {
                         cat("Found items which do not have one 'global' and one domain-specific parameter. Cancel reshaping.\n")
                     }  else  {
                         itemVera[,"dummy"] <- car::recode ( itemVera[,"domain"], "'GL'='GL'; else = 'domain'")
                         colsValid <- c("lh", "trennschaerfe", "logit", "infit", "bista", "kstufe")
                         colsValid <- colsValid[which(colsValid %in% colnames(itemVera))]
                         long      <- reshape2::melt ( itemVera, id.vars = c("iqbitem_id", "dummy"), measure.vars = colsValid, na.rm=TRUE)
                         itemVera  <- suppressWarnings(eatTools::asNumericIfPossible(reshape2::dcast ( long , iqbitem_id ~ dummy + variable, value.var = "value"), force.string = FALSE))
                     }
                 }
            }
       }
       if ( "q3" %in% results[,"par"]) {
            itemVera   <- addQ3(dfr=itemVera, results=results, q3bound=q3bound)
       }
       return(itemVera) }

addQ3 <- function (dfr, results, q3bound) {
       q3  <- q3FromRes(results, out="long")
       q3  <- do.call(plyr::rbind.fill, lapply(names(q3), FUN = function (nq3) {
              x <- q3[[nq3]][which(abs(q3[[nq3]][,"value"]) > q3bound),]
              x <- suppressWarnings(eatTools::asNumericIfPossible(do.call(plyr::rbind.fill, by(x, INDICES = x[,"var1"], FUN = function (y) {
                   mat <- matrix(as.vector(unlist(t(y[,-1]))), nrow=1)
                   colnames(mat) <- paste(rep(c("q3item","q3value"), times=ncol(mat)/2), rep(1:(ncol(mat)/2), each=2), sep="_")
                   ret <- data.frame ( domain = nq3, iqbitem_id = unique(y[,"var1"]),mat, stringsAsFactors=FALSE)
                   return(ret)})), force.string=FALSE))
              return(x)}))
       dfr <- merge(dfr, q3, by=c("domain", "iqbitem_id"), all=TRUE)
       return(dfr)}
















eapFromRes <- function ( resultsObj, idVarName = NULL, verbose = TRUE ) {
          eapRo<- setdiff(intersect( which(resultsObj[,"par"] == "eap"),which(resultsObj[,"indicator.group"] == "persons")), which(resultsObj[,"derived.par"] == "rel"))
          id   <- unique(resultsObj[intersect(which(resultsObj[,"type"] == "tech"), which(resultsObj[,"par"] == "ID")),"derived.par"])
          id   <- getIdVarName(id, idVarName, verbose=verbose)
          if ( length ( eapRo ) == 0 ) {
               warning("'resultsObj' does not contain any eap values.")
               return ( NULL )
          }  else  {
             sel  <- resultsObj[eapRo,]
             sel  <- do.call("rbind", by(sel, INDICES = sel[,c("model", "group")], FUN = function ( gr ) {
                     res  <- reshape2::dcast ( gr , model+group+var1~derived.par, value.var = "value")
                     colnames(res)[-c(1:2)] <- c(id, "EAP", "SE.EAP")
                     weg  <- match(c("model", id), colnames(res))
                     res  <- data.frame ( res[,c("model", id)], dimension = as.character(gr[1,"group"]), res[,-weg,drop=FALSE], stringsAsFactors = FALSE)
                     return(res)}))
             return(sel)
          }  }

regcoefFromRes <- function (resultsObj, digits = NULL){
          regRo<- which(resultsObj[,"type"] == "regcoef")
          if(length(regRo)==0) {
              cat("No regression coefficients found in results object.\n")
              return(NULL)
          }  else  {
              re <- resultsObj[regRo,]
              re <- by(re, INDICES = re[,c("model", "group")], FUN = function (m) {
                    m[,"derived.par"] <- car::recode(m[,"derived.par"], "NA='est'")
                    mw <- reshape2::dcast(m, var1~derived.par, value.var="value")
                    mw[,"p"]   <- 2*(1-pnorm(abs(mw[,"est"] / mw[,"se"])))
                    mw[,"sig"] <- eatTools::num.to.cat(mw[,"p"], cut.points = c(0.001, 0.01, 0.05, 0.1), cat.values = c("***", "**", "**", ".", ""))
                    colnames(mw)[1] <- "parameter"
                    if(!is.null(digits)) {mw <- eatTools::roundDF(mw, digits =digits)}
                    return(mw)}, simplify=FALSE)
              nam<- eatTools::facToChar(expand.grid(attr(re, "dimnames")[["model"]], attr(re, "dimnames")[["group"]]))
              names(re) <- paste0("model: '",nam[,1],"', group: '",nam[,2],"'")
              re <- re[lengths(re) != 0]
              return(re)
          }}

correlationFromRes <- function (resultsObj, digits = NULL){
          corRo<- which(resultsObj[,"par"] == "correlation")
          if(length(corRo)==0) {
              cat("No correlation coefficients found in results object.\n")
              return(NULL)
          }  else  {
              re <- resultsObj[corRo,]
              re <- by(re, INDICES = re[,"model"], FUN = function (m) {
                    mw <- reshape2::dcast(m, var1~var2, value.var="value")
                    if(!is.null(digits)) {mw <- eatTools::roundDF(mw, digits =digits)}
                    return(mw)}, simplify=FALSE)
              return(re)
          }}

pvFromRes  <- function ( resultsObj, toWideFormat = TRUE, idVarName = NULL, verbose=TRUE) {
          pvRow<- intersect( which(resultsObj[,"par"] == "pv"),which(resultsObj[,"indicator.group"] == "persons"))
          if ( length ( pvRow ) == 0 ) {
               warning("'resultsObj' does not contain any pv values.")
               return ( NULL )
          }  else  {
             sel  <- resultsObj[pvRow, ]
             id   <- unique(resultsObj[intersect(which(resultsObj[,"type"] == "tech"), which(resultsObj[,"par"] == "ID")),"derived.par"])
             id   <- getIdVarName(id, idVarName, verbose=verbose)
             if (toWideFormat == TRUE ) {
                 sel  <- do.call("rbind", by(sel, INDICES = sel[,c("model","group")], FUN = function ( gr ) {
                         res  <- reshape2::dcast ( gr , model+var1~derived.par, value.var = "value")
                         colnames(res)[2] <- id
                         weg  <- match(c("model", id), colnames(res))
                         res  <- data.frame ( res[,c("model", id)], dimension = as.character(gr[1,"group"]), res[,-weg,drop=FALSE], stringsAsFactors = FALSE)
                         return(res)}))
             }  else  {
                 sel  <- sel[,c("model", "var1", "group", "derived.par", "value")]
                 recSt<- paste ( "'var1'='",id,"'; 'derived.par'='imp'",sep="")
                 colnames(sel) <- car::recode ( colnames(sel), recSt)
             }
             return(sel)
         }  }

getIdVarName <- function ( id, idVarName, verbose=TRUE) {
          if (length( id ) == 0 ) {
              if ( is.null(idVarName)) { new <- "idstud"} else { new <- idVarName}
              if(verbose){warning(paste0("Cannot identify student identifier variable (possibly because 'resultsObj' was created by an older version of 'eatModel'). student id variable will be defaulted to '",new,"'."))}
              id <- new
          }
          return(id)}

itemFromRes<- function ( resultsObj ) {
          res <- do.call(plyr::rbind.fill, by ( data = resultsObj, INDICES = resultsObj[,"model"], FUN = function ( mod ) {
                 sel  <- mod[intersect( which(mod[,"par"] %in% c("est", "estSlope", "Nvalid", "itemP", "ptBis", "itemDiscrim", "offset")),which(mod[,"indicator.group"] == "items")),]
                 if (nrow(sel)==0) {
                     return(NULL)
                 }  else  {
                     isDif<- intersect(which(mod[,"type"] == "tech"), which(mod[,"par"] == "DIF.var"))
                     if ( length( isDif ) > 0 ) {
                           vars     <- mod[intersect(which(mod[,"type"] == "tech"),which(mod[,"par"] == "variablen")),"derived.par"]
                           itemList <- do.call("rbind", lapply ( vars, FUN = function ( v ) {
                                       ind <- grep( paste0("_",v,"_"), sel[,"var1"])
                                       it  <- sort ( unique ( sel[ind,"var1"]))
                                       if(length(it)>2) {
                                          warning(paste0("DIF variable '",mod[isDif,"derived.par"],"' seems to have more than two categories. To date, this is not supported by 'eatModel'."))
                                       }
                                       return ( data.frame ( item = v, dif = it[1], weg = it[length(it)] , stringsAsFactors = FALSE) ) }))
                           weg      <- eatTools::whereAre ( itemList[,"weg"], sel[,"var1"], verbose=FALSE)
                           forDif   <- eatTools::whereAre ( itemList[,"dif"], sel[,"var1"], verbose=FALSE)
                           stopifnot(length( intersect(weg, forDif)) == 0 )
                           selForDif<- sel[forDif, ]
                           sel      <- sel[-c(weg, forDif) , ]
                           sel      <- sel[which ( sel[,"par"] != "ptBis" ) , ]
                           selDIF   <- do.call("rbind", by(selForDif, INDICES = selForDif[,"group"], FUN = function ( gr ) {
                                       res  <- reshape2::dcast ( gr , model+var1~par+derived.par, value.var = "value")
                                       mat  <- lapply( vars, FUN = function ( v ) { grep(paste0("_",v,"_"), res[,"var1"])})
                                       stopifnot (  all ( sapply(mat, length) == 1) )
                                       res[unlist(mat),"item"]  <- vars
                                       colnames(res) <- car::recode ( colnames(res) , "'est_infit'='infitDif'; 'est_se'='seDif'; 'est_NA'='estDif'")
                                       res[,"absDif"]<- abs ( res[,"estDif"]  * 2 )
                                       pval <- intersect(intersect(which(mod[,"type"] == "tech"), which(mod[,"par"] == "dif")), which(mod[,"derived.par"] == "p.value"))
                                       stopifnot (length(pval) == 1)
                                       pval <- mod[pval, "value"]
                                       adb  <- mod[intersect(intersect(which(mod[,"type"] == "tech"), which(mod[,"par"] == "dif")), which(mod[,"derived.par"] == "abs.dif.bound")),"value"]
                                       sdb  <- mod[intersect(intersect(which(mod[,"type"] == "tech"), which(mod[,"par"] == "dif")), which(mod[,"derived.par"] == "sig.dif.bound")),"value"]
                                       res[,paste("CI__", pval ,"__lb",sep="")] <- res[,"absDif"] - 2*abs(qnorm(0.5*(1-pval))) * res[,"seDif"]
                                       res[,paste("CI__", pval ,"__ub",sep="")] <- res[,"absDif"] + 2*abs(qnorm(0.5*(1-pval))) * res[,"seDif"]
                                       res  <- data.frame ( res, do.call("rbind", apply(res[,c("absDif", "seDif", paste("CI__",pval,"__lb",sep=""), paste("CI__",pval,"__ub",sep=""))], MARGIN = 1, FUN = function ( d ) {
                                                           check <- all ( !is.na(d) )
                                                           if(check == TRUE) {
                                                              crit1 <- d[["absDif"]] > adb
                                                              crit2 <- !all ( sort ( c ( d[[paste("CI__",pval,"__lb",sep="")]], sdb , d[[paste("CI__",pval,"__ub",sep="")]]), index.return = TRUE)$ix == 1:3 )
                                                              if ( crit1 == TRUE & crit2 == TRUE) { res <- 1 }  else { res <- 0}
                                                              ets   <- "A"
                                                              ets1  <- d[["absDif"]] > 0.43
                                                              ets2  <- !all ( sort ( c ( d[[paste("CI__",pval,"__lb",sep="")]], 0 , d[[paste("CI__",pval,"__ub",sep="")]]), index.return = TRUE)$ix == 1:3 )
                                                              if ( ets1 == TRUE & ets2 == TRUE) { ets <- "B" }
                                                              etsC1 <- d[["absDif"]] > 0.64
                                                              etsC2 <- !all ( sort ( c ( d[[paste("CI__",pval,"__lb",sep="")]], 0.43 , d[[paste("CI__",pval,"__ub",sep="")]]), index.return = TRUE)$ix == 1:3 )
                                                              if ( etsC1 == TRUE & etsC2 == TRUE) { ets <- "C" }
                                                              res   <- data.frame(difIndex = res, ETS = ets )
                                                           }  else  {
                                                              res   <- data.frame(difIndex = NA, ETS = NA )
                                                           }
                                                           return(res)}) ) )
                                       return(res)}))
                     }
                     sel  <- do.call(plyr::rbind.fill, by(sel, INDICES = sel[,"group"], FUN = function ( gr ) {
                             sfp  <- intersect ( which ( gr[,"par"] == "itemP"), which ( !is.na(gr[,"var2"])))
                             if ( length ( sfp ) > 0 ) {
                                  res  <- reshape2::dcast ( gr[-sfp,] , model+var1~par+derived.par, value.var = "value")
                                  sfp  <- reshape2::dcast ( gr[sfp,] , model+var1~par+var2, value.var = "value")
                                  res  <- merge ( res, sfp, by = c("model", "var1"), all = TRUE)
                             }  else  {
                                  res  <- reshape2::dcast ( gr , model+var1~par+derived.par, value.var = "value")
                             }
                             colnames(res) <- car::recode ( colnames(res) , "'var1'='item'; 'est_infit'='infit'; 'est_outfit'='outfit'; 'est_se'='se'; 'est_NA'='est'; 'estSlope_se'='seSlope'; 'estSlope_NA'='estSlope'; 'offset_NA'='estOffset'; 'Nvalid_NA'='Nvalid'; 'ptBis_NA'='ptBis'; 'itemP_NA'='itemP'; 'itemDiscrim_NA'='itemDiscrim'")
                             cols <- c("Nvalid", "itemP", "itemDiscrim", "est", "estOffset", "se", "estSlope", "seSlope", "infit","outfit", "ptBis")
                             drin1<- which(cols %in% colnames(res))
                             drin2<- grep("ptBis_", colnames(res))
                             drin3<- grep("itemP", colnames(res))
                             res  <- data.frame ( res[,c("model", "item")], dimension = as.character(gr[1,"group"]), res[,c(cols[drin1], colnames(res)[drin2] , setdiff (colnames(res)[drin3], cols[drin1])),drop=FALSE], stringsAsFactors = FALSE)
                             return(res)}))
                     if ( length( isDif ) > 0 ) {
                           ciCo<- colnames(selDIF)[grep("^CI__", colnames(selDIF))]
                           sel <- merge(sel, selDIF[,c("item", "model", "estDif", "seDif", "infitDif", "absDif", ciCo, "difIndex", "ETS")], by=c("item","model"), all=TRUE)
                     }
                     return(sel)
                 }
          }))
          return (res )}

q3FromRes<- function ( resultsObj, out = c("wide", "long" ), triangular = FALSE ) {
       out   <- match.arg(arg = out, choices = c("wide", "long" ))
       selM  <- by(data = resultsObj, INDICES = resultsObj[,"model"], FUN = function ( mr ) {
                sel  <- mr[which(mr[,"par"] == "q3"),]
                if ( nrow(sel)>0) {
                     if ( out == "wide") {
                          sel  <- reshape2::dcast( sel, var1~var2, value.var = "value")
                          if (triangular ) {sel <- eatTools::makeTria(sel)}
                     }  else  {
                          sel  <- sel[,c("var1", "var2", "value")]
                     }
                } else  { sel <- NULL }
                return(sel)})
       return(selM)}

wleRelFromRes <- function(resultsObj) {
          ret <- resultsObj[intersect(which(resultsObj[,"derived.par"] == "rel"), which(resultsObj[,"par"] == "wle")),c("model", "group", "value")]
          colnames(ret) <- car::recode(colnames(ret), "'value'='rel'; 'group'='domain'")
          return(ret)}

eapRelFromRes <- function(resultsObj) {
          ret <- resultsObj[intersect(which(resultsObj[,"derived.par"] == "rel"), which(resultsObj[,"par"] == "eap")),c("model", "group", "value")]
          colnames(ret) <- car::recode(colnames(ret), "'value'='rel'; 'group'='domain'")
          return(ret)}


wleFromRes <- function ( resultsObj , idVarName = NULL, verbose=TRUE) {
          wleRo<- setdiff(intersect( which(resultsObj[,"par"] %in% c("wle","NitemsSolved", "NitemsTotal")),which(resultsObj[,"indicator.group"] == "persons")), which(resultsObj[,"derived.par"] == "rel"))
          if(length(wleRo) == 0 ) {
             warning("'resultsObj' does not contain any WLE values.")
             return(NULL)
          }  else  {
             sel  <- resultsObj[wleRo,]
             sel  <- do.call("rbind", by(sel, INDICES = sel[,c("model", "group")], FUN = function ( gr ) {
                     res  <- reshape2::dcast ( gr , model+var1~par+derived.par, value.var = "value")
                     id   <- resultsObj[intersect(intersect(which(resultsObj[,"model"] == gr[1,"model"]),which(resultsObj[,"type"] == "tech")), which(resultsObj[,"par"] == "ID")),"derived.par"]
                     id   <- getIdVarName(id, idVarName, verbose=verbose)
                     recSt<- paste("'var1'='",id,"'; 'NitemsSolved_NA'='NitemsSolved'; 'NitemsTotal_NA'='NitemsTotal'",sep="")
                     colnames(res) <- car::recode ( colnames(res) , recSt)
                     weg  <- match(c("model", id), colnames(res))
                     res  <- data.frame ( res[,c("model", id),], dimension = as.character(gr[1,"group"]), res[,-weg,drop=FALSE], stringsAsFactors = FALSE)
                     return(res)}))
             return(sel)
          }  }

get.plausible <- function(file, quiet = FALSE, forConquestResults = FALSE)  {
                 input           <- scan(file,what="character",sep="\n",quiet=TRUE)
                 input           <- strsplit(eatTools::crop(gsub("-"," -",input) ) ," +")
                 n.spalten       <- max ( sapply(input,FUN=function(ii){ length(ii) }) )
                 input           <- data.frame( matrix( t( sapply(input,FUN=function(ii){ ii[1:n.spalten] }) ),length(input), byrow = FALSE), stringsAsFactors = FALSE)
                 pv.pro.person   <- sum (input[-1,1]==1:(nrow(input)-1) )
                 n.person        <- nrow(input)/(pv.pro.person+3)
                 weg             <- c(1, as.numeric( sapply(1:n.person,FUN=function(ii){((pv.pro.person+3)*ii-1):((pv.pro.person+3)*ii+1)}) ) )
                 cases           <- input[(1:n.person)*(pv.pro.person+3)-(pv.pro.person+2),1:2]
                 input.sel       <- input[-weg,]
                 n.dim <- dim(input.sel)[2]-1
                 if(quiet == FALSE) {cat(paste(n.person,"persons and",n.dim,"dimensions(s) found.\n"))
                               cat(paste(pv.pro.person,"plausible values were drawn for each person on each dimension.\n"))}
                 ID              <- input[  (pv.pro.person + 3) *  (1:n.person) - (pv.pro.person + 2) ,2]
                 colnames(input.sel) <- c("PV.Nr", paste("dim.",1:(ncol(input.sel)-1),sep=""))
                 input.sel[,1]   <- gsub( " ", "0", formatC(input.sel[,1],width = max(nchar(input.sel[,1]))))
                 input.sel$ID    <- rep(ID, each = pv.pro.person)
                 is.na.ID        <- FALSE
                 if(is.na(input.sel$ID[1])) {
                    is.na.ID        <- TRUE
                    input.sel$ID    <- rep( 1: n.person, each = pv.pro.person)
                 }
                 input.melt      <- reshape2::melt(input.sel, id.vars = c("ID", "PV.Nr") , stringsAsFactors = FALSE)
                 input.melt[,"value"] <- as.numeric(input.melt[,"value"])
                 input.wide      <- data.frame( case = gsub(" ", "0",formatC(as.character(1:n.person),width = nchar(n.person))) , reshape2::dcast(input.melt, ... ~ variable + PV.Nr) , stringsAsFactors = FALSE)
                 colnames(input.wide)[-c(1:2)] <- paste("pv.", paste( rep(1:pv.pro.person,n.dim), rep(1:n.dim, each = pv.pro.person), sep = "."), sep = "")
                 weg.eap         <- (1:n.person)*(pv.pro.person+3) - (pv.pro.person+2)
                 input.eap    <- input[setdiff(weg,weg.eap),]
                 input.eap    <- na.omit(input.eap[,-ncol(input.eap),drop=FALSE])
                 stopifnot(ncol(input.eap) ==  n.dim)
                 input.eap    <- lapply(1:n.dim, FUN=function(ii) {matrix(unlist(as.numeric(input.eap[,ii])), ncol=2,byrow = TRUE)})
                 input.eap    <- do.call("data.frame",input.eap)
                 colnames(input.eap) <- paste(rep(c("eap","se.eap"),n.dim), rep(paste("Dim",1:n.dim,sep="."),each=2),sep="_")
                 PV           <- data.frame(input.wide,input.eap, stringsAsFactors = FALSE)
                 numericColumns <- grep("pv.|eap_|case",colnames(PV))
                 if(is.na.ID == TRUE) {PV$ID <- NA}
                 for (ii in numericColumns) {PV[,ii] <- as.numeric(as.character(PV[,ii]))  }
                 if(  forConquestResults == TRUE ) {
                      return(list ( pvWide = PV, pvLong = input.melt, eap = input.eap))
                 }  else {
                 return(PV)}}

get.wle <- function(file)      {
            input <- eatTools::crop(scan(file, what = "character", sep = "\n", quiet = TRUE))
            input <- strsplit(input," +")
            n.spalten <- max ( sapply(input,FUN=function(ii){ length(ii) }) )
            n.wle <- floor((n.spalten-1) / 4)
            input <- suppressWarnings(eatTools::asNumericIfPossible(data.frame( matrix( t( sapply(input,FUN=function(ii){ ii[1:n.spalten] }) ),length(input),byrow = FALSE), stringsAsFactors = FALSE), force.string = FALSE))
            valid <- na.omit(input)
            cat(paste("Found valid WLEs of ", nrow(valid)," person(s) for ", n.wle, " dimension(s).\n",sep=""))
            if (nrow(valid) != nrow(input)) { cat(paste("    ",nrow(input)-nrow(valid)," persons with missings on at least one latent dimension.\n",sep="")) }
            namen1<- c(rep ( x = c("n.solved", "n.total"), times = n.wle), rep ( x = c("wle", "std.wle"), times = n.wle))
            namen2<- rep(rep ( paste(".", 1:n.wle, sep=""), each = 2),2)
            colnames(valid)[(ncol(valid)-length(namen2)):1] <- c("ID","case")[1:(ncol(valid)-length(namen2))]
            colnames(valid)[(ncol(valid)-length(namen2)+1):ncol(valid)] <- paste(namen1,namen2,sep="")
            return(valid)}

get.shw <- function(file, dif.term, split.dif = TRUE, abs.dif.bound = 0.6, sig.dif.bound = 0.3, p.value = 0.9) {
            all.output<- list();   all.terms <- NULL
            input.all <- scan(file,what="character",sep="\n",quiet=TRUE)
            rowToFind <- c("Final Deviance","Total number of estimated parameters")
            rowToFind <- sapply(rowToFind, FUN = function(ii) {
                         row.ii <- grep(ii,input.all)
                         stopifnot(length(row.ii) == 1)
                         row.ii <- as.numeric(unlist(lapply (strsplit(input.all[row.ii], " +"), FUN=function(ll) {ll[length(ll)]}) ))
                         return(row.ii)})
            ind       <- grep("TERM",input.all)
            grenzen   <- grep("An asterisk",input.all)
            if(length(ind)==0) {stop(paste("No TERM-statement found in file ",file,".\n",sep=""))}
            for (i in 1:length(ind)) {
                 term <- eatTools::crop(unlist(strsplit(input.all[ind[i]], ":"))[2])
                 cat(paste0("Found TERM ",i,": '",term,"' \n"))
                 all.terms <- c(all.terms,term)
                 bereich <- (ind[i]+6) : (grenzen[i] -2)
                 namen   <- gsub("\\^","",c("No.", strsplit(input.all[bereich[1]-2]," +")[[1]][-1]))
                 namen   <- rep(namen, car::recode(namen, "'CI'=2; else=1"))    ### Wenn ein "CI" als Spaltenname erscheint, muessen daraus im R-Dataframe zwei Spalten werden!
                 inp.sel <- gsub("\\(|)|,"," ",eatTools::crop(input.all[bereich]))
                 inp.sel <- gsub("\\*    ", "  NA", inp.sel)
                 if(!is.null(attr(file, "allNames"))) {
                     inp.sel <- eatTools::gsubAll(inp.sel, old = attr(file, "allNames")[["variablen"]], new = paste0(attr(file, "allNames")[["variablen"]], " "))
                 }
                 foo        <- strsplit(inp.sel," +")
                 maxColumns <- max(sapply(foo, FUN=function(ii){ length(ii)}))
                 nDifferentColumns <- length( table(sapply(foo, FUN=function(ii){ length(ii)  })))
                 maxColumns <- which( sapply(foo, FUN=function(ii){ length(ii) == maxColumns  }) )[1]
                 foo.2      <- which( sapply(1:nchar(inp.sel[maxColumns]),FUN=function(ii){u <- substr(inp.sel[maxColumns],ii,ii); b <- u==" "  }) )
                 foo.3      <- diff(foo.2)
                 foo.3      <- foo.2[foo.3 !=1]
                 ESTIMATE   <- which( sapply(1:nchar(input.all[ind[i] + 4] ),FUN=function(ii){u <- substr(input.all[ind[i] + 4],ii,ii+7); b <- u=="ESTIMATE"  }) )
                 foo.3      <- foo.3[foo.3>(ESTIMATE-3)]
                 if(nDifferentColumns>1) {
                    if(length(foo.3)>0) {
                       for (ii in 1:length(inp.sel)) {
                            for (iii in 1:length(foo.3)) {
                                 if(substr( inp.sel[ii], foo.3[iii] + 2 , foo.3[iii] + 2 ) == " ") {inp.sel[ii] <- paste(substr(inp.sel[ii],1,foo.3[iii]), "NA", substring(inp.sel[ii],foo.3[iii]+3) , sep="")}}}}
                    if(length(foo.3)==0) {cat(paste("There seem to be no values in any columns behind 'ESTIMATE'. Check outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep=""))}}
                 inp.sel    <- strsplit(inp.sel," +")
                 if(length(inp.sel[[1]]) == 0 ) {cat(paste("There seem to be no valid values associated with term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep=""))
                                                   all.terms <- all.terms[-i]}
                 if(length(inp.sel[[1]]) > 0 ) {
                    referenzlaenge <- max (sapply( inp.sel, FUN=function(ii ){  length(ii)    }) )
                    if(referenzlaenge < length(namen) ) {
                       cat(paste("Several columns seem to be empty for term '",all.terms[length(all.terms)],"' in file: '",file,"'.\n",sep=""))
                       head <- eatTools::crop(input.all[bereich[1]-2])
                       leerz<- gregexpr(" ", head)[[1]]
                       leerd<- which ( diff ( leerz) > 1 )[2]
                       vgl  <- length(strsplit ( eatTools::crop(substr(input.all[bereich[1]], 1, leerd)), split = " +")[[1]])
                       if ( vgl == 4 ) {
                            namen <- c(namen[1:2], "add.column1", namen[3:(referenzlaenge-1)])
                       }  else  {
                            referenzlaenge <- length(namen)
                       }
                    }
                    if(referenzlaenge > length(namen) ) {
                       if(referenzlaenge == length(namen) + 1) {
                          cat(paste("There seem to be one more column than columns names. Expect missing column name before 'ESTIMATE'. \nCheck outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep=""))
                          ind.name <- which(namen == "ESTIMATE")
                          namen    <- c(namen[1:ind.name-1], "add.column",namen[ind.name:length(namen)])}
                       if(referenzlaenge >  length(namen) + 1) {
                          cat(paste("There seem to be more columns than names for it. Check outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep=""))
                          namen<- c(namen, rep("add.column",referenzlaenge-length(namen) )) }}
                    inp.sel  <- t(sapply(inp.sel, FUN=function(ii){ c(ii, rep(NA,referenzlaenge-length(ii))) }))
                    colnames(inp.sel) <- namen
                    inp.sel  <- suppressWarnings(eatTools::asNumericIfPossible(data.frame( gsub("NNA", NA, gsub("NA", NA, gsub("\\*","",inp.sel))), stringsAsFactors = FALSE), force.string = FALSE))
                    results.sel<- data.frame(inp.sel,filename=as.character(file),stringsAsFactors = FALSE)
                    if(is.na(as.numeric(results.sel$ESTIMATE[1]))) {cat(paste("'ESTIMATE' column in Outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"' does not seem to be a numeric value. Please check!\n",sep=""))}
                    if(!missing(dif.term)) {
                       if(all.terms[length(all.terms)] == dif.term) {
                          cat(paste("Treat '",all.terms[length(all.terms)],"' as DIF TERM.\n",sep=""))
                          results.sel <- data.frame(results.sel,abs.dif = 2*results.sel$ESTIMATE,stringsAsFactors=FALSE)
                          konfNiveau  <- round(100*p.value)
                          results.sel[,paste("KI.",konfNiveau,".u",sep="")] <- results.sel$abs.dif-2*abs(qnorm(0.5*(1-p.value)))*results.sel$ERROR
                          results.sel[,paste("KI.",konfNiveau,".o",sep="")] <- results.sel$abs.dif+2*abs(qnorm(0.5*(1-p.value)))*results.sel$ERROR
                          results.sel[,paste("sig.",konfNiveau,sep="")] <- ifelse(abs(results.sel[,"abs.dif"])>abs.dif.bound & abs(results.sel[,paste("KI.",konfNiveau,".u",sep="")])>sig.dif.bound & abs(results.sel[,paste("KI.",konfNiveau,".o",sep="")])>sig.dif.bound,1,0)
                          results.sel$filename <- file
                          if(split.dif==TRUE) {results.sel <- results.sel[1:(dim(results.sel)[1]/2),]
                                               if(dim(results.sel)[1]!=dim(results.sel)[1]) {warning("missing variables in DIF table.")}}}}
                 all.output[[i]] <- results.sel}}
              if(!missing(dif.term)) {if(sum(all.terms==dif.term)==0) {cat(paste("Term declarated as DIF: '",dif.term,"' was not found in file: '",file,"'. \n",sep=""))  }}
              names(all.output) <- all.terms
            	regrStart <- grep("REGRESSION COEFFICIENTS", input.all) + 2
              isRegression <- length(regrStart) > 0
            	if ( isRegression)   {
                  regrEnd <- grep("An asterisk next", input.all)
              		regrEnd <- regrEnd[which(regrEnd > regrStart)][1] - 2
              		if ( is.na(regrEnd)) {
              		     warning(paste0("Regression coefficients seems to be missing or corrupted in '",file,"'."))
              		} else {
                    	 regrInput <- eatTools::crop(input.all[regrStart:regrEnd])
                  		 zeileDimensions <- grep("Regression Variable",input.all)
                       stopifnot(length(zeileDimensions) ==1)
                  		 nameDimensions  <- unlist(strsplit(input.all[zeileDimensions], "  +"))[-1]
                  		 regrRows <- grep("CONSTANT",input.all)
                       regrRows <- regrRows[regrRows<=regrEnd][1]
                  		 regrNamen <- unlist(lapply(strsplit(input.all[regrRows:regrEnd],"  +"), FUN=function(ii) {unlist(ii)[1]} ))
                       regrInputSel <- eatTools::crop(input.all[regrRows:regrEnd])
                       regrInputSel <- gsub("\\(","",regrInputSel)
                  		 regrInputSel <- gsub(")","",regrInputSel)
                  		 regrInputSel <- gsub("\\*","  NA",regrInputSel)
                  		 regrInputSel <- unlist( strsplit(regrInputSel," +") )
                  		 nDimensions  <- (length(  regrInputSel ) / length(regrNamen) - 1 )/2
                       cat(paste("Found ",nDimensions," dimension(s): ",paste(nameDimensions,collapse=", "),"\n",sep=""))
                       cat(paste("Found ",length(regrNamen)-1," regressor(s).\n",sep=""))
                       regrInputSel <- data.frame(matrix(regrInputSel, ncol=2*nDimensions+1, byrow=T),stringsAsFactors=F)
                       for (ii in 2:ncol(regrInputSel))  {regrInputSel[,ii] <- as.numeric(regrInputSel[,ii])}
                       colnames(regrInputSel) <- c("reg.var", paste(rep(c("coef","error"),nDimensions), rep(nameDimensions,each=2),sep="_") )
                       regrInputSel$filename <- file
                  		 all.output$regression <- regrInputSel
              		}
              }
              korStart <- grep("COVARIANCE/CORRELATION MATRIX", input.all)
              korEnd   <- grep("An asterisk next", input.all)
              korEnd   <- min(korEnd[korEnd > korStart])
              korStriche <- grep("-----",input.all)
              korStriche <- korStriche[korStriche > korStart & korStriche < korEnd]
              if(length(korStriche) == 2) {
                 varRow    <- grep("Variance", input.all)
                 variance  <- as.numeric( unlist( lapply(strsplit(input.all[varRow]," +"), FUN=function(ll) {ll[length(ll)]}) ) )
                 names(variance) <- "variance"
                 all.output$cov.structure <- variance
              }
              if(length(korStriche) > 2) {
                 bereich     <- input.all[ (min(korStriche) + 1) : (max(korStriche) - 1 ) ]
                 bereich     <- bereich[ -grep("----",bereich)]
                 bereich     <- strsplit(eatTools::crop(bereich),"  +")
                 for (ii in 2:(length(bereich)-1) )  {
                     if(ii <= length(bereich[[ii]]) )  {
                        bereich[[ii]] <- c(bereich[[ii]][1:(ii-1)], NA, bereich[[ii]][ii:length(bereich[[ii]])])
                     }
                     if(ii > length(bereich[[ii]]) )  {
                        bereich[[ii]] <- c(bereich[[ii]][1:(ii-1)], NA)
                     }
                 }
                 bereich.data.frame <- suppressWarnings(eatTools::asNumericIfPossible(data.frame(do.call("rbind", bereich[-1]),stringsAsFactors=FALSE), force.string = FALSE))
                 colnames(bereich.data.frame) <- bereich[[1]]
                 all.output$cov.structure <- bereich.data.frame
              }
            all.output$final.deviance <- rowToFind
            i1   <- grep("Dimension: \\(Dimension", input.all)
            all.output$reliability <- do.call("rbind", lapply(i1, FUN = function (z ) {
                    stopifnot(substr(input.all[z+3], 2,35) == "WLE Person separation RELIABILITY:")
                    stopifnot(substr(input.all[z+4], 2,20) == "EAP/PV RELIABILITY:")
                    return(data.frame ( dim = eatTools::crop(eatTools::crop(substring(input.all[z], 13)), ")"),
                           wle.rel = as.numeric(eatTools::crop(substring(input.all[z+3], 36))),
                           eap.rel = as.numeric(eatTools::crop(substring(input.all[z+4], 36))),
                           stringsAsFactors = FALSE))}))
            return(all.output)}

get.prm <- function(file)   {
            input <- scan(file,what="character",sep="\n",quiet=TRUE)
            input <- strsplit( gsub("\\\t"," ",eatTools::crop(input)), "/\\*")
            ret   <- data.frame ( do.call("rbind", strsplit( eatTools::crop(unlist(lapply(input, FUN = function ( l ) {l[1]}))), " +")), stringsAsFactors = FALSE)
            nameI <- eatTools::crop(eatTools::removePattern ( eatTools::crop( eatTools::crop(unlist(lapply(input, FUN = function ( l ) {l[length(l)]}))), char = "item"), pattern = "\\*/"))
            ret   <- data.frame ( Case= as.numeric(ret[,1]), item = nameI, parameter= as.numeric(ret[,2]) ,stringsAsFactors = FALSE)
            return(ret)}

get.itn <- function(file)  {
            input <- scan(file, what = "character", sep="\n", quiet = TRUE)
            ind.1 <- grep("==========",input)
            items <- grep( "item:", input )
            diff.last <- ind.1[length(ind.1)-1] - items[length(items)] + 4
            items <- cbind(1:length(items),items,c(diff(items),diff.last))
            ind.2 <- gregexpr(":", input[items[,2]])
            ind.3 <- unlist(ind.2)
            ind.3 <- matrix(ind.3,length(ind.2),byrow=T)
            item.namen <- substr(input[items[,2]], ind.3[,dim(ind.3)[2]]+1+nchar(as.character(items[,1])),100)
            item.namen <- gsub(" ","",item.namen)
            item.namen <- gsub("\\)","",item.namen); item.namen <- gsub("\\(","",item.namen)
            if(dim(ind.3)[2]>1)
              {stopifnot(length(table(ind.3[,1]))==1)
               dif.name <- rep(substr(input[items[,2]], 1, ind.3[,1]-1),(items[,3]-11))
               dif.value <- rep(as.numeric(substr(input[items[,2]], ind.3[,1]+1, ind.3[,1]+1)),(items[,3]-11))}
            zeilen <- list(); reihe <- NULL
            for (i in 1:dim(items)[1])
                {zeilen[[i]] <- (items[i,2]+7) : (items[i,2]+ (items[i,3]-5) )
                 cases       <- gsub("NA ","NA",input[zeilen[[i]]])
                 cases <- gsub("_BIG_ ","NA",cases)
                 cases <- gsub("_BIG_","NA",cases)
                 if(length(table(sapply(1:length(cases),FUN=function(ii){length(unlist(strsplit(cases[ii]," +"))) }) ) )>1 )
                   {cases <- gsub("          ","    NA    ",cases)}
                 cases       <- data.frame( matrix ( unlist( strsplit(eatTools::crop(gsub(" +"," ", cases))," ") ), nrow=length(zeilen[[i]]),byrow=T ) , stringsAsFactors=F)
                 ind         <- grep("\\)",cases[1,]); cases[,ind] <- gsub("\\)","",cases[,ind] )
                 cases       <- data.frame(cases[,1:(ind-1)],matrix(unlist(strsplit(cases[,6],"\\(")),nrow=length(zeilen[[i]]),byrow=T),cases[,-c(1:ind)],stringsAsFactors=F)
                 for(jj in 1:ncol(cases)) {cases[,jj] <- as.numeric(cases[,jj])}
                 colnames(cases) <- c("Label","Score","Abs.Freq","Rel.Freq","pt.bis","t.value","p.value",paste(rep(c("PV1.Avg.","PV1.SD."),((ncol(cases)-7)/2) ),rep(1:((ncol(cases)-7)/2),each=2),sep=""))
                 threshold.zeile   <- input[items[i,2]+2]; threshold <- NULL; delta <- NULL
                 bereich <- ifelse( (items[i,3]-12)<1,1,(items[i,3]-12))
                 if((items[i,3]-12)<1) {cat(paste("Item",i,"hat nur eine Antwortkategorie.\n"))}
                 for (j in 1: bereich )
                     {threshold  <- c(threshold ,as.numeric(substr(threshold.zeile,  6*j+16,6*j+21)))
                      delta      <- c(delta,     as.numeric(substr(input[items[i,2]+3],6*j+13,6*j+18)))}
                 while(length(threshold) < nrow(cases)) {threshold <- c(threshold,NA)}
                 while(length(delta) < nrow(cases)) {delta <- c(delta,NA)}
                 item.p <- NA
                 valid.p <- which(is.na(cases$Score))
                 if(length(valid.p) == 0)
                    {item.p <- cases[which(cases$Score == max(cases$Score)),"Abs.Freq"] / sum(cases$Abs.Freq)}
                 sub.reihe   <- data.frame(item.nr=i, item.name=item.namen[i], cases[,1:2], n.valid = sum(cases$Abs.Freq), cases[,3:4], item.p = item.p, diskrim=as.numeric(substr(input[items[i,2]+1],45,55)),cases[,-c(1:4)], threshold, delta, stringsAsFactors=F)
                 reihe <- rbind(reihe,sub.reihe)}
             if(dim(ind.3)[2]>1)
               {reihe <- data.frame(dif.name,dif.value,reihe,stringsAsFactors=FALSE)}
             return(reihe)}

get.dsc <- function(file) {
            input     <- scan(file,what="character",sep="\n",quiet=TRUE)
            n.gruppen    <- grep("Group: ",input)
            gruppennamen <- unlist( lapply( strsplit(input[n.gruppen]," ") , function(ll) {paste(ll[-1],collapse=" ")} ) )
            cat(paste("Found ",length(n.gruppen)," group(s) in ",file,".\n",sep=""))
            trenner.1 <- grep("------------------",input)
            trenner.2 <- grep("\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.",input)
            stopifnot(length(trenner.1) == length(trenner.2))
            daten     <- lapply(1:(length(trenner.1)/2), FUN=function(ii) {
                 dat <- strsplit(input[(trenner.1[2*ii]+1):(trenner.2[2*ii-1]-1)]," +")
                 dat <- data.frame(matrix(unlist(lapply(dat, FUN=function(iii) {  c(paste(iii[1:(length(iii)-4)],collapse=" "),iii[-c(1:(length(iii)-4))])  })), ncol=5,byrow=T) , stringsAsFactors=F)
                 dat <- data.frame(group.name = gruppennamen[ii], dat, stringsAsFactors = FALSE)
                 colnames(dat) <- c("group.name","dimension","N","mean","std.dev","variance")
                 for (iii in 3:ncol(dat)) {dat[,iii] <- as.numeric(dat[,iii])}
                 desc <- strsplit(input[(trenner.2[2*ii-1]+1):(trenner.2[2*ii]-1)]," +")
                 desc <- data.frame(matrix(unlist(lapply(desc, FUN=function(iii) {  c(paste(iii[1:(length(iii)-3)],collapse=" "),iii[-c(1:(length(iii)-3))])  })), ncol=4,byrow=T) , stringsAsFactors=F)
                 colnames(desc) <- c("dimension","mean","std.dev","variance")
                 for (iii in 2:ncol(desc)) {desc[,iii] <- as.numeric(desc[,iii])}
                 dat.list <- list( single.values=dat, aggregates=desc)
                 return(dat.list) } )
            names(daten) <- gruppennamen
            n.dim        <- names(table(unlist(lapply(1:length(daten), FUN=function(ii) {length( grep("Error", daten[[ii]]$aggregates$dimension))}) ) ))
            stopifnot(length(n.dim)==1)
            cat(paste("Found ",n.dim," dimension(s) in ",file,".\n",sep=""))
            return(daten)}


get.equ <- function(file)  {
            input       <- scan(file,what="character",sep="\n",quiet = TRUE)
            dimensionen <- grep("Equivalence Table for",input)
            cat(paste("Found ",length(dimensionen), " dimension(s).\n",sep=""))
            ende        <- grep("================",input)
            ende        <- sapply(dimensionen, FUN=function(ii) {ende[ende>ii][1]})
            tabellen    <- lapply(1:length(dimensionen), FUN=function(ii)
                           {part <- eatTools::crop(input[(dimensionen[ii]+6):(ende[ii]-1)])
                            part <- data.frame(matrix(as.numeric(unlist(strsplit(part," +"))),ncol=3,byrow=T),stringsAsFactors=F)
                            colnames(part) <- c("Score","Estimate","std.error")
                            return(part)})
            regr.model  <- grep("The regression model",input)
            item.model  <- grep("The item model",input)
            stopifnot(length(regr.model) == length(item.model))
            name.dimensionen <- unlist( lapply(dimensionen,FUN=function(ii) {unlist(lapply(strsplit(input[ii], "\\(|)"),FUN=function(iii){iii[length(iii)]}))}) )
            model       <- lapply(1:length(regr.model), FUN=function(ii) {rbind ( eatTools::crop(gsub("The regression model:","",input[regr.model[ii]])), eatTools::crop(gsub("The item model:","",input[item.model[ii]])) ) })
            model       <- do.call("data.frame",args=list(model,row.names=c("regression.model","item.model"),stringsAsFactors=F))
            colnames(model) <- name.dimensionen
            tabellen$model.specs <- model
            names(tabellen)[1:length(dimensionen)] <- name.dimensionen
            return(tabellen)}






isLetter <- function ( string ) {
            splt <- strsplit(string, "")
            isL  <- lapply(splt, FUN = function ( x ) {
                    ind <- which ( x %in% c( letters , LETTERS ))
                    x[setdiff(1:length(x),ind)] <- " "
                    x <- eatTools::crop(paste(x, sep="", collapse=""))
                    x <- unlist ( strsplit(x, " +") )
                    return(x)  } )
            return(isL)}


getConquestVersion <- function ( path.conquest , path.temp , asDate = TRUE ) {
    wd <- path.temp
		f <- file.path ( wd , "delete.cqc" )
		write ( "quit;" , f )
		f <- normalizePath ( f )
		path.conquest <- normalizePath ( path.conquest )
		cmd <- paste ( "\"", path.conquest, "\" \"", f , "\"" , sep ="")
		r <- NULL
		suppressWarnings(try ( r <- system ( command = cmd , intern = TRUE ) , silent = TRUE ))
		file.remove ( f )
		if ( !is.null ( r ) ) {
				r <- r[1]
				r <- sub ( "ConQuest build: " , "" , r )
				r <- gsub ( "\\s+" , "-" , r )
				if ( asDate ) r <- date::as.date(r)
		}
		return (r)}





prepRep <- function ( calibT2, bistaTransfT1, bistaTransfT2, makeIdsUnique = TRUE) {
           if ( !inherits(calibT2, "transfBista" )) { stop("'calibT2' object must be of class 'transfBista'.\n")}
           if ( !inherits(bistaTransfT1, "transfBista" )) { stop("'bistaTransfT2' object must be of class 'transfBista'.\n")}
           if ( !inherits(bistaTransfT2, "transfBista") ) { stop("'bistaTransfT2' object must be of class 'transfBista'.\n")}
           if (!nrow(calibT2[["itempars"]]) < nrow(bistaTransfT1[["itempars"]])) { stop("Mismatch between 'calibT2' and 'bistaTransfT1'. \n")}
           if (!nrow(calibT2[["itempars"]]) < nrow(bistaTransfT2[["itempars"]])) { stop("Mismatch between 'calibT2' and 'bistaTransfT2'. \n")}
           idT1<- unique(bistaTransfT1[["all.Names"]][which(bistaTransfT1[["all.Names"]][,"par"] == "ID"),"derived.par"])
           idT2<- unique(bistaTransfT2[["all.Names"]][which(bistaTransfT2[["all.Names"]][,"par"] == "ID"),"derived.par"])
           stopifnot(length(idT1)==1, length(idT2)==1)
           if ( idT1 != idT2 ) {
                warning(paste0("ID variables do not match between t1 and t2. ID for t1: '",idT1,"'. ID for t2: '",idT2,"'. \n    IDs will be unified with '",idT1,"'."))
                recStat <- paste ( "'", idT2 , "' = '", idT1, "'", sep="")
                colnames ( bistaTransfT2[["personpars"]] ) <- car::recode ( colnames ( bistaTransfT2[["personpars"]] ), recStat)
           }
           lc  <- colnames( calibT2[["personpars"]] ) [grep("^linking", colnames(calibT2[["personpars"]]) )]
           if(length(lc)==0) { stop("No columns with linking error information found in 'calibT2'.\n")}
           lcn <- paste("trend", eatTools::removePattern(string = lc, pattern = "linking"), sep="")
           colnames( calibT2[["personpars"]] ) [grep("^linking", colnames(calibT2[["personpars"]]) )] <- lcn
           merg<- c("group", "imp", "traitLevel", "dimension")
           frms<- list ( calibT2=calibT2, bistaTransfT1=bistaTransfT1, bistaTransfT2=bistaTransfT2 )
           toM <- unique(unlist(lapply ( names(frms), FUN = function ( l.Name ) {
                  l    <- frms[[l.Name]]
                  drin <- merg %in% colnames(l[["personpars"]])
                  fehlt<- merg[which(drin==FALSE)]
                  if (!all(drin == TRUE)) { warning(paste0("Column(s) '",paste(fehlt, collapse = "', '"), "' are unexpectedly missing in '",l.Name,"'."))}
                  keep <- merg[which(drin==TRUE)]
                  return(keep)})))
           if(length(toM)==0) { stop("Merging impossible.\n")}
           red <- calibT2[["personpars"]][,c(toM,  lcn)]
           red <- red[!duplicated(red),]
           dat1<- data.frame ( trend = "T1" , merge ( bistaTransfT1[["personpars"]], red, by = toM, all = TRUE))
           stopifnot ( nrow(dat1) == nrow(bistaTransfT1[["personpars"]]))
           dat2<- data.frame ( trend = "T2" , merge ( bistaTransfT2[["personpars"]], red, by = toM, all = TRUE))
           stopifnot ( nrow(dat2) == nrow(bistaTransfT2[["personpars"]]))
           if ( makeIdsUnique == TRUE ) {
                dat1[, paste(idT1, "unique", sep="_")] <- paste(dat1[, "trend"], dat1[, idT1], sep="_")
                dat2[, paste(idT1, "unique", sep="_")] <- paste(dat2[, "trend"], dat2[, idT1], sep="_")
           }
           return(rbind ( dat1, dat2))}

plotICC <- function ( resultsObj, defineModelObj, item = NULL, personPar = c("WLE", "EAP", "PV"), personsPerGroup = 30, pdfFolder = NULL, smooth = 7 ) {
           personPar  <- match.arg(arg = toupper(personPar), choices = c("WLE", "EAP", "PV"))
           if (smooth<5) {smooth <- 5}
           it  <- itemFromRes ( resultsObj )
           if ( !"est" %in% colnames(it) ) { it[,"est"] <- NA }
           if ( !"estOffset" %in% colnames(it) ) { it[,"estOffset"] <- NA }
           it[,"est"] <- rowSums(it[,c("est", "estOffset")], na.rm = TRUE)
           if ( !"estSlope" %in% colnames(it) ) { it[,"estSlope"] <- 1 }
           if ( length(which(is.na(it[,"estSlope"]))) > 0) { it[which(is.na(it[,"estSlope"])), "estSlope"] <- 1 }
           eapA<- eapFromRes (resultsObj)
           if ( personPar == "WLE") {
                eapA <- wleFromRes(resultsObj)
                colnames(eapA) <- car::recode(colnames(eapA), "'wle_est'='EAP'")
           }
           if ( personPar == "PV") {
                eapA <- pvFromRes(resultsObj, toWideFormat = TRUE)
                colnames(eapA) <- car::recode(colnames(eapA), "'pv1'='EAP'")
           }
           cat("Note: To date, only 1pl/2pl dichotomous models are supported.\n"); flush.console()
           if ( is.null(item) & is.null(pdfFolder)) {stop("If ICCs for more than one item should be displayed, please specify an output folder for pdf.\n")}
           if ( !is.null(pdfFolder)) { grDevices::pdf(file = pdfFolder, width = 10, height = 7.5) }
           if ( !is.null ( item ) )  {
                if ( !item %in% it[,"item"]) { stop (paste("Item '",item,"' was not found in 'resultsObj'.\n",sep=""))}
                it <- it[which(it[,"item"] == item),]
           }
           pl  <- by ( data = it, INDICES = it[,c("model", "item")], FUN = function ( i ) {
                  xlm <- c(i[["est"]]+2, i[["est"]]-2)
                  anf <- -6
                  ende<- 6
                  x   <- seq ( anf, ende, l = 400)
                  y   <- exp( i[["estSlope"]]*x - i[["est"]] ) / (1+exp( i[["estSlope"]]*x - i[["est"]] ))
                  plot (x, y, type = "l", main = paste("Item '",as.character(i[["item"]]),"'\n\n",sep=""), xlim = c(-6,6), ylim = c(0,1), xlab = "theta", ylab = "P(X=1)", col = "darkred", cex = 8, lwd = 2)
                  graphics::mtext( paste("Model = ",i[["model"]],"  |  Dimension = ",i[["dimension"]], "  |  difficulty = ",round(i[["est"]], digits = 3),"  |  Infit = ",round(i[["infit"]], digits = 3),"\n",sep=""))
                  eap <- eapA[intersect ( which (eapA[,"dimension"] == i[["dimension"]]) , which (eapA[,"model"] == i[["model"]])),]
                  if ( inherits(defineModelObj, "defineMultiple")) {
                       woIst<- which ( lapply ( defineModelObj, FUN = function ( g ) {   g[["analysis.name"]] == i[["model"]] }) == TRUE)
                       stopifnot(length(woIst) == 1)
                       dat  <-defineModelObj[[woIst]][["daten"]]
                  }  else  {
                       dat  <- defineModelObj[["daten"]]
                  }
                  id  <- unique(resultsObj[intersect(which(resultsObj[,"type"] == "tech"), which(resultsObj[,"par"] == "ID")),"derived.par"])
                  stopifnot(length(id)==1)
                  prbs<- na.omit ( merge ( dat[,c( "ID", as.character(i[["item"]]))], eap[,c( id, "EAP")], by.x = "ID", by.y = id))
                  anz <- round ( nrow(prbs) / personsPerGroup ) + 1
                  if ( anz < 3 ) { anz <- 3 }
                  if ( anz > smooth) { anz <- round(smooth)}
                  eapQ<- quantile ( prbs[,"EAP"], probs = seq(0,1,l = anz))
                  prbs[,"gr"] <- eatTools::num.to.cat ( x = prbs[,"EAP"], cut.points = eapQ[-c(1,length(eapQ))])
                  prbs<- do.call("rbind", by ( data = prbs, INDICES = prbs[,"gr"], FUN = function ( g ) {
                         g[,"mw"] <- mean(g[,"EAP"])
                         g[,"anz"]<- length(g[,"EAP"])
                         g[,"lh"] <- mean(g[, as.character(i[["item"]]) ])
                         return(g)}))
                  matr<- prbs[!duplicated(prbs[,c("mw", "lh")]),c("mw", "lh")]
                  matr<- data.frame(matr[sort(matr[,"mw"],decreasing=FALSE,index.return=TRUE)$ix,])
                  graphics::points ( x = matr[,"mw"], y = matr[,"lh"], cex = 1, pch = 21, bg = "darkblue")
                  graphics::lines ( x = matr[,"mw"], y = matr[,"lh"], col = "blue", lty = 3, lwd = 3) } )
           if ( !is.null(pdfFolder)) { grDevices::dev.off() } }




