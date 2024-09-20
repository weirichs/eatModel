



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




