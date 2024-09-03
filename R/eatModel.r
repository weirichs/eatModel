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
      B.obj  <- array(unlist(lapply(2:ncol(qMatrix),FUN = function ( col) {data.frame ( Cat0 = 0, Cat1 = qMatrix[,col])})), dim = c(nrow(qMatrix), 2, ncol(qMatrix)-1), dimnames = list(qMatrix[,"item"], c("Cat0", "Cat1"), paste0("Dim0",1:(ncol(qMatrix)-1)) ))
      tamObj <- list ( AXsi = xsi.obj, B = B.obj, resp = resp, Y=Y, pid = pid)
      class(tamObj) <- c("list", "tamBayes")
      return(tamObj)}


nObsItemPairs <- function ( responseMatrix, q3MinType) {
                 spl <- data.frame ( combinat::combn(colnames(responseMatrix),2), stringsAsFactors = FALSE)
                 splM<- do.call("rbind", lapply ( spl, FUN = function ( y ) {
                        if ( q3MinType == "singleObs" ) {
                             minVal <- min ( table(data.frame ( responseMatrix[,y])))
                        }  else  {
                             minVal <- min(c(rowSums(table(data.frame ( responseMatrix[,y]))), colSums(table(data.frame ( responseMatrix[,y])))))
                        }
                        ret <- data.frame ( Var1 = sort(y)[1], Var2 = sort(y)[2], minValue = minVal)
                        return(ret)}))
                 return(splM)}


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


getResults <- function ( runModelObj, overwrite = FALSE, Q3 = TRUE, q3theta = c("pv", "wle", "eap"), q3MinObs = 0, q3MinType = c("singleObs", "marginalSum"), omitFit = FALSE, omitRegr = FALSE, omitWle = FALSE, omitPV = FALSE, abs.dif.bound = 0.6, sig.dif.bound = 0.3, p.value = 0.9,
              nplausible = NULL, ntheta = 2000, normal.approx = FALSE, samp.regr = FALSE, theta.model=FALSE, np.adj=8, group = NULL, beta_groups = TRUE, level = .95, n.iter = 1000, n.burnin = 500, adj_MH = .5, adj_change_MH = .05, refresh_MH = 50, accrate_bound_MH = c(.45, .55),	sample_integers=FALSE, theta_init=NULL, print_iter = 20, verbose = TRUE, calc_ic=TRUE, omitUntil=1, seed=NA) {
            q3MinType<- match.arg(q3MinType)
            q3theta  <- match.arg(q3theta )
            if(inherits(runModelObj, "runMultiple")) {                          
                if(is.null ( attr(runModelObj, "split")[["nCores"]] ) || attr(runModelObj, "split")[["nCores"]] == 1 ) {
                   res <- lapply( runModelObj, FUN = function ( r ) {           
                          do  <- paste ( "getResults ( ", paste(names(formals(getResults)), car::recode(names(formals(getResults)), "'runModelObj'='r'"), sep =" = ", collapse = ", "), ")",sep="")
                          ret <- eval(parse(text=do))
                          return(ret)})
                   }  else  {
                          doIt<- function (laufnummer,  ... ) {
                                 if(!exists("getResults"))  { library(eatModel) }
                                 if(!exists("tam.mml") &  length(grep("tam.", class(runModelObj[[1]])))>0 ) {library(TAM, quietly = TRUE)}
                                 do  <- paste ( "getResults ( ", paste(names(formals(getResults)), car::recode(names(formals(getResults)), "'runModelObj'='runModelObj[[laufnummer]]'"), sep =" = ", collapse = ", "), ")",sep="")
                                 ret <- eval(parse(text=do))
                                 return(ret)}
                          beg <- Sys.time()
                          if ( attr(runModelObj, "split")[["mcPackage"]] == "parallel") {
                               cl  <- makeCluster(attr(runModelObj, "split")[["nCores"]], type = "SOCK")
                          }  else  {
                               cl  <- future::makeClusterPSOCK(attr(runModelObj, "split")[["nCores"]], verbose=FALSE)
                          }
                          res <- clusterApply(cl = cl, x = 1:length(runModelObj), fun = doIt , overwrite = overwrite, omitFit = omitFit, omitRegr = omitRegr, omitWle = omitWle, omitPV = omitPV, abs.dif.bound = abs.dif.bound, sig.dif.bound = sig.dif.bound, p.value = p.value)
                          stopCluster(cl)
                          cat(paste ( "Results of ",length(runModelObj), " analyses processed: ", sep="")); print( Sys.time() - beg, digits=3)
                   }
               res <- do.call("rbind", res )
               class(res) <- c("data.frame", "multipleResults")
               rownames(res) <- NULL
               return(res)
            }  else {                                                           
               if ( is.null(runModelObj)) {return(NULL)}
               isTa  <- FALSE
               if(inherits(runModelObj, "runConquest")) {                       
                    if ( isTRUE(Q3) ) {
                        if ( ncol ( runModelObj[["qMatrix"]]) !=2 ) {
                            cat("Q3 is only available for unidimensional models. Estimation will be skipped.\n")
                            Q3 <- FALSE
                        }
                    }
                    do    <- paste ( "res <- getConquestResults ( ", paste(names(formals(getConquestResults)), car::recode(names(formals(getConquestResults)), "'path'='runModelObj$dir'; 'analysis.name'='runModelObj$analysis.name'; 'model.name'='runModelObj$model.name'; 'qMatrix'='runModelObj$qMatrix'; 'all.Names'='runModelObj$all.Names'; 'deskRes'='runModelObj$deskRes'; 'discrim'='runModelObj$discrim'; 'daten'='runModelObj$daten'"), sep =" = ", collapse = ", "), ")",sep="")
                    eval(parse(text=do))                                        
                    dir <- runModelObj[["dir"]]                                 
                    name<- runModelObj[["analysis.name"]]                       
                    allN<- runModelObj[["all.Names"]]                           
               }  else  {                                                       
                    isTa<- TRUE                                                 
                    if ( Q3 ) {
                        if ( ncol ( attr(runModelObj, "defineModelObj")[["qMatrix"]]) !=2 ) {
                            cat("Q3 is only available for unidimensional models. Estimation will be skipped.\n")
                            Q3 <- FALSE
                        }
                    }
                    if(!is.null(nplausible)) { attr(runModelObj, "defineModelObj")[["n.plausible"]] <- nplausible }  else  { nplausible <- attr(runModelObj, "defineModelObj")[["n.plausible"]] }
                    do    <- paste ( "res <- getTamResults ( ", paste(names(formals(getTamResults)), car::recode(names(formals(getTamResults)),"'pvMethod'='attr(runModelObj, \"defineModelObj\")[[\"pvMethod\"]]'"),  sep =" = ", collapse = ", "), ")",sep="")
                    eval(parse(text=do))
                    dir <- attr(runModelObj, "defineModelObj")[["dir"]]         
                    name<- attr(runModelObj, "defineModelObj")[["analysis.name"]]
                    allN<- attr(runModelObj, "defineModelObj")[["all.Names"]]   
               }
               if(!is.null(res)) {                                              
                    stopifnot ( length(unique(res[,"model"])) == 1)             
                    alln<- do.call("rbind", lapply(names(allN), FUN = function ( x ) {
                           if ( length( allN[[x]] ) > 0 ) {
                                res <- data.frame ( type = "tech", par = x, derived.par = allN[[x]])
                           }  else  {
                                res <- NULL
                           }
                           return(res)}))
                    res <- plyr::rbind.fill ( res, data.frame ( res[1,c("model", "source")], alln, stringsAsFactors = FALSE) )
                    difS<- list (abs.dif.bound = abs.dif.bound, sig.dif.bound = sig.dif.bound, p.value = p.value)
                    resD<- data.frame ( res[1,c("model", "source")], type = "tech", par = "dif", derived.par = names(difS), value = unlist(difS), stringsAsFactors = FALSE)
                    if ( inherits(runModelObj, "runConquest")) {                
                         resN<- data.frame ( res[1,c("model", "source")], type = "tech", par = c("method",rep("nodes", 3)), derived.par = c(runModelObj[["method"]],"nodes", "p.nodes", "f.nodes"), value = c(1,runModelObj[["nodes"]], runModelObj[["p.nodes"]], runModelObj[["f.nodes"]]), stringsAsFactors = FALSE)
                    }  else  {                                                  
                         nNod<- length(attr(runModelObj, "defineModelObj")[["control"]][["nodes"]])
                         resN<- data.frame ( res[1,c("model", "source")], type = "tech", par = c("QMC",rep("nodes", 1+nNod)), derived.par = c(NA,rep("discrete.theta",nNod), "snodes"), value = c(attr(runModelObj, "defineModelObj")[["control"]][["QMC"]],attr(runModelObj, "defineModelObj")[["control"]][["nodes"]], attr(runModelObj, "defineModelObj")[["control"]][["snodes"]]), stringsAsFactors = FALSE)
                    }
                    res <- plyr::rbind.fill ( res, resD, resN )
                    id  <- unique(res[intersect(which(res[,"type"] == "tech"), which(res[,"par"] == "ID")),"derived.par"])
                   if(!is.null(dir)) {
                        stopifnot(length(id)==1)
                        item<-itemFromRes ( res )
                        if ( file.exists(file.path(dir, paste(name, "_items.csv",sep=""))) & overwrite == FALSE) {
                             cat(paste("Item results cannot be saved, file '",  file.path(dir, paste(name, "_items.csv",sep="")),"' already exists.\n    Please remove/rename existing file or use 'overwrite=TRUE'.\n",sep=""))
                        }  else  {
                             write.csv2(item, file.path(dir, paste(name, "_items.csv",sep="")), na="", row.names = FALSE)
                        }                                                       
                        txt <- capture.output ( wle <- wleFromRes(res) )        
                        if (!is.null ( wle ) ) {
                             wleL<- reshape2::melt ( wle, id.vars = c(id, "dimension"), measure.vars = c("wle_est", "wle_se"), na.rm = TRUE)
                             form<- as.formula ( paste ( id, "~dimension+variable",sep=""))
                             wleW<- reshape2::dcast ( wleL, form, value.var = "value" )
                        }
                        txt <- capture.output ( pv  <- pvFromRes(res) )
                        if(!is.null(pv)) {
                             pvL <- reshape2::melt ( pv, id.vars = c( id , "dimension"), na.rm = TRUE)
                             form<- as.formula ( paste ( id, "~dimension+variable",sep=""))
                             pvW <- reshape2::dcast ( pvL, form, value.var = "value" )
                        }
                        txt <- capture.output ( eap <- eapFromRes(res) )
                        if(!is.null(eap)) {
                             eapL<- reshape2::melt ( eap, id.vars = c(id, "dimension"), measure.vars = c("EAP", "SE.EAP"), na.rm = TRUE)
                             form<- as.formula ( paste ( id, "~dimension+variable",sep=""))
                             eapW<- reshape2::dcast ( eapL, form, value.var = "value" )
                        }                                                       
                        alls<- list ( wle, pv, eap )                            
                        allP<- NULL                                             
                        notN<- which ( unlist(lapply ( alls, FUN = function ( x ) { !is.null(x)})) )
                        if ( length( notN ) >= 1 ) { allP <- alls[[notN[1]]] }
                        if ( length( notN ) > 1 )  {
                             for ( u in notN[-1] )   {
                                   allP <- merge ( allP, alls[[u]], by = c ( id, "dimension"), all = TRUE)
                             }
                        }
                        if ( !is.null(allP)) {
                              if ( file.exists(file.path(dir, paste(name, "_persons.csv",sep=""))) & overwrite == FALSE) {
                                   cat(paste("Person estimates cannot be saved, file '",  file.path(dir, paste(name, "_persons.csv",sep="")),"' already exists.\n    Please remove/rename existing file or use 'overwrite=TRUE'.\n",sep=""))
                              }  else  {
                                   write.csv2(allP, file.path(dir, paste(name, "_persons.csv",sep="")), na="", row.names = FALSE)
                              }
                        }
                        if ( Q3 == TRUE ) {
                              q3m <- q3FromRes ( res )                          
                              stopifnot(length(q3m)==1)
                              q3m <- q3m[[1]]
                              if ( file.exists(file.path(dir, paste(name, "_q3.csv",sep=""))) & overwrite == FALSE) {
                                   cat(paste("Item results cannot be saved, file '",  file.path(dir, paste(name, "_q3.csv",sep="")),"' already exists.\n    Please remove/rename existing file or use 'overwrite=TRUE'.\n",sep=""))
                              }  else  {
                                   write.csv2(q3m, file.path(dir, paste(name, "_q3.csv",sep="")), na="", row.names = FALSE)
                              }
                        }
                   }
               }
               rownames(res) <- NULL
               return(res)
               }}

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

equat1pl<- function ( results , prmNorm , item = NULL, domain = NULL, testlet = NULL, value = NULL, excludeLinkingDif = TRUE, difBound = 1, iterativ = FALSE, method = c("Mean.Mean", "Haebara", "Stocking.Lord", "robust", "Haberman"),
           itemF = NULL, domainF = NULL, testletF = NULL, valueF = NULL, estimation=c("OLS", "BSQ", "HUB", "MED", "LTS", "L1", "L0"), b_trim=Inf, lts_prop=.5) {
           estimation <- match.arg(estimation)
           method     <- match.arg(method)
           isRunM<- all(c("model" , "source" , "var1" , "var2" , "type" , "indicator.group", "group", "par", "derived.par", "value") %in% names(results))
           if ( isRunM) {
                nMods <- table(results[,"model"])
                cat(paste("Found ", length(nMods), " model(s).\n   Equating is executed for each dimension in each model separately.\n",sep=""))
                dims  <- unique(unlist(by ( data = results, INDICES = results[,"model"], FUN = function ( x ) { names(table(as.character(itemFromRes(x)[,"dimension"])))})))
                if ( is.null(dims)) {
                     dims <- unique(na.omit(results[,"group"]))
                     warning(paste0("Cannot extract dimensions from 'results' object. This should only occur for bayesian plausible values imputation. Assume following dimensions: \n    '",paste(dims, collapse = "', '"),"'."))
                }
           }  else  {
                resList <- transformItemParListIntoResults (results = results, itemF = itemF, domainF = domainF, testletF = testletF, valueF = valueF)
                results <- resList[["results"]]
                dims    <- nMods <- resList[["dims"]]
           }
           if ( missing ( prmNorm) ) {                                          
                if ( isFALSE(isRunM) ) { stop("No norm parameter defined ('prmNorm' is missing).\n")}
                cat("No norm parameter defined ('prmNorm' is missing). Treat current sample as drawn from the reference population.\n")
                items <- by ( data = results, INDICES = results[,"model"], FUN = function ( d ) {
                         dimN <- buildEmptyResultsObject(d=d, method = method, results=results)
                         return(dimN)}, simplify = FALSE)
                ret   <- list(items = items, results = results)                 
                class(ret) <- c("eq2tom", class(ret))                           
                return(ret)
           }  else {                                                            
                prmNorm<- eatTools::makeDataFrame(prmNorm)
                allN   <- checkItemParLists(prmNorm =prmNorm, item = item, domain = domain, testlet = testlet, value = value, dims=dims)
                items <- by ( data = results, INDICES = results[,"model"], FUN = function ( d ) {
                       if ( isRunM  ) {
                            it  <- itemFromRes(d)
                       }  else  {
                            it  <- d
                       }
                       if ( "estOffset" %in% colnames ( it ) ) {
                            cat(paste("W A R N I N G:  Model '",d[1,"model"],"' was estimated with (at least partially) anchored items parameters. Equating seems questionable.\n",sep=""))
                            d[,"par"] <- car::recode ( d[,"par"], "'offset'='est'")
                            it <- itemFromRes(d)
                       }
                       dimN <- by ( data = it, INDICES = it[,"dimension"], FUN = function ( prmDim ) {
                               if(!is.null(allN[["domain"]]) ) {
                                   prmM<- prmNorm [ which(prmNorm[,allN[["domain"]]] %in% unique(it[,"dimension"])) ,]
                               }  else  {
                                   prmM<- prmNorm
                               }
                               mess1 <- NULL                                    
                               if (!is.null(allN[["testlet"]])) {
                                   mrge <- merge(prmDim[,c("item", "dimension")], prmM, by="item", all=FALSE)
                                   stopifnot(nrow(mrge)>0)
                                   if(length(which(is.na(mrge[,allN[["testlet"]]]))) > 0) {
                                       mess1 <- paste0("Domain '",prmDim[1,"dimension"],"': Found ",length(which(is.na(mrge[,allN[["testlet"]]]))), " missing values in '",allN[["testlet"]],"' column of 'prmNorm'. Withdraw from incorporating testlets into linking error computation.")
                                       allN[["testlet"]] <- NULL                
                                   }
                               }
                               if ( length(prmDim[, "item"]) != length(unique(prmDim[, "item"])) ) {  stop(paste("Items are not unique for model '",as.character(d[1,"model"]),"'.\n",sep="")) }
                               eq  <- equAux ( x = prmDim[ ,c("item", "est")], y = prmM[,c(allN[["item"]], allN[["value"]], allN[["testlet"]])] )
                               if ( eq[["descriptives"]][["N.Items"]] > 0) {
                                     if ( method == "robust") {
                                         prm<- merge(prmM[,c(allN[["item"]], allN[["value"]], allN[["testlet"]])], prmDim[ ,c("item", "est")], by.y="item", by.x = allN[["item"]], all=FALSE)
                                         eqr<- sirt::linking.robust(prm)        
                                     }                                          
                                     if ( method == "Haberman") {               
                                         prm<- rbind ( data.frame ( study = "norm", item = prmM[,allN[["item"]]], a=1, b = prmM[,allN[["value"]]], stringsAsFactors = FALSE), data.frame ( study = "focus", item = prmDim[ ,"item"], a=1, b = prmDim[ ,"est"], stringsAsFactors = FALSE))
                                         eqh<- sirt::linking.haberman(prm, progress = FALSE, estimation = estimation, b_trim=b_trim, lts_prop=lts_prop)
                                     }
                               }
                               if ( method %in% c("Mean.Mean", "Haebara", "Stocking.Lord")) {
                                     dif <- eq[["anchor"]][,"TransfItempar.Gr1"] - eq[["anchor"]][,"Itempar.Gr2"]
                                     prbl<- which ( abs ( dif ) > difBound )
                               }  else  {
                                   if ( eq[["descriptives"]][["N.Items"]] == 0) { eqr <- eqh <- NULL}
                               }
                               foo <- printToConsole(d=d, nMods=nMods, it=it, prmDim=prmDim, eq=eq, allN=allN, method=method, estimation=estimation, eqh=eqh, eqr=eqr, mess1=mess1)
                               if ( method != "robust" && method != "Haberman" && length( prbl ) > 0 ) {
                                    eld  <- handleLinkingDif(prmDim=prmDim,prbl=prbl, eq=eq, difBound=difBound, dif=dif, method=method, excludeLinkingDif=excludeLinkingDif, iterativ=iterativ,prmM=prmM, allN=allN)
                               }  else  {                                       
                                    eld  <- noLinkingDif(method=method, eq=eq, eqr=eqr, eqh=eqh)
                               }
                               if( isFALSE(iterativ) && excludeLinkingDif && !is.null(eld[["info2"]])) {
                                   cat("\nItems with DIF:\n")
                                   print(eatTools::roundDF(eld[["info2"]][,1:2], digits = 3)); flush.console()
                               }
                               cat("\n")
                               print(eatTools::roundDF(eld[["info"]])); flush.console()
                               cat("\n")
                               if(!is.null(mess1)) {message(mess1); cat("\n")}
                               if ( method %in% c("robust", "Haberman")) {
                                    eld <- createOutput(method=method, eqr=eqr, prm=prm, eqh=eqh, info=eld[["info"]])
                               }
                               ret <- list ( eq = eld[["eq"]], items = prmDim, info = eld[["info"]], method = method )
                               return ( ret ) }, simplify =FALSE)
                       return(dimN) }, simplify = FALSE)
              ret  <- list ( items = items, results = results)                  
              class(ret) <- c("eq2tom", class(ret))                             
              return(ret)                                                       
              }  }

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

transformToBista <- function ( equatingList, refPop, cuts, weights = NULL, defaultM = 500, defaultSD = 100, q3bound = .20, roman = FALSE, vera = TRUE, idVarName = NULL, years = NULL ) {
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
       isRunM<- all(c("model" , "source" , "var1" , "var2" , "type" , "indicator.group", "group", "par", "derived.par", "value") %in% names(equatingList[["results"]]))
       if (isRunM) {
           it     <- itemFromRes(equatingList[["results"]])
       }  else  {
           it     <- equatingList[["results"]]
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
       if (isRunM) {
           id     <- unique(equatingList[["results"]][intersect(which(equatingList[["results"]][,"type"] == "tech"), which(equatingList[["results"]][,"par"] == "ID")),"derived.par"])
           if(length(id)!=1) {
               id   <- getIdVarName(id=NULL, idVarName, verbose=TRUE)
           }
           refList<- lapply ( dims, FUN = function (dimname) {
                     rex  <- pvFromRes(equatingList[["results"]][unique(c(which(equatingList[["results"]][,"group"] == dimname),which(equatingList[["results"]][,"type"] == "tech"))), ], toWideFormat = FALSE, idVarName=idVarName, verbose=FALSE)
                     if (is.null(rex)) {return(NULL)}                               
                     if ( is.null(weights) ) {
                          txt <- capture.output ( msd <- eatRep::repMean ( datL = rex, ID = id, imp = "imp", dependent = "value", na.rm = TRUE))
                          msd <- adaptEatRepVersion(msd)                            
                     }  else  {                                                     
                        rex <- eatTools::mergeAttr ( rex, weights , by.x = id, by.y = colnames(weights)[1], all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "cases", xName = paste0("plausible values for dimension ",dimname), yName = "weights", verbose = c("match", "dataframe"))
                        mis <- which(is.na(rex[,colnames(weights)[2]]))
                        if ( length(mis) > 0 ) {                                    
                             if(length(mis) == nrow(rex)) {stop(paste("Mergin of weights and plausible values for '", dimname, "' failed. No common units."))}
                             cat(paste ( "Found ",length(mis)," missing values in the 'weights' frame.\n    Cases with missing values on weighting variable will be ignored for transformation.\n",sep=""))
                             rex <- rex[-mis,]
                        }
                        txt <- capture.output ( msd <- eatRep::repMean ( datL = rex, ID = id, imp = "imp", wgt = colnames(weights)[2], dependent = "value", na.rm = TRUE) )
                        msd <- adaptEatRepVersion(msd)
                     }
                     rp <- data.frame ( domain = dimname , m = msd[intersect(which(msd[,"parameter"] == "mean"), which(msd[,"coefficient"] == "est")),"value"], sd = msd[intersect(which(msd[,"parameter"] == "sd"), which(msd[,"coefficient"] == "est")),"value"])
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
                   refPop[unique(mat[mis,"x"]), unique(mat[mis,"y"])] <- ref[unique(mat[mis,"x"]), unique(mat[mis,"y"])]
              }
           }
       }
       if(ncol ( refPop ) == 3) {
           cat ( paste("The 'refPop' data.frame does not include information about reference population mean/SD on Bista metric. Values will be defaulted to ",defaultM,"/",defaultSD,".\n",sep=""))
           refPop[,4] <- defaultM; refPop[,5] <- defaultSD
       }  else  {
           if ( ncol ( refPop) != 5 ) { stop ( "Invalid 'refPop'.\n") }
       }
       if(!all(dims %in% refPop[,1]))  {warning(paste0("Dimension names in the 'results' object '",paste(dims, collapse="', '"), "' do not match to names in the first columns of 'refPop': '",paste(refPop[,1], collapse="', '"),"."))}
       if(!all(dims %in% names(cuts))) {warning(paste0("Dimension names in the 'results' object '",paste(dims, collapse="', '"), "' do not match to names in the cut scores object: '",paste(names(cuts), collapse="', '"),"."))}
       modN<- lapply(nam1, FUN = function ( mod ) {                             
              nam2 <- names(equatingList[["items"]][[mod]])
              dimN <- lapply(nam2, FUN = function ( dims ) {                    
                      if (isRunM) {
                          resMD<- equatingList[["results"]][unique(c(intersect(which(equatingList[["results"]][,"model"] == mod), which(equatingList[["results"]][,"group"] == dims)),  which(equatingList[["results"]][,"type"] == "tech"))),]
                          rex  <- pvFromRes(resMD, toWideFormat = TRUE, idVarName = idVarName, verbose=FALSE)
                          if (!is.null(rex)) {
                              if ( length ( rex[,id]) != unique(length ( rex[,id])) ) {
                                   stop(paste( "Model '",mod,"', Dimension '",dims,"': cases according to '", id,"' variable are not unique.\n",sep=""))
                              }
                          }
                          offSet  <- grep("offset", as.character(resMD[,"par"]))
                          if(length(offSet)>0) {  resMD[,"par"] <- car::recode ( resMD[,"par"], "'offset'='est'") }
                          itFrame <- itemFromRes(resMD)
                      }  else  {
                          itFrame <- it[intersect(which(it[,"model"] == mod), which(it[,"dimension"] == dims)),]
                      }
                      if ( !is.null(itFrame) && !itFrame[1,"dimension"] %in% refPop[,1] ) {
                            cat(paste("Cannot found dimension '",itFrame[1,"dimension"],"' in the first column of the 'refPop' argument. Skip transformation ... \n",sep=""))
                            return ( list ( itempars = NULL, personpars = NULL, rp = NULL))
                      }  else  {
                            if ( is.null ( itFrame ) ) {
                                cat(paste0("Model '",mod,"', dimension '",dims,"': No item parameters found. This should only occur for bayesian plausible values imputation. Transformation of item parameters will be skipped.\n"))
                            }  else  {
                                if ( isFALSE(cutsMis) ) {
                                     if ( !itFrame[1,"dimension"] %in% names(cuts) ) { stop(paste("Cannot found dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.",sep=""))}
                                     mat1<- match( itFrame[1,"dimension"], names(cuts))
                                     if ( !"values" %in% names(cuts[[mat1]]) ) { stop(paste("'cuts' must be a named list. Cannot found 'values' element for dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
                                     if ( length(cuts[[mat1]])>1) {
                                         if ( !"labels" %in% names(cuts[[mat1]]) ) { stop(paste("'cuts' must be a named list. Cannot found 'labels' element for dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
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
                                         txt <- capture.output ( msdF <- eatRep::repMean ( datL = pv, ID = id, imp = "imp", dependent = "valueTransfBista", na.rm = TRUE))
                                         msdF<- adaptEatRepVersion(msdF)
                                    }  else  {
                                         pvF <- eatTools::mergeAttr ( pv, weights , by.x = id, by.y = colnames(weights)[1], all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "cases", xName = paste0("plausible values for dimension ",dims), yName = "weights", verbose = c("match", "dataframe"))
                                         mis <- which(is.na(pvF[,colnames(weights)[2]]))
                                         if ( length(mis) > 0 ) {                       
                                              cat(paste ( "Found ",length(mis)," missing values in the 'weights' frame.\n    Cases with missing values on weighting variable will be ignored for transformation.\n",sep=""))
                                              pvF <- pvF[-mis,]
                                         }
                                         txt <- capture.output ( msdF <- eatRep::repMean ( datL = pvF, ID = id, imp = "imp", wgt = colnames(weights)[2], dependent = "valueTransfBista", na.rm = TRUE) )
                                         msdF<- adaptEatRepVersion(msdF)
                                    }
                                    msdFok <- c(msdF[intersect(which(msdF[,"parameter"] == "mean"), which(msdF[,"coefficient"] == "est")),"value"], msdF[intersect(which(msdF[,"parameter"] == "sd"), which(msdF[,"coefficient"] == "est")),"value"])
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
                       			                    del<- data.frame ( traitLevel = attr(traitLevel, "cat.values")[l], linkingErrorTraitLevel = del )
                       			                    return(del)}))
                                         ori <- colnames(itFrame)
                                         chk <- unique(le[,"traitLevel"]) %in% unique(itFrame[,"traitLevel"])
                                         if ( length( which(chk == FALSE)) > 0) {
                                             warning(paste("Model '",unique(itFrame[,"model"]),"', dimension '",unique(itFrame[,"dimension"]),"': No items on trait level(s) '",paste( unique(le[,"traitLevel"])[which(chk == FALSE)], collapse = "', '"), "'.", sep=""))
                                         }
                                         itFrame <- eatTools::mergeAttr ( itFrame, le, by = "traitLevel", sort = FALSE, all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "trait levels", xName = "item parameter list", yName = "linking error list", verbose = c("match"))
                                         itFrame <- itFrame[,c(ori, "linkingErrorTraitLevel")]
                                    }
                                    itFrame[,"refMean"]        <- refPop[mat,2]
                                    itFrame["refSD"]           <- refPop[mat,3]
                                    itFrame[,"refTransfMean"]  <- refPop[mat,4]
                                    itFrame[,"refTransfSD"]    <- refPop[mat,5]
                                }
                                if(!exists("mat1") ) {mat1 <- match(dims, names(cuts)); stopifnot(length(mat1)==1)}
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
                                        pv[,"linkingError"] <- equatingList[["items"]][[mod]][[dims]][["eq"]][["descriptives"]][["linkerror"]]
                                        pv[,"linkingErrorTransfBista"] <- unique ( itFrame[,"linkingErrorTransfBista"])
                                    }
                                    ori <- colnames(pv)                             
                                    if ( cutsMis == FALSE && !is.null ( equatingList[["items"]]) && exists("le") ) {
                                         pv  <- eatTools::mergeAttr ( pv, le, by = "traitLevel", sort = FALSE, all.x = TRUE, all.y = FALSE, setAttr = FALSE, unitName = "trait levels", xName = "plausible values", yName = "linking error list", verbose = c("match"))
                                         pv  <- pv[,c(ori, "linkingErrorTraitLevel")]
                                    }
                                    if (!is.null(weights)) {
                                         pv  <- merge ( pv, weights , by.x = id, by.y = colnames(weights)[1], all.x = TRUE, all.y = FALSE)
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
           itemVera <- createItemVeraObj(itempars=itempars, roman=roman, results = equatingList[["results"]], q3bound=q3bound)
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
       ret        <- list ( itempars = itempars, personpars = personpars, refPop = refPop, means = rp, all.Names = context, itemparsVera = itemVera, linkingErrors = leo)
       class(ret) <- c("list", "transfBista")
       return( ret ) }


runModel <- function(defineModelObj, show.output.on.console = FALSE, show.dos.console = TRUE, wait = TRUE) {
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
                   return(mod)  }  }   }

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


defineModel <- function(dat, items, id, splittedModels = NULL, irtmodel = c("1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", "2PL.groups", "GPCM.design", "3PL"),
               qMatrix=NULL, DIF.var=NULL, HG.var=NULL, group.var=NULL, weight.var=NULL, anchor = NULL, domainCol=NULL, itemCol=NULL, valueCol=NULL,check.for.linking = TRUE,
               minNperItem = 50, removeMinNperItem = FALSE, boundary = 6, remove.boundary = FALSE, remove.no.answers = TRUE, remove.no.answersHG = TRUE, remove.missing.items = TRUE, remove.constant.items = TRUE,
               remove.failures = FALSE, remove.vars.DIF.missing = TRUE, remove.vars.DIF.constant = TRUE, verbose=TRUE, software = c("conquest","tam"), dir = NULL,
               analysis.name, schooltype.var = NULL, model.statement = "item",  compute.fit = TRUE, pvMethod = c("regular", "bayesian"), fitTamMmlForBayesian = TRUE, n.plausible=5, seed = NULL, conquest.folder=NULL,
               constraints=c("cases","none","items"),std.err=c("quick","full","none"), distribution=c("normal","discrete"), method=c("gauss", "quadrature", "montecarlo", "quasiMontecarlo"),
               n.iterations=2000,nodes=NULL, p.nodes=2000, f.nodes=2000,converge=0.001,deviancechange=0.0001, equivalence.table=c("wle","mle","NULL"), use.letters=FALSE,
               allowAllScoresEverywhere = TRUE, guessMat = NULL, est.slopegroups = NULL, fixSlopeMat = NULL, slopeMatDomainCol=NULL, slopeMatItemCol=NULL, slopeMatValueCol=NULL,
               progress = FALSE, Msteps = NULL, increment.factor=1 , fac.oldxsi=0, export = list(logfile = TRUE, systemfile = FALSE, history = TRUE, covariance = TRUE, reg_coefficients = TRUE, designmatrix = FALSE) )   {
                  software <- match.arg(arg = software, choices = c("conquest","tam"))
                  if ( software == "conquest" && is.null(conquest.folder) ) {   
                       conquest.folder <- identifyConquestFolder()
                  }
                  ismc <- attr(dat, "multicore")                                
                  dat  <- eatTools::makeDataFrame(dat, name = "dat")
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
                     irtmodel <- match.arg(irtmodel)
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
                            qMatrix <- qMatrix[-match(notInDat, qMatrix[,1]),]
                            if(nrow(qMatrix) == 0) { stop("No common items in Q matrix and data.\n")}
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
                          return ( ret )    }   }  }

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

checkContextVars <- function(x, varname, type = c("weight", "DIF", "group", "HG"), itemdata, suppressAbort = FALSE, internal = FALSE)   {
                     type <- match.arg(arg = type, choices = c("weight", "DIF", "group", "HG"))
                     stopifnot(length(x) == nrow(itemdata))
                     if(missing(varname))  {varname <- "ohne Namen"}            
                     if(!inherits(x, c("numeric", "integer")) && isTRUE(internal))  {
                        if (type == "weight") {stop(paste(type, " variable has to be 'numeric' necessarily. Automatic transformation is not recommended. Please transform by yourself.\n",sep=""))}
                        cat(paste(type, " variable has to be 'numeric'. Variable '",varname,"' of class '",class(x),"' will be transformed to 'numeric'.\n",sep=""))
                        x <- suppressWarnings(unlist(eatTools::asNumericIfPossible(x = data.frame(x, stringsAsFactors = FALSE), transform.factors = TRUE, maintain.factor.scores = FALSE, force.string = FALSE)))
                        if(!inherits(x, "numeric"))  {                          
                           x <- as.numeric(as.factor(x))
                        }
                        cat(paste("    '", varname, "' was converted into numeric variable of ",length(table(x))," categories. Please check whether this was intended.\n",sep=""))
                        if(length(table(x)) < 12 ) { cat(paste("    Values of '", varname, "' are: ",paste(names(table(x)), collapse = ", "),"\n",sep=""))}
                     }
                     toRemove<- NULL
                     mis     <- length(unique(x))
                     if(mis == 0 )  {
                        if ( suppressAbort == FALSE ) {
                             stop(paste("Error: ",type," Variable '",varname,"' without any values.",sep=""))
                        }  else  {
                             cat(paste0("Warning: ", type," Variable '",varname,"' without any values. '",varname,"' will be removed.\n"))
                             toRemove <- varname
                        }
                     }
                     if(mis == 1 )  {
                        if ( suppressAbort == FALSE ) {
                             stop(paste("Error: ",type," Variable '",varname,"' is a constant.",sep=""))
                        }  else  {
                             cat(paste0(type," Variable '",varname,"' is a constant. '",varname,"' will be removed.\n"))
                             toRemove <- varname
                        }
                     }
                     if(type == "DIF" | type == "group") {
                        if(mis > 10 && isTRUE(internal))   {warning(paste0(type," Variable '",varname,"' with more than 10 categories. Recommend recoding."))}
                     }
                     wegDifMis <- NULL; wegDifConst <- NULL; char <- 1; weg <- which(is.na(1:12)); info <- NULL
                     if ( is.null(toRemove)) {
                          char    <- max(nchar(as.character(na.omit(x))))
                          weg     <- which(is.na(x))
                          if(length(weg) > 0 ) {warning(paste0("Found ",length(weg)," cases with missing on ",type," variable '",varname,"'. Conquest probably will collapse unless cases are not deleted.\n"))}
                          if(type == "DIF" ) {
                                        if(mis > 2 && isTRUE(internal))   {cat(paste(type, " Variable '",varname,"' does not seem to be dichotomous.\n",sep=""))}
                                        y       <- paste0("V", x)               
                                        n.werte <- lapply(itemdata, FUN=function(iii){by(iii, INDICES=list(y), FUN=table, simplify=FALSE)})
                                        completeMissingGroupwise <- data.frame(t(sapply(n.werte, function(ll){lapply(ll, FUN = function (uu) { length(uu[uu>0])}  )})), stringsAsFactors = FALSE)
                                        for (iii in seq(along=completeMissingGroupwise)) {
                                             missingCat.i <- which(completeMissingGroupwise[,iii] == 0)
                                             if(length(missingCat.i) > 0) {
                                                cat(paste("Warning: Following ",length(missingCat.i)," items with no values in ",type," variable '",varname,"', group ",substring(colnames(completeMissingGroupwise)[iii],2),": \n",sep=""))
                                                wegDifMis <- c(wegDifMis, rownames(completeMissingGroupwise)[missingCat.i] )
                                                cat(paste0("   ", paste(rownames(completeMissingGroupwise)[missingCat.i],collapse=", "), "\n"))
                                                info      <- plyr::rbind.fill(info, data.frame ( varname = varname, varlevel = substring(colnames(completeMissingGroupwise)[iii],2), nCases = table(y)[colnames(completeMissingGroupwise)[iii]], type = "missing", vars =rownames(completeMissingGroupwise)[missingCat.i], stringsAsFactors = FALSE))
                                             }
                                             constantCat.i<- which(completeMissingGroupwise[,iii] == 1)
                                             if(length(constantCat.i) > 0) {
                                                cat(paste("Warning: Following ",length(constantCat.i)," items are constants in ",type," variable '",varname,"', group ",substring(colnames(completeMissingGroupwise)[iii],2),":\n",sep=""))
                                                wegDifConst<- c(wegDifConst, rownames(completeMissingGroupwise)[constantCat.i] )
                                                values    <- n.werte[rownames(completeMissingGroupwise)[constantCat.i]]
                                                values    <- lapply(values, FUN = function(v){v[[colnames(completeMissingGroupwise)[iii]]]})
                                                cat(paste0("   ", paste(rownames(completeMissingGroupwise)[constantCat.i],collapse=", "), "\n"))
                                                info      <- plyr::rbind.fill(info, data.frame ( varname = varname, varlevel = substring(colnames(completeMissingGroupwise)[iii],2), nCases = table(y)[colnames(completeMissingGroupwise)[iii]], type = "constant", vars =names(values), value =  sapply(values, names), nValue = unlist(values), stringsAsFactors = FALSE))
                                             }
                                        }
                          }
                     }
                     return(list(x = x, char = char, weg = weg, varname=varname, wegDifMis = wegDifMis, wegDifConst = wegDifConst, toRemove = toRemove, info=info))}


prepAnchorTAM <- function (ank, allNam) {
        if(!is.null(ank)) {
            stopifnot(ncol(ank) == 2 )                                          
            notInData   <- setdiff(ank[,1], allNam[["variablen"]])              
            if(length(notInData)>0)  {ank <- ank[-match(notInData, ank[,1]),]}
            ank[,1]    <- match(as.character(ank[,1]), allNam[["variablen"]])
        }
        return(ank)}

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
            
checkItemConsistency <- function(dat, allNam, remove.missing.items, verbose, removeMinNperItem, minNperItem, remove.constant.items, model.statement){
          namen.items.weg <- NULL                                               
          is.NaN <- do.call("cbind", lapply(dat[,allNam[["variablen"]], drop = FALSE], FUN = function (uu) { is.nan(uu) } ) )
          if(sum(is.NaN) > 0 ) {
             cat(paste("Found ",sum(is.NaN)," 'NaN' values in the data. Convert 'NaN' to 'NA'.\n",sep=""))
             for ( j in allNam[["variablen"]]) {
                   weg <- which ( is.nan(dat[,j] ))
                   if(length(weg)>0) {  dat[weg,j] <- NA }
             }
          }
          n.werte <- eatTools::tableUnlist(dat[,allNam[["variablen"]], drop = FALSE])
          zahl    <- grep("[[:digit:]]", names(n.werte))                        ### sind das alles Ziffern? (auch wenn die Spalten als "character" klassifiziert sind)
          noZahl  <- setdiff(1:length(n.werte), zahl)
          if (length( zahl ) == 0 )  { stop("Please use numeric values for item responses.\n")}
          if (length( noZahl ) > 0 ) { cat(paste(" W A R N I N G !  Found ",sum(n.werte[noZahl])," non-numeric values in the item responses. These values will be treated as missing responses!\n",sep="")) }
          klasse  <- unlist( lapply(dat[,allNam[["variablen"]], drop = FALSE], class) )
          if(any(unlist(lapply(dat[,allNam[["variablen"]], drop = FALSE], inherits, what=c("integer", "numeric"))) == FALSE)) {
               cat(paste(" W A R N I N G !  Found unexpected class type(s) in item response columns: '",paste(setdiff(klasse, c("numeric", "integer")), collapse = "', '"), "'\n",sep=""))
               cat("                  All item columns will be transformed to be 'numeric'. Recommend to edit your data manually prior to analysis.\n")
               for ( uu in allNam[["variablen"]] ) { dat[,uu] <- as.numeric(as.character(dat[,uu]))}
          }
          values  <- lapply(dat[,allNam[["variablen"]], drop = FALSE], FUN = function ( ii ) { table(ii)})
          isDichot<- unlist(lapply(values, FUN = function ( vv ) { identical(c("0","1"), names(vv)) }))
          n.werte <- sapply(values, FUN=function(ii) {length(ii)})
          n.mis   <- which(n.werte == 0)
          if(length(n.mis) >0) {
             cat(paste("Serious warning: ",length(n.mis)," testitems(s) without any values.\n",sep=""))
             if(verbose == TRUE) {cat(paste("    ", paste(names(n.mis), collapse=", "), "\n", sep=""))}
             if(remove.missing.items == TRUE) {
                 cat(paste("Remove ",length(n.mis)," variable(s) due to solely missing values.\n",sep=""))
                 namen.items.weg <- c(namen.items.weg, names(n.mis))
             }
          }
          if ( removeMinNperItem == TRUE ) {                                    
               nValid <- unlist(lapply(dat[,allNam[["variablen"]], drop = FALSE], FUN = function ( ii ) { length(na.omit ( ii )) }))
               below  <- which ( nValid < minNperItem )
               if ( length ( below ) > 0 ) {
                    cat (paste ( "Found ", length(below), " items with less than ", minNperItem, " valid responses. These items will be removed.\n", sep=""))
                    namen.items.weg <- unique ( c(namen.items.weg, names(below)))
               }
          }
          constant <- which(n.werte == 1)
          if(length(constant) >0) {
             cat(paste("Warning: ",length(constant)," testitems(s) are constants.\n",sep=""))
             if(verbose == TRUE) {foo <- lapply(names(constant),FUN=function(ii) {cat(paste(ii,": ",names(table(dat[,ii])), " ... ",length(na.omit(dat[,ii]))," valid responses", sep="")); cat("\n")})}
             if(remove.constant.items == TRUE) {
                 cat(paste("Remove ",length(constant)," variable(s) due to solely constant values.\n",sep=""))
                 namen.items.weg <- c(namen.items.weg, names(constant))
             }
          }
          n.rasch  <- which( !isDichot )                                        
          if(length(n.rasch) >0 )   {                                           
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
          return(list(dat=dat,allNam=allNam, namen.items.weg=namen.items.weg))}
          
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

checkQmatrixConsistency <-  function(qmat) {
             qmat  <- eatTools::makeDataFrame(qmat, name = "Q matrix")
             if(!inherits(qmat[,1], "character")) { qmat[,1] <- as.character(qmat[,1])}
             nClass<- sapply(qmat[,-1,drop=FALSE], inherits, what=c("numeric", "integer"))
             if ( !all(nClass)) {
                  warning(paste0("Found non-numeric indicator column(s) in the Q matrix. Transform column(s) '",paste(colnames(qmat)[ which(nClass==FALSE)+1], collapse = "', '") ,"' into numeric format."))
                  qmat <- data.frame ( qmat[,1,drop=FALSE], eatTools::asNumericIfPossible(qmat[,-1,drop=FALSE]), stringsAsFactors = FALSE)
             }
             werte <- eatTools::tableUnlist(qmat[,-1,drop=FALSE], useNA="always")
             if(length(setdiff( names(werte) , c("0","1", "NA")))<0) {stop("Q matrix must not contain entries except '0' and '1'.\n")}
             if(werte[match("NA", names(werte))] > 0) {stop("Missing values in Q matrix.\n")}
             wertes<- lapply(qmat[,-1,drop=FALSE], FUN = function (col)  {all ( col == 0)})
             konst <- which(wertes == TRUE)
             if ( length(konst)>0) {
                  cat(paste0("Column(s) '",paste(names(konst), collapse = "', '"),"' in Q matrix are konstant with value 0. Delete column(s).\n"))
                  qmat <- qmat[,-match(names(konst), colnames(qmat)), drop=FALSE]
             }
             doppel<- which(duplicated(qmat[,1]))
             if(length(doppel)>0) {
                cat("Found duplicated elements in the item id column of the q matrix. Duplicated elements will be removed.\n")
                chk  <- table(qmat[,1])
                chk  <- chk[which(chk > 1)]
                chkL <- lapply(names(chk), FUN = function ( ch ) {
                        qChk <- qmat[which(qmat[,1] == ch),]
                        pste <- apply(qChk, 1, FUN = function ( x ) { paste(x[-1], collapse="")})
                        if( !all ( pste == pste[1] )) { stop("Inconsistent q matrix.\n")}
                        })
                qmat <- qmat[!duplicated(qmat[,1]),]
             }
             zeilen<- apply(qmat, 1, FUN = function ( y ) { all ( names(table(y[-1])) == "0")  })
             weg   <- which(zeilen == TRUE)
             if(length(weg)>0) {
                cat(paste("Note: Following ",length(weg)," item(s) in Q matrix do not belong to any dimension. Delete these item(s) from Q matrix.\n",sep=""))
                cat("    "); cat(paste(qmat[weg,1],collapse=", ")); cat("\n")
                qmat  <- qmat[-weg,]
             }
             return(qmat)}

converged<- function (dir, logFile) {
            isConv <- TRUE
            if (!file.exists(file.path ( dir, logFile ))) {
                 warning(paste0("Model seems not to have converged. Cannot find log file '",file.path ( dir, logFile ),"'."))
                 isConv <- FALSE
            }  else  {
                 logF  <- scan(file = file.path ( dir, logFile ), what="character",sep="\n",quiet=TRUE)
                 if(length(logF) == 0 ) {
                    warning(paste0("Model seems not to have converged. Log file '",file.path ( dir, logFile ),"' is empty."))
                    isConv <- FALSE
                 }  else  {
                    last  <- logF[length(logF)]
                    if ( ! eatTools::crop(last) == "=>quit;" ) {
                       if ( length( grep("quit;" , last)) == 0 ) {
                           warning(paste0("Model seems not to have converged. Log file unexpectedly finishs with '",last,"'.\nReading in model output might fail."))
                           isConv <- FALSE
                       }  }  }  }
            return(isConv)  }

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

getConquestShw <- function (model.name, qMatrix, qL, shw, altN){
         shw1 <- data.frame ( model = model.name, source = "conquest", var1 = shw$item[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = shw$item[,"dimensionName"], par = "est",  derived.par = NA, value = as.numeric(shw$item[,"ESTIMATE"]), stringsAsFactors = FALSE)
         shw2 <- data.frame ( model = model.name, source = "conquest", var1 = shw$item[,"item"], var2 = NA , type = "fixed", indicator.group = "items",group = shw$item[,"dimensionName"], par = "est",  derived.par = "se", value = as.numeric(shw$item[,"ERROR"]), stringsAsFactors = FALSE)
         toOff<- shw2[ which(is.na(shw2[,"value"])), "var1"]
         if(length(toOff)>0) {
            shw1[match(toOff, shw1[,"var1"]), "par"] <- "offset"
            shw2  <- shw2[-which(is.na(shw2[,"value"])),]                       ### entferne Zeilen aus shw2, die in der "value"-Spalte NA haben
         }
         return(list(shw1=shw1, shw2=shw2))}
         
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

getConquestDiscrim <- function (model.name, discrim , qMatrix, qL){
         if( is.null(discrim) )  {return(NULL)}
         discR<- merge(discrim, qL[,-match("value", colnames(qL))], by.x = "item.name", by.y = colnames(qMatrix)[1], all=TRUE)
         shw5 <- data.frame ( model = model.name, source = "conquest", var1 = discR[,"item.name"], var2 = NA , type = "fixed", indicator.group = "items", group = discR[,"dimensionName"], par = "itemDiscrim",  derived.par = NA, value = discR[,"item.diskrim"], stringsAsFactors = FALSE)
         return(shw5)}

getConquestInfit <- function (model.name,  shw){
         res <- rbind(data.frame ( model = model.name, source = "conquest", var1 = shw[["item"]][,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = shw$item[,"dimensionName"], par = "est",  derived.par = "infit", value = as.numeric(shw$item[,"MNSQ.1"]), stringsAsFactors = FALSE),
                      data.frame ( model = model.name, source = "conquest", var1 = shw[["item"]][,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = shw$item[,"dimensionName"], par = "est",  derived.par = "outfit", value = as.numeric(shw$item[,"MNSQ"]), stringsAsFactors = FALSE) )
         return(res)}


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
                      warning(paste0("Cannot identify the group the term '",i,"' in file '",shwFile,"' belongs to. Insert 'NA' to the 'group' column."))
                      gr <- NA
                   }  else {
                      gr <- colnames(qMatrix)[2]
                   }
                   vars<- c("ESTIMATE", "MNSQ", "MNSQ.1", "ERROR")
                   cls <- sapply(shw[[i]][,vars], inherits, what=c("numeric", "integer"))
                   if ( !all(cls) ) {
                        warning(paste0("Expect column(s) '",paste(vars[which(cls==FALSE)],collapse= "', '"), "' in file '",shwFile,"' (statement '",i,"') to be numeric. Current column format is: '",paste(sapply(shw[[i]][,vars[which(cls==FALSE)]],class), collapse="', '"),"'. Column will be transformed."))
                        shw[[i]] <- eatTools::set.col.type(shw[[i]], col.type = list("numeric.if.possible" = names(cls[which(cls==FALSE)])), maintain.factor.scores = TRUE)
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
         if ( !is.null(q3MinObs) ) {                                            
              if ( q3MinObs > 1 ) { nObs <- nObsItemPairs ( responseMatrix = daten[,all.Names[["variablen"]]], q3MinType = q3MinType ) }
         }
         matL  <- reshapeQ3 (mat = q3.res$q3.matrix, q3MinObs = q3MinObs, nObs = nObs)
         if( nrow(matL)== 0) { return(NULL)}
         res   <- data.frame ( model = model.name, source = "conquest", var1 = matL[,"Var1"],  var2 = matL[,"Var2"] , type = "fixed",indicator.group = "items",group = paste(names(table(shw1[,"group"])), collapse="_"), par = "q3", derived.par = NA, value = matL[,"value"] , stringsAsFactors = FALSE)
         return(res)}

getConquestDeviance <- function ( path, analysis.name, omitUntil = omitUntil) {
         cqc  <- scan(file.path ( path, paste0(analysis.name, ".cqc")),what="character",sep="\n",quiet=TRUE)
         such <- c("method", "nodes")
         ret  <- lapply(such, FUN = function ( su ) {
                 indm <- grep(paste0(su, "="), cqc)
                 if ( length(indm)>1) {                                         
                    hf   <- grep("f_nodes", cqc)
                    indm <- setdiff(indm, hf)
                 }
                 if(length(indm) != 1) {
                    cat(paste("Cannot identify '",su,"'from cqc file.\n",sep=""))
                    met <- NULL
                 }  else  {
                    pos1<- nchar(unlist(strsplit(cqc[indm], su))[1])            
                    pos2<- which(sapply(1:nchar(cqc[indm]), FUN = function(x){ substr(cqc[indm],x,x) == ","}))
                    pos2<- min(pos2[which(pos2>pos1)])
                    met <- eatTools::removePattern(substr(cqc[indm], pos1+1, pos2-1), paste0(su,"="))
                 }
                 return(met)})                                                  
         tme  <- file.info ( file.path ( path, paste0(analysis.name, ".shw")))[["mtime"]] - file.info ( file.path ( path, paste0(analysis.name, ".cqc")))[["mtime"]]
         grDevices::pdf(file = file.path ( path, paste0(analysis.name, "_dev.pdf")), width = 10, height = 7.5)
         plotDevianceConquest ( logFile = list ( path=path, analysis.name=analysis.name, ret=ret, tme=tme), omitUntil = omitUntil)
         grDevices::dev.off() }


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

getTamItempars    <- function(runModelObj, qL, qMatrix, leseAlles) {
         if(leseAlles == FALSE) {return(NULL)}                                  
         if ( is.null(attr(runModelObj, "defineModelObj")[["all.Names"]][["DIF.var"]])) {
              xsis <- merge(data.frame ( item = rownames(runModelObj[["xsi"]]), runModelObj[["xsi"]], stringsAsFactors = FALSE), qL[,-match("value", colnames(qL))],  by.x = "item", by.y = colnames(qMatrix)[1], all = TRUE)
         }  else  {                                                             
              xsis <- mergeDimensionIfDIF (dat = data.frame ( item = rownames(runModelObj[["xsi"]]), runModelObj[["xsi"]], stringsAsFactors = FALSE), qmatLong = qL[,-match("value", colnames(qL))], datMergingVar="item", remove = "toMerge")
         }
         shw1 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = xsis[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = xsis[,"dimensionName"], par = "est",  derived.par = NA, value = xsis[,"xsi"], stringsAsFactors = FALSE)
         shw2 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = xsis[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = xsis[,"dimensionName"], par = "est",  derived.par = "se", value = xsis[,"se.xsi"], stringsAsFactors = FALSE)
         if ( !is.null(attr(runModelObj, "defineModelObj")[["all.Names"]][["DIF.var"]])) {          
               shw1 <- renameDifParameters (dat=shw1, qmatLong = qL[,-match("value", colnames(qL))])
               shw2 <- renameDifParameters (dat=shw2, qmatLong = qL[,-match("value", colnames(qL))])
         }
         toOff<- shw2[ which(shw2[,"value"] == 0 ), "var1"]                     
         if(length(toOff)>0) {
            shw1[match(toOff, shw1[,"var1"]), "par"] <- "offset"
            shw2  <- shw2[-which(shw2[,"value"] == 0 ),] }                      ### entferne Zeilen aus shw2, die in der "value"-Spalte NA haben, danach: p-Werte einfuegen
         return(list ( shw1=shw1, shw2=shw2))}

getTamDescriptives    <- function(runModelObj, qL, qMatrix, leseAlles) {
         if(leseAlles == FALSE || is.null ( attr(runModelObj, "defineModelObj")[["deskRes"]] )) {return(NULL)}
         deskR<- merge(attr(runModelObj, "defineModelObj")[["deskRes"]], qL[,-match("value", colnames(qL))],  by.x = "item.name", by.y = colnames(qMatrix)[1], all = TRUE)
         shw3 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = as.character(deskR[,"item.name"]), var2 = NA , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "itemP",  derived.par = NA, value = deskR[,"item.p"], stringsAsFactors = FALSE)
         shw4 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = as.character(deskR[,"item.name"]), var2 = NA , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "Nvalid",  derived.par = NA, value = deskR[,"valid"], stringsAsFactors = FALSE)
         cols <- setdiff ( colnames(deskR)[grep("^item.p", colnames(deskR))], "item.p")
         if ( length ( cols ) > 0 ) {
              colsR <- data.frame ( original = cols, reduziert = eatTools::removePattern ( string = cols, pattern = "item.p.") , stringsAsFactors = FALSE)
              shw31 <- do.call("rbind", apply ( colsR, MARGIN = 1, FUN = function ( zeile ) { data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = as.character(deskR[,"item.name"]), var2 = zeile[["reduziert"]] , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "itemP",  derived.par = NA, value = deskR[,zeile[["original"]]], stringsAsFactors = FALSE) }))
              return(rbind(shw3, shw4, shw31))
         }  else  {
              return(rbind(shw3, shw4))
         } }

getTamDiscrim    <- function(runModelObj, qL, qMatrix, leseAlles) {
         if(leseAlles == FALSE || is.null ( attr(runModelObj, "defineModelObj")[["discrim"]] )) {return(NULL)}
         discR<- merge( attr(runModelObj, "defineModelObj")[["discrim"]] , qL[,-match("value", colnames(qL))],  by.x = "item.name", by.y = colnames(qMatrix)[1], all = TRUE)
         shw5 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = discR[,"item.name"], var2 = NA , type = "fixed", indicator.group = "items", group = discR[,"dimensionName"], par = "itemDiscrim",  derived.par = NA, value = discR[,"item.diskrim"], stringsAsFactors = FALSE)
         return(shw5)}

getTam2plDiscrim <- function(runModelObj, qMatrix, leseAlles, regr, omitRegr) {
         if(leseAlles == FALSE || !attr(runModelObj, "defineModelObj")[["irtmodel"]] %in% c("2PL", "2PL.groups", "GPCM", "3PL") ) {return(NULL)}
         shw6 <- do.call("rbind", lapply (  1 : length ( colnames( qMatrix ) [-1] ) , FUN = function ( dims ) {
                 if ( isFALSE(omitRegr) ) {                                     
                      obj <- regr[["B"]]                                        
                 } else {                                                       
                      obj <- as.data.frame ( runModelObj[["B"]])                
                      colnames(obj) <- paste0("B.", gsub("Dim0", "Dim", colnames(obj)))
                      obj[,"item"]  <- rownames(obj)                            
                      isNull        <- which(sapply(obj, FUN = function ( x ) { all(x==0)})==TRUE)
                      if (length (isNull)>0) {
                          obj <- obj[,-isNull]
                      }
                 }
                 cols  <- grep(paste0(".Dim",dims,"$" ), colnames(obj), value=TRUE)
                 tamMat<- obj[,c("item",cols)]
                 weg   <- which(tamMat[,2] == 0)
                 if(length(weg)>0) {tamMat <- tamMat[-weg,]}
                 shw6D <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = tamMat[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = colnames(qMatrix)[dims+1], par = "estSlope",  derived.par = NA, value = tamMat[,2], stringsAsFactors = FALSE)
                 if (ncol(tamMat) == 3 ) {
                     shw6se<- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = tamMat[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = colnames(qMatrix)[dims+1], par = "estSlope",  derived.par = "se", value = tamMat[,3], stringsAsFactors = FALSE)
                 }  else  {
                     shw6se<- NULL
                 }
                 return(rbind(shw6D, shw6se)) }))
         return(shw6)}

mergeDimensionIfDIF <- function(dat, qmatLong, datMergingVar, remove) {
         dat[,datMergingVar]    <- as.character(dat[,datMergingVar])
         dat[,"toMerge"] <- eatTools::halveString(dat[,datMergingVar], ":", first=TRUE)[,1]
         dat  <- merge(dat, qmatLong[,c("item", "dimensionName")], by.x = "toMerge", by.y = "item", all.x = TRUE)
         dat  <- dat[,-match(remove, colnames(dat))]
         return(dat)}

renameDifParameters <- function(dat, qmatLong) {
         indD5<- setdiff( 1:nrow(dat), grep(":DIF", dat[,"var1"]))
         indD5<- setdiff( dat[indD5,"var1"], qmatLong[,"item"])
         to   <- eatTools::removePattern(string = indD5, pattern = "DIF_")
         to1  <- eatTools::removeNumeric(to)
         to2  <- eatTools::removeNonNumeric(to)
         indD5<- data.frame ( from = indD5, to = paste(to1, to2, sep="_"), stringsAsFactors = FALSE)
         recSt<- paste("'",indD5[,"from"] , "' = '" , indD5[,"to"],"'", collapse="; ",sep="")
         indD <- grep(":DIF", dat[,"var1"])
         indD2<- eatTools::halveString(dat[indD,"var1"], pattern = ":DIF_", first=TRUE)
         indD3<- eatTools::removeNumeric(indD2[,2])
         indD4<- eatTools::removeNonNumeric(indD2[,2])
         dat[indD,"var1"] <- paste("item_", indD2[,1], "_X_", paste(indD3, indD4, sep="_"), sep="")
         dat[,"var1"]     <- car::recode(dat[,"var1"], recSt)
         return(dat)}

getTamInfit    <- function(runModelObj, qL, qMatrix, leseAlles, omitFit, seed) {
         if(leseAlles == FALSE || omitFit == TRUE ) {return(NULL)}
         infit<- tam.fit(runModelObj, progress=FALSE, seed=seed)                
         fits <- merge(infit[["itemfit"]], qL[,-match("value", colnames(qL))],  by.x = "parameter", by.y = colnames(qMatrix)[1], all = TRUE)
         if ( is.null(attr(runModelObj, "defineModelObj")[["all.Names"]][["DIF.var"]])) {           
              ret  <- rbind(data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = fits[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = fits[,"dimensionName"], par = "est",  derived.par = "infit", value = fits[,"Infit"], stringsAsFactors = FALSE),
                            data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = fits[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = fits[,"dimensionName"], par = "est",  derived.par = "outfit", value = fits[,"Outfit"], stringsAsFactors = FALSE))
         }  else  {                                                             
              ret  <- rbind(data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = infit$itemfit[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = NA, par = "est",  derived.par = "infit", value = infit$itemfit[,"Infit"], stringsAsFactors = FALSE),
                            data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = infit$itemfit[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = NA, par = "est",  derived.par = "outfit", value = infit$itemfit[,"Outfit"], stringsAsFactors = FALSE) )
              ret  <- mergeDimensionIfDIF(dat=ret, qmatLong=qL[,-match("value", colnames(qL))], datMergingVar="var1", remove = c("group", "toMerge"))
              colnames(ret) <- car::recode(colnames(ret), "'dimensionName'='group'")
              ret  <- renameDifParameters(dat=ret, qmatLong=qL[,-match("value", colnames(qL))])
         }
         return(ret)}

getTamPopPar    <- function(runModelObj, qMatrix, leseAlles) {
         if(leseAlles == FALSE ) {return(NULL)}
         if(ncol(qMatrix) == 2) {                                               
            ret  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = colnames(qMatrix)[2], var2 = NA , type = "distrpar", indicator.group = NA, group = "persons", par = "var",  derived.par = NA, value = runModelObj[["variance"]][1,1] , stringsAsFactors = FALSE)
         }  else  {                                                             
            cov1 <- runModelObj[["variance"]]
            colnames(cov1) <- colnames(qMatrix)[-1]
            rownames(cov1) <- colnames(qMatrix)[-1]
            cor1 <- cov2cor(cov1)
            for (ii in 1:nrow(cor1))   {                                        
                 cor1[ii,ii:ncol(cor1)] <- NA}
            cor1 <- reshape2::melt(cor1, measure.vars = colnames(cor1), na.rm=TRUE)
            vars <- Matrix::diag(cov1)
            ret  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = c(names(vars),as.character(cor1[,"Var1"])) , var2 = c(rep(NA, length(vars)), as.character(cor1[,"Var2"])) , type = "random", indicator.group = NA, group = "persons", par = c(rep("var",length(vars)), rep("correlation", nrow(cor1))) ,  derived.par = NA, value = c(unlist(vars), cor1[,"value"]), stringsAsFactors = FALSE)
         }
         return(ret)}

getTamRegPar    <- function(runModelObj, leseAlles, qMatrix, omitRegr, regr) {
         if(leseAlles == FALSE || omitRegr == TRUE ) {return(NULL)}
         if( !isTRUE(all.equal ( dim(runModelObj$beta) , c(1,1))))  {           
             regr <- data.frame ( reg.var = rownames(regr$beta), regr$beta, stringsAsFactors = FALSE)
             regr <- reshape2::melt(regr, id.vars = "reg.var", na.rm=TRUE)
             regr2<- data.frame ( par = "est", derived.par = car::recode(unlist(lapply(strsplit(as.character(regr[,"variable"]),"\\."), FUN = function ( l ) {l[1]})), "'se'='se'; else=NA"), group = colnames(qMatrix)[as.numeric(eatTools::removePattern( string = unlist(lapply(strsplit(as.character(regr[,"variable"]),"\\."), FUN = function ( l ) {l[2]})), pattern = "Dim")) + 1], regr, stringsAsFactors = FALSE)
             regr3<- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = regr2[,"reg.var"], var2 = NA , type = "regcoef", indicator.group = NA, group = regr2[,"group"], par = regr2[,"par"],  derived.par = regr2[,"derived.par"], value = regr2[,"value"] , stringsAsFactors = FALSE)
         }  else {
             return(NULL)
         }
         return(regr3)}

getTamModInd    <- function(runModelObj, leseAlles) {
         if(leseAlles == FALSE ) {return(NULL)}
         return(data.frame ( model = attr(runModelObj,"defineModelObj")[[ "analysis.name"]], source = "tam", var1 = NA, var2 = NA , type = "model", indicator.group = NA, group = NA, par = c("deviance", "Npar", "AIC", "BIC"), derived.par = NA, value = unlist(runModelObj[["ic"]][c("deviance", "Npars", "AIC", "BIC")]), stringsAsFactors = FALSE))}

getTamWles    <- function(runModelObj, qMatrix, leseAlles, omitWle) {
         if(leseAlles == FALSE || omitWle == TRUE ) {return(NULL)}
         txt  <- capture.output(wle  <- tam.wle(runModelObj, progress = FALSE)) 
         eind1<- ncol(wle) == 7                                                 
         if(isTRUE(eind1)) {
            cols1<- grep("theta$", colnames(wle))
            cols2<- grep("error$", colnames(wle))
            stopifnot(length(cols1) ==1, length(cols2) ==1)
            colnames(wle)[c(cols1,cols2)] <- paste(colnames(wle)[c(cols1,cols2)], ".Dim1", sep="")
         }
         weg1 <- grep("WLE.rel", colnames(wle))
         wleL <- reshape2::melt(wle, id.vars = "pid", measure.vars = colnames(wle)[-c(1:2,weg1)], na.rm=TRUE)
         wleL[,"group"] <- colnames(qMatrix)[as.numeric(eatTools::removePattern(string = unlist(lapply(strsplit(as.character(wleL[,"variable"]),"\\."), FUN = function (l) {l[2]})), pattern = "Dim"))+1]
         trans<- na.omit(unique(data.frame ( original = unlist(lapply(strsplit(as.character(wleL[,"variable"]),"\\."), FUN = function (l) {l[2]})), uebersetzt = wleL[,"group"], stringsAsFactors = FALSE)))
         wleL[,"par"]   <- car::recode(unlist(lapply(strsplit(as.character(wleL[,"variable"]),"\\."), FUN = function (l) {l[1]})), "'PersonScores'='NitemsSolved'; 'PersonMax'='NitemsTotal'; 'theta'='wle'; 'error'='wle'")
         wleL[,"derived.par"] <- car::recode(unlist(lapply(strsplit(as.character(wleL[,"variable"]),"\\."), FUN = function (l) {l[1]})), "'theta'='est'; 'error'='se';else=NA")
         rel  <- reshape2::melt(as.data.frame ( wle)[1,weg1], na.rm = TRUE)     
         if ( "variable" %in% colnames(rel)) {                                  
               rel[,"original"] <- unlist(lapply(strsplit(as.character(rel[,"variable"]),"\\."), FUN = function (l) {l[length(l)]}))
               rel  <- merge(rel, trans, by="original", all=TRUE)
         }  else  {
               rel[,"uebersetzt"] <- colnames(qMatrix)[-1]
         }
         res  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = c(wleL[,"pid"],rep(NA,nrow(rel))), var2 = NA , type = c(rep("indicator",nrow(wleL)), rep("tech", nrow(rel))), indicator.group = "persons", group = c(wleL[,"group"],rel[,"uebersetzt"]), par = c(wleL[,"par"],rep("wle", nrow(rel))),  derived.par = c(wleL[,"derived.par"],rep("rel",nrow(rel))), value = c(wleL[,"value"] ,rel[,"value"]), stringsAsFactors = FALSE)
         return(res)}

getTamPVs <- function ( runModelObj, qMatrix, leseAlles, omitPV, pvMethod, tam.pv.arguments) {
         if(omitPV == TRUE ) {return(NULL)}
         for ( i in names( tam.pv.arguments )) { assign(i, tam.pv.arguments[[i]]) }
         if(leseAlles == TRUE ) {
            if ( pvMethod == "regular" ) {                                      
                 do   <- paste ( "tam.pv ( ", paste(names(formals(tam.pv)), car::recode ( names(formals(tam.pv)), "'tamobj'='runModelObj'"), sep =" = ", collapse = ", "), ")",sep="")
            } else {                                                            
                 if ( is.null ( attr(runModelObj, "Y") ) ) {
                      warning("Conditioning model was not defined ('Y' is NULL).")
                      Y1 <- NULL
                 } else {
                      Y1 <- data.frame ( intercpt = 1, attr(runModelObj, "Y"))
                 }
                 do   <- paste ( "tam.pv.mcmc ( ", paste(names(formals(tam.pv.mcmc)), car::recode ( names(formals(tam.pv.mcmc)), "'tamobj'='runModelObj'; 'Y'='Y1'"), sep =" = ", collapse = ", "), ")",sep="")
            }
         }  else  {                                                             
            class(runModelObj) <- "list"                                        
            stopifnot ( pvMethod == "bayesian")
            do   <- paste ( "tam.pv.mcmc ( ", paste(names(formals(tam.pv.mcmc)), car::recode ( names(formals(tam.pv.mcmc)), "'tamobj'='runModelObj'; 'Y'='runModelObj[[\"Y\"]]'; 'nplausible'='attr(runModelObj, \"defineModelObj\")[[\"n.plausible\"]]'"), sep =" = ", collapse = ", "), ")",sep="")
         }
         pv   <- eval(parse(text=do))
         pvL  <- reshape2::melt(pv$pv, id.vars = "pid", na.rm=TRUE)
         pvL[,"PV.Nr"] <- as.numeric(eatTools::removePattern(string = unlist(lapply(strsplit(as.character(pvL[,"variable"]),"\\."), FUN = function (l) {l[1]})), pattern = "PV"))
         pvL[,"group"] <- colnames(qMatrix)[as.numeric(eatTools::removePattern(string = unlist(lapply(strsplit(as.character(pvL[,"variable"]),"\\."), FUN = function (l) {l[2]})), pattern = "Dim"))+1]
         res <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = pvL[,"pid"], var2 = NA , type = "indicator", indicator.group = "persons", group = pvL[,"group"], par = "pv",  derived.par = paste("pv", pvL[,"PV.Nr"],sep=""), value = pvL[,"value"] , stringsAsFactors = FALSE)
         return(res)}
         
getTamEAPs <- function ( runModelObj, qMatrix, leseAlles = leseAlles) {
         if(leseAlles == FALSE ) {return(NULL)}
         eaps <- runModelObj[["person"]]                                        
         eind1<- ncol(eaps) == 7                                                
         if(eind1 == TRUE) {                                                    
            cols <- grep("EAP$", colnames(eaps))                                
            stopifnot(length(cols) == 2)                                        
            colnames(eaps)[cols] <- paste(colnames(eaps)[cols], ".Dim1", sep="")
         }
         eaps <- reshape2::melt(eaps, id.vars = "pid", measure.vars = grep("EAP", colnames(eaps)), na.rm=TRUE)
         eaps[,"tam"]      <- eatTools::halveString(string = as.character(eaps[,"variable"]), pattern = "\\.", first = FALSE)[,"X2"]
         eaps[,"dimnumber"]<- as.numeric(eatTools::removePattern ( eaps[,"tam"], "Dim"))
         eaps[,"group"]    <- colnames(qMatrix)[eaps[,"dimnumber"] + 1]         
         checkmate::assert_numeric(unique(eaps[,"dimnumber"]), len = ncol(qMatrix)-1, any.missing = FALSE)
         eaps[,"par"]      <- "est"
         eaps[grep("^SD.",as.character(eaps[,"variable"])),"par"]   <- "se"
         res  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = eaps[,"pid"], var2 = NA , type = "indicator", indicator.group = "persons", group = eaps[,"group"], par = "eap", derived.par = eaps[,"par"], value = eaps[,"value"] , stringsAsFactors = FALSE)
         if (ncol(qMatrix)>2) {
             grp <-  eatTools::recodeLookup(names(runModelObj[["EAP.rel"]]), unique(eaps[,c("tam", "group")]))
         } else {
             stopifnot(length(unique(eaps[,"group"])) == 1)
             grp <- unique(eaps[,"group"])
         }
         res  <- rbind(res, data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = NA, var2 = NA , type = "tech", indicator.group = "persons", group = grp, par = "eap", derived.par = "rel", value = runModelObj[["EAP.rel"]] , stringsAsFactors = FALSE) )
         return(res)}

getTamQ3 <- function(runModelObj, leseAlles, shw1, Q3, q3MinObs, q3MinType){
         if(leseAlles == FALSE || Q3 == FALSE) {return(NULL)}
         nObs <- NULL
         if ( !is.null(q3MinObs) ) {                                            
              if ( q3MinObs > 1 ) {
                   nObs <- nObsItemPairs ( responseMatrix = runModelObj[["resp"]], q3MinType = q3MinType )
              }
              mat  <- tam.modelfit ( tamobj = runModelObj, progress = FALSE )
              matL <- reshapeQ3 (mat = mat$Q3.matr, q3MinObs = q3MinObs, nObs = nObs)
              if( nrow(matL)>0) {
                  res  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "tam", var1 = matL[,"Var1"], var2 = matL[,"Var2"] , type = "fixed",indicator.group = "items", group = paste(names(table(shw1[,"group"])), collapse="_"), par = "q3", derived.par = NA, value = matL[,"value"] , stringsAsFactors = FALSE)
              }  else  {
                  res  <- NULL
              }
         }
         return(res)}


getTamResults <- function(runModelObj, omitFit, omitRegr, omitWle, omitPV, nplausible , ntheta , normal.approx, samp.regr, theta.model, np.adj, Q3=Q3, q3MinObs =  q3MinObs, q3MinType = q3MinType,
                     pvMethod , group, beta_groups , level , n.iter , n.burnin, adj_MH , adj_change_MH , refresh_MH, accrate_bound_MH,	sample_integers, theta_init, print_iter , verbose, calc_ic, seed) {
         qMatrix<- attr(runModelObj, "defineModelObj")[["qMatrix"]]
         qL     <- reshape2::melt(qMatrix, id.vars = colnames(qMatrix)[1], variable.name = "dimensionName", na.rm=TRUE)
         qL     <- qL[which(qL[,"value"] != 0 ) , ]
         varName<- colnames(qMatrix)[1]                                         
         if( omitRegr == FALSE && !inherits(runModelObj, "tamBayes")) {
             beg <- Sys.time()
             txt <- capture.output ( regr <- tam.se(runModelObj))               
             stopifnot ( nrow(regr$beta) == ncol(attr(runModelObj, "Y") )+1)    
             rownames(regr$beta) <- c("(Intercept)", colnames(attr(runModelObj, "Y")))
         } else {
             regr <- NULL
         }
         if ( !inherits(runModelObj, "tamBayes") ) {leseAlles <- TRUE} else {leseAlles <- FALSE}
         ret    <- NULL                                                         
         resItem<- getTamItempars(runModelObj=runModelObj, qL=qL, qMatrix=qMatrix, leseAlles = leseAlles)
         ret    <- rbind(ret, resItem[["shw1"]], resItem[["shw2"]])
         ret    <- rbind(ret, getTamDescriptives(runModelObj=runModelObj, qL=qL, qMatrix=qMatrix, leseAlles = leseAlles))
         ret    <- rbind(ret, getTamDiscrim(runModelObj=runModelObj, qL=qL, qMatrix = qMatrix, leseAlles = leseAlles))
         ret    <- rbind(ret, getTam2plDiscrim(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, regr = regr, omitRegr=omitRegr))
         beg    <- Sys.time()
         ret    <- rbind(ret, getTamInfit(runModelObj=runModelObj, qL=qL, qMatrix = qMatrix, leseAlles = leseAlles, omitFit = omitFit, seed=seed))
         ret    <- rbind(ret, getTamPopPar(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles))
         ret    <- rbind(ret, getTamRegPar(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, omitRegr = omitRegr, regr=regr))
         ret    <- rbind(ret, getTamModInd(runModelObj=runModelObj, leseAlles = leseAlles))
         beg    <- Sys.time()
         ret    <- rbind(ret, getTamWles(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, omitWle = omitWle))
         beg    <- Sys.time()
         tamArg <- as.list(match.call(definition = getTamResults))
         weg    <- which(names(tamArg) %in% c("tamobj", "Y", "runModelObj", "qMatrix", "leseAlles", "omitPV", "pvMethod", "omitFit", "omitRegr", "omitWle"))
         if ( length(weg)>0) {tamArg <- tamArg[-weg]}
         tamarg <- list()                                                       
         for ( i in 2:length(tamArg)) {
              tamarg[[names(tamArg)[i]]] <- eval(tamArg[[i]])
         }
         retPVs <- getTamPVs ( runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, omitPV = omitPV, pvMethod = pvMethod, tam.pv.arguments = tamarg)
         ret    <- rbind(ret, retPVs)
         ret    <- rbind(ret, getTamEAPs(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles))
         ret    <- rbind(ret, getTamQ3(runModelObj=runModelObj, leseAlles = leseAlles, shw1 = resItem[["shw1"]], Q3=Q3, q3MinObs=q3MinObs, q3MinType=q3MinType))
         return(ret)}


reshapeQ3 <- function ( mat, q3MinObs, nObs ) {
             for (ii in 1:(nrow(mat)-1)) { mat[ii,ii:ncol(mat)] <- NA}          
             matL <- reshape2::melt ( mat , na.rm = TRUE)                       
             if ( !is.null(nObs)) {                                             
                   check<- do.call("rbind", apply(matL[,-ncol(matL)], MARGIN = 1, FUN = function ( y ) { ret <- sort ( y); ret <- data.frame ( Var1 = ret[1], Var2 = ret[2], stringsAsFactors = FALSE); return(ret)}))
                   matL <- data.frame ( check, value = matL[,"value"], stringsAsFactors = FALSE)
                   matL <- merge ( matL, nObs, by = c("Var1", "Var2"), all = TRUE)
                   matL <- matL[which(!is.na(matL[,"value"])),]
                   weg  <- which(matL[,"minValue"] < q3MinObs)
                   if (length(weg)>0) { matL <- matL[-weg,]}
             }
             if ( nrow(matL) == 0 ) {
                   cat("No observations left in Q3 matrix.\n")
             }
             return(matL)}

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

normalize.path <- function(string)
                  {string <- gsub("//","/",string)
                   string <- gsub("/","//",string)
                   string <- gsub("//","\\\\",string)
                   return(string)}

gen.syntax     <- function(Name,daten, all.Names, namen.all.hg = NULL, all.hg.char = NULL, var.char, model = NULL, anchored, constraints=c("cases","none","items"), pfad=NULL, Title=NULL,n.plausible=5,std.err=c("quick","full","none"), compute.fit ,
                           distribution=c("normal","discrete"), method=c("gauss", "quadrature", "montecarlo"), n.iterations=200, nodes=NULL, p.nodes=2000, f.nodes=2000, converge=0.001,deviancechange=0.0001, equivalence.table=c("wle","mle","NULL"), use.letters=use.letters, model.statement=model.statement, conquest.folder = NULL, allowAllScoresEverywhere,
                           seed , export = list(logfile = TRUE, systemfile = FALSE, history = TRUE, covariance = TRUE, reg_coefficients = TRUE, designmatrix = FALSE) )  {
                   if(is.null(anchored)) {anchored <- FALSE} else {anchored <- TRUE}
                   export.default <- list(logfile = TRUE, systemfile = FALSE, history = TRUE, covariance = TRUE, reg_coefficients = TRUE, designmatrix = FALSE)
                   mustersyntax <- c("title = ####hier.title.einfuegen####;",
                   "export logfile >> ####hier.name.einfuegen####.log;",
                   "datafile ####hier.name.einfuegen####.dat;",
                   "Format pid ####hier.id.einfuegen####",
                   "group",
                   "codes ####hier.erlaubte.codes.einfuegen####;",
                   "labels  << ####hier.name.einfuegen####.lab;",
                   "import anchor_parameters << ####hier.name.einfuegen####.ank;",
                   "caseweight",
                   "set constraints=####hier.constraints.einfuegen####;",
                   "set warnings=no,update=yes,n_plausible=####hier.anzahl.pv.einfuegen####,p_nodes=####hier.anzahl.p.nodes.einfuegen####,f_nodes=####hier.anzahl.f.nodes.einfuegen####;",
                   "set seed=####hier.seed.einfuegen####;",
                   "export par    >> ####hier.name.einfuegen####.prm;",
                   "regression",
                   "model ####hier.model.statement.einfuegen####;",
                   "estimate ! fit=####hier.fitberechnen.einfuegen####,method=####hier.method.einfuegen####,iter=####hier.anzahl.iterations.einfuegen####,nodes=####hier.anzahl.nodes.einfuegen####,converge=####hier.converge.einfuegen####,deviancechange=####hier.deviancechange.einfuegen####,stderr=####hier.std.err.einfuegen####,distribution=####hier.distribution.einfuegen####;",
                   "Itanal >> ####hier.name.einfuegen####.itn;",
                   "show cases! estimates=latent >> ####hier.name.einfuegen####.pvl;",
                   "show cases! estimate=wle >> ####hier.name.einfuegen####.wle;",
                   "equivalence ####hier.equivalence.table.einfuegen#### >> ####hier.name.einfuegen####.equ;",
                   "show >> ####hier.name.einfuegen####.shw;",
                   "export history >> ####hier.name.einfuegen####.his;",
									 "export covariance >> ####hier.name.einfuegen####.cov;",
									 "export reg_coefficients >> ####hier.name.einfuegen####.reg;",
									 "export designmatrix >> ####hier.name.einfuegen####.mat;",
                   "put >> ####hier.name.einfuegen####.cqs;  /* export systemfile */",
                   "descriptives !estimates=pv >> ####hier.name.einfuegen####_pvl.dsc;",
                   "descriptives !estimates=wle >> ####hier.name.einfuegen####_wle.dsc;",
                   "quit;")
                   if(is.null(Title))   {                                       
                      all.inf  <- Sys.getenv()
                      Title    <- paste("Analysis name: ",Name, ", User: ",all.inf["USERNAME"],", Computername: ",all.inf["COMPUTERNAME"],", ", R.version$version.string , ", Time: ",date(),sep="")}
                   converge <- paste("0",substring(as.character(converge+1),2),sep="")
                   deviancechange <- paste("0",substring(as.character(deviancechange+1),2),sep="")
                   syntax    <- gsub("####hier.title.einfuegen####",Title,mustersyntax)
                   if(is.null(n.plausible))   {n.plausible <- 0}  ; if(is.na(n.plausible))     {n.plausible <- 0}
                   if(n.plausible == 0 )     {                                  
                      syntax    <- gsub("n_plausible=####hier.anzahl.pv.einfuegen####,","",syntax) } else {
                      syntax    <- gsub("####hier.anzahl.pv.einfuegen####",n.plausible,syntax)
                   }
                   syntax    <- gsub("####hier.name.einfuegen####",Name,syntax)
                   ID.char   <- max(as.numeric(names(table(nchar(daten[,"ID"])))))
                   syntax    <- gsub("####hier.id.einfuegen####",paste("1-",as.character(ID.char)," ",sep="" ) ,syntax)
                   syntax    <- gsub("####hier.anzahl.iterations.einfuegen####",n.iterations,syntax)
                   syntax    <- gsub("####hier.anzahl.p.nodes.einfuegen####",p.nodes,syntax)
                   syntax    <- gsub("####hier.anzahl.f.nodes.einfuegen####",f.nodes,syntax)
                   syntax    <- gsub("####hier.converge.einfuegen####",converge,syntax)
                   syntax    <- gsub("####hier.deviancechange.einfuegen####",deviancechange,syntax)
                   if(!is.null(seed)) {syntax    <- gsub("####hier.seed.einfuegen####",seed,syntax)}
                   syntax    <- gsub("####hier.constraints.einfuegen####",match.arg(constraints),syntax)
                   compute.fit  <- if(compute.fit == TRUE ) compute.fit <- "yes" else compute.fit <- "no"
                   syntax    <- gsub("####hier.fitberechnen.einfuegen####",compute.fit,syntax)
                   syntax    <- gsub("####hier.anzahl.nodes.einfuegen####",nodes,syntax)
                   syntax    <- gsub("####hier.std.err.einfuegen####",match.arg(std.err),syntax)
                   syntax    <- gsub("####hier.distribution.einfuegen####",match.arg(distribution),syntax)
                   syntax    <- gsub("####hier.equivalence.table.einfuegen####",match.arg(equivalence.table),syntax)
                   syntax    <- gsub("####hier.model.statement.einfuegen####",tolower(model.statement),syntax)
                   erlaubte.codes <- paste(gsub("_","",sort(gsub(" ","_",formatC(names(eatTools::tableUnlist(daten[, all.Names[["variablen"]], drop=FALSE ])),width=var.char)),decreasing=TRUE)),collapse=",")
                   syntax    <- gsub("####hier.erlaubte.codes.einfuegen####",erlaubte.codes, syntax )
                   ind       <- grep("Format pid",syntax)
                   beginn    <- NULL                                            
                   if(length(namen.all.hg)>0)    {                              
                     all.hg.char.kontroll <- all.hg.char
                     all.hg.char <- sapply(namen.all.hg, FUN=function(ii) {max(nchar(as.character(na.omit(daten[,ii]))))})
                     stopifnot(all(all.hg.char == all.hg.char.kontroll))        
                     for (ii in 1:length(namen.all.hg))  {
                          if(is.null(beginn)) {beginn <- ID.char+1}
                          ende   <- beginn-1+all.hg.char[ii]
                          if (beginn != ende) {syntax[ind] <- paste(syntax[ind],namen.all.hg[ii], " ", beginn,"-",ende," ",sep="")}
                          if (beginn == ende) {syntax[ind] <- paste(syntax[ind],namen.all.hg[ii], " ", beginn," ",sep="")}
                          beginn  <- ende+1 }
                   }
                   if(length(all.Names[["DIF.var"]])>0)   {                     ### in folgender Schleife ueberschrieben und dann in der Schleife "if(!is.null(HG.var))" ergaenzt, nicht neu geschrieben
                      if(model.statement != "item") {
                        cat(paste("Caution! DIF variable was specified. Expected model statement is: 'item - ",tolower(all.Names[["DIF.var"]])," + item*",tolower(all.Names[["DIF.var"]]),"'.\n",sep=""))
                        cat(paste("However, '",tolower(model.statement),"' will used as 'model statement' to accomplish your will.\n",sep=""))
                      }
                      if(model.statement == "item") {
                         ind.model <- grep("model item", syntax)                
                         stopifnot(length(ind.model)==1)
                         syntax[ind.model] <- paste("model item - ",paste(tolower(all.Names[["DIF.var"]]),collapse=" - ") ," + ", paste("item*",tolower(all.Names[["DIF.var"]]),collapse=" + "), ";",sep="")
                      }
                   }
                   if(length(all.Names[["HG.var"]])>0)  {
                      ind.2   <- grep("^regression$",syntax)
                      syntax[ind.2] <- paste(eatTools::crop(paste( c(syntax[ind.2], tolower(all.Names[["HG.var"]])), collapse=" ")),";",sep="")
                      if(method == "gauss") {warning("Gaussian quadrature is only available for models without latent regressors.\n         Use 'Bock-Aiken quadrature' for estimation.")
                                             method <- "quadrature"} }          ### method muss "quadrature" oder "montecarlo" sein
                   syntax    <- gsub("####hier.method.einfuegen####",method,syntax)
                   if(length(all.Names[["weight.var"]])>0)  {                   
                      ind.4   <- grep("caseweight",syntax)
                      syntax[ind.4] <- paste( syntax[ind.4], " ", tolower(all.Names[["weight.var"]]),";",sep="") }
                   if(length(all.Names[["group.var"]])>0) {
                       ind.3   <- grep("^group$",syntax)
                       stopifnot(length(ind.3) == 1)
                       syntax[ind.3] <- paste(eatTools::crop(paste( c(syntax[ind.3], tolower(all.Names[["group.var"]])), collapse=" ")),";",sep="")
                       add.syntax.pv  <- as.vector(sapply(all.Names[["group.var"]], FUN=function(ii) {paste("descriptives !estimates=pv, group=",tolower(ii)," >> ", Name,"_",tolower(ii),"_pvl.dsc;",sep="")} ))
                       add.syntax.wle <- as.vector(sapply(all.Names[["group.var"]], FUN=function(ii) {paste("descriptives !estimates=wle, group=",tolower(ii)," >> ", Name,"_",tolower(ii),"_wle.dsc;",sep="")} ))
                       ind.3    <- grep("quit",syntax)
                       stopifnot(length(ind.3)==1)
                       syntax   <- c(syntax[1:(ind.3-1)],add.syntax.pv, add.syntax.wle, syntax[ind.3:length(syntax)]) }
                   if(is.null(beginn)) {beginn <- ID.char+1}
                   syntax[ind] <- paste(syntax[ind], "responses ",beginn,"-",beginn-1+var.char*ncol(data.frame(daten[,all.Names[["variablen"]]],stringsAsFactors = FALSE)),";",sep="")
                   if(var.char>1)  {                                            
                      syntax[ind] <- paste(gsub(";","",syntax[ind]), " (a",var.char,");",sep="")}
                   score.statement <- .writeScoreStatementMultidim (data=daten, itemCols=all.Names[["variablen"]], qmatrix=model, columnItemNames = 1 ,use.letters=use.letters, allowAllScoresEverywhere = allowAllScoresEverywhere )
                   expected.nodes  <- nodes^(ncol(model)-1)
                   if(expected.nodes>3500 & method != "montecarlo") {cat(paste("Specified model probably will use ",expected.nodes," nodes. Choosen method ",method," may not appropriate. Recommend to use 'montecarlo' instead.\n",sep=""))}
                   ind <- grep("labels ",syntax)
                   stopifnot(length(ind)==1)
                   syntax <- c(syntax[1:ind],score.statement,syntax[(ind+1):length(syntax)])
                   if(length(all.Names[["HG.var"]])==0) {                       
                      ind.2 <- grep("^regression$",syntax)
                      stopifnot(length(ind.2)==1)
                      syntax <- syntax[-ind.2]
                      ind.3 <- grep("export reg_coefficients",syntax)
                      stopifnot(length(ind.3)==1)
                      syntax <- syntax[-ind.3] }
                   if(length(all.Names[["group.var"]]) ==0) {                   
                      ind.3 <- grep("^group$",syntax)
                      stopifnot(length(ind.3)==1)
                      syntax <- syntax[-ind.3]}
                   if(length(all.Names[["weight.var"]]) ==0) {                  
                      ind.4 <- grep("^caseweight$",syntax)
                      stopifnot(length(ind.4)==1)
                      syntax <- syntax[-ind.4]}
                   if(match.arg(equivalence.table) == "NULL") {                 
                      ind.5   <- grep("^equivalence",syntax)
                      stopifnot(length(ind.5)==1)
                      syntax <- syntax[-ind.5]}
                   if(is.null(seed)) {                                          
                      ind.7   <- grep("^set seed",syntax)
                      stopifnot(length(ind.7)==1)
                      syntax <- syntax[-ind.7]}
                   if(n.plausible == 0)     {                                   
                      ind.6   <- grep("^show cases! estimates=latent", syntax)
                      stopifnot(length(ind.6) == 1)
                      syntax  <- syntax[-ind.6]}
                   if(anchored == FALSE) {ind.2 <- grep("anchor_parameter",syntax)
                                        syntax <- syntax[-ind.2]}
                   if(anchored == TRUE)  {ind.2 <- grep("^set constraints",syntax)# wenn ANKER gesetzt, setze constraints auf "none"
                                        if(match.arg(constraints) != "none") { cat("Anchorparameter were defined. Set constraints to 'none'.\n")}
                                        syntax[ind.2]  <- "set constraints=none;"}
                   if(!all(sapply(export, inherits, what="logical"))) {stop("All list elements of argument 'export' have to be of class 'logical'.")}
                   export <- as.list(userSpecifiedList ( l = export, l.default = export.default ))
                   weg <- names(export[which(export == FALSE)])
                   if(length(weg)>0)    {                                       
                      for (ii in seq(along=weg) ) {
                           ind.x <- grep(paste("export ", weg[ii], sep=""), syntax)
                           stopifnot(length(ind.x) == 1)
                           syntax <- syntax[-ind.x]}}
                   write(syntax,file.path(pfad,paste(Name,".cqc",sep="")),sep="\n")}

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

isLetter <- function ( string ) {
            splt <- strsplit(string, "")
            isL  <- lapply(splt, FUN = function ( x ) {
                    ind <- which ( x %in% c( letters , LETTERS ))
                    x[setdiff(1:length(x),ind)] <- " "
                    x <- eatTools::crop(paste(x, sep="", collapse=""))
                    x <- unlist ( strsplit(x, " +") )
                    return(x)  } )
            return(isL)}

.writeScoreStatementMultidim <- function(data, itemCols, qmatrix, columnItemNames = 1 ,columnsDimensions = -1, use.letters=use.letters , allowAllScoresEverywhere) {
            n.dim      <- (1:ncol(qmatrix) )[-columnItemNames]                  
            stopifnot(length( which( rowSums(qmatrix[,n.dim,drop = FALSE]) == 0))==0)
      	    if(length(setdiff(names(eatTools::tableUnlist(qmatrix[,-1, drop = FALSE])), c("0","1"))) > 0 )  {
               cat("Found unequal factor loadings for at least one dimension. This will result in a 2PL model.\n")
               for (u in 2:ncol(qmatrix)) {qmatrix[,u] <- as.character(round(qmatrix[,u], digits = 3))}
            }                                                                   
            stopifnot(all(qmatrix[,1] == itemCols))                             
            cat(paste("Q matrix specifies ",length(n.dim)," dimension(s).\n",sep=""))
            stopifnot(length(setdiff(colnames(data[,itemCols]),  qmatrix[,columnItemNames]) )==0)
            unique.patter <- qmatrix[which(!duplicated(do.call("paste", qmatrix[,-1, drop = FALSE] ))), -1, drop = FALSE]
            colnames(unique.patter) <- paste("Var",1:ncol(unique.patter), sep="")
            score.matrix  <- data.frame(score=1, unique.patter, matrix(NA, nrow= nrow(unique.patter), ncol=length(itemCols), dimnames=list(NULL, paste("X",1:length(itemCols),sep=""))),stringsAsFactors = FALSE)
            scoreColumns  <- grep("^Var",colnames(score.matrix))
            for (i in 1:length(itemCols))  {                                    
               qmatrix.i    <- qmatrix[qmatrix[,columnItemNames] == itemCols[i],]
               matchRow     <- which(sapply ( 1:nrow(score.matrix) , function(ii) {all ( as.numeric(qmatrix.i[,n.dim]) == as.numeric(score.matrix[ii,scoreColumns])) }))
               stopifnot(length(matchRow) == 1)
               matchColumn  <- min(which(is.na(score.matrix[matchRow,])))       
               stopifnot(length(matchColumn) == 1)
               score.matrix[matchRow,matchColumn] <- i
		        }
            rowsToDelete <- which(is.na(score.matrix[, max(scoreColumns) + 1])) 
            if(length(rowsToDelete)>0) {score.matrix <- score.matrix[-rowsToDelete, ]}
            for (ii in 1:nrow(score.matrix)) {score.matrix[,ii] <- as.character(score.matrix[,ii])}
            score.matrix <- fromMinToMax(dat = data[,itemCols, drop = FALSE], score.matrix = score.matrix, qmatrix = qmatrix, allowAllScoresEverywhere = allowAllScoresEverywhere, use.letters = use.letters)
            kollapse <- lapply(1:nrow(score.matrix), FUN=function(ii) {na.omit(as.numeric(score.matrix[ii,-c(1,scoreColumns)]))})
            kollapse.diff   <- lapply(kollapse,FUN=function(ii) {c(diff(ii),1000)})
            kollapse.ascend <- lapply(kollapse.diff, FUN=function(ii) {unique(c(0, which(ii!=1)))})
            kollapse.string <- list()
            for (a in 1:length(kollapse.ascend))  {
                string   <- list()
                for (i in 2:length(kollapse.ascend[[a]]))   {
                    string.i <- unique( c(kollapse[[a]][kollapse.ascend[[a]][i-1]+1], kollapse[[a]][kollapse.ascend[[a]][i]]))
                    string.i <- ifelse(length(string.i) == 2,paste(string.i[1],"-",string.i[2],sep=""),as.character(string.i))
                    string[[i]] <- string.i
				        }
                string <- paste(unlist(string),collapse=", ")
                kollapse.string[[a]] <- string
			      }
            control <- lapply(kollapse.string,FUN=function(ii) {eval(parse(text=paste("c(",gsub("-",":",ii),")",sep="")))})
            if (!all(unlist(lapply(1:length(control), FUN=function(ii) {all(kollapse[[ii]] == control[[ii]])})))) {
                cat("Error in creating score statement.\n")
			      }
            score.matrix <- data.frame(prefix="score",score.matrix[,c(1,scoreColumns)],items="! items(",kollapse.string=unlist(kollapse.string),suffix=");",stringsAsFactors=F)
            score.statement <- sapply(1:nrow(score.matrix), FUN=function(ii) { paste(score.matrix[ii,],collapse=" ")})
            return(score.statement) }


fromMinToMax <- function(dat, score.matrix, qmatrix, allowAllScoresEverywhere, use.letters)    {
                all.values <- plyr::alply(as.matrix(score.matrix), .margins = 1, .fun = function(ii) {sort(names(eatTools::tableUnlist(dat[,na.omit(as.numeric(ii[grep("^X", names(ii))])), drop = FALSE])) ) })
                if ( length(all.values) > 1) {                                  
                     if ( all ( outer ( all.values, all.values, Vectorize(identical))) == FALSE ) {
                          cat(paste("Found different values for dimensions: \n",sep=""))
                          for ( u in 1:length(all.values)) {
                               cat(paste0("   Dimension ", u, ": values '",paste(all.values[[u]], collapse= "', '"), "' \n"))
                          }
                          if ( allowAllScoresEverywhere == TRUE ) {
                               all.values <- lapply(all.values, FUN = function ( ii ) { sort(unique( unlist ( all.values ) ))})
                               cat(paste("Following value definition was done according to 'allowAllScoresEverywhere == TRUE': \n",sep=""))
                               for ( u in 1:length(all.values)) {
                                    cat(paste0("   Dimension ", u, ": values '",paste(all.values[[u]], collapse= "', '"), "' \n"))
                               }
                          }
                     }
                }
                if(use.letters == TRUE )  {minMaxRawdata  <- unlist ( lapply( all.values, FUN = function (ii) {paste("(",paste(LETTERS[which(LETTERS == ii[1]) : which(LETTERS == ii[length(ii)])], collapse=" "),")") } ) ) }
                if(use.letters == FALSE ) {minMaxRawdata  <- unlist ( lapply( all.values, FUN = function (ii) {paste("(",paste(ii[1] : ii[length(ii)],collapse = " "),")")  } ) ) }
                scoring <- unlist( lapply( minMaxRawdata , FUN = function(ii) { paste("(", paste( 0 : (length(unlist(strsplit(ii, " ")))-3), collapse = " "),")")}) )
                stopifnot(length(scoring) == length( minMaxRawdata ), length(scoring) == nrow(score.matrix )  )
                for (i in 1:nrow(score.matrix))    {
                    score.matrix$score[i] <- minMaxRawdata[i]
                    targetColumns         <- suppressWarnings(intersect ( grep("Var",colnames(score.matrix)), which(as.numeric(score.matrix[i,]) != 0 ) ))
                    stopifnot(length(targetColumns) > 0 )
                    score.matrix[i,targetColumns]  <- suppressWarnings(unlist(lapply(score.matrix[i,targetColumns], FUN = function ( y ) {paste( "(", paste(as.numeric(y) * na.omit(as.numeric(unlist(strsplit(scoring[i]," ")))), collapse = " "), ")")})))
                    nonTargetColumns      <- suppressWarnings(intersect ( grep("Var",colnames(score.matrix)), which(as.numeric(score.matrix[i,]) == 0 ) ))
                    if ( length ( nonTargetColumns ) > 0 )    {
                       score.matrix[i,nonTargetColumns]  <- "()"
                    }
                }
                return(score.matrix)}

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

userSpecifiedList <- function ( l, l.default ) {
		if ( !is.null ( names ( l ) ) ) {
				names ( l ) <- match.arg ( names(l) , names(l.default) , several.ok = TRUE )
		} else {
        if(length(l) > length(l.default) )  {
           stop("Length of user-specified list with more elements than default list.\n")
        }
				names ( l ) <- names ( l.default )[seq(along=l)]
		}
		if ( length(l) < length(l.default) ) {
				l <- c ( l , l.default )
				l <- l[!duplicated(names(l))]
				l <- l[match ( names (l) , names(l.default) )]
		}
		return(l)}

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

item.diskrim <- function(daten, itemspalten, streng = TRUE) {
                 if(!missing(itemspalten))  {daten <- daten[,itemspalten]}      
                 trenn <- suppressWarnings(eatTools::pwc(daten))                
                 if(streng) {return(data.frame(item.name=trenn[,"item"],item.diskrim = trenn[,"partWholeCorr"],stringsAsFactors = FALSE))} else {return(data.frame(item.name=trenn[,"item"],item.diskrim = trenn[,"corr"],stringsAsFactors = FALSE))}}
                 
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

plotDevianceConquest <- function ( logFile, omitUntil = 1, reverse = TRUE, change = TRUE ) {
           if ( inherits(logFile, "character")) {lf <- logFile}  else  { lf <- file.path(logFile[["path"]], paste0(logFile[["analysis.name"]], ".log"))}
           input<- scan(lf,what="character",sep="\n",quiet=TRUE)
           ind  <- grep("eviance=", input)
           dev  <- unlist(lapply(input[ind], FUN = function (x) {               
                   brace <- grep("\\(", x)                                      
                   if(length(brace)>0) {
                       weg <- grep("\\(", unlist(strsplit(x, "")))              
                       x   <- substr(x, 1, weg-1)
                   }
                   return(x)}))
           dev  <- data.frame ( lapply(data.frame ( eatTools::halveString(dev, "\\."), stringsAsFactors = FALSE), eatTools::removeNonNumeric), stringsAsFactors = FALSE)
           mat  <- data.frame ( iter = 1:length(ind), dev = as.numeric(paste(dev[,1], dev[,2], sep=".")), stringsAsFactors = FALSE)
           if(omitUntil>0)  {
              dc<- mat[-c(1:omitUntil),2]                                       
           } else {
              dc<- mat[,2]
           }
           if ( change ){
              dc<- diff(dc)
              yl<- "Deviance Change"                                            
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
           plot ( dc[,1], dc[,2], type="o",
                main=titel,  xlab="Iteration",
                xlim=c(min(dc[,1]),max(dc[,1])),  xaxp=c(0,xm,xt),
                ylab=yl, pch=20, cex=cex, lwd=0.75 )
           si   <- devtools::session_info(pkgs = "eatModel")
           si   <- si[["packages"]][which(si[["packages"]][,"package"] == "eatModel"),]
           sysi <- Sys.info()
           stri <- paste0("'eatModel', version ", si[["loadedversion"]], ", build ",si[["date"]], ", user: ",sysi[["user"]], " (", sysi[["sysname"]],", ",sysi[["release"]], ", ",sysi[["version"]], ")")
           if (class(logFile) == "list") {
               stri <- paste0("Method = '",logFile[["ret"]][[1]],"'  |  nodes = ",logFile[["ret"]][[2]],"  |  ",capture.output(logFile[["tme"]]), "\n",stri)
           }
           graphics::mtext(stri)
           graphics::abline( a=0, b=0 )
           dcr  <- dc[dc[,2]<0,]
           graphics::points( dcr[,1], dcr[,2], pch=20, cex=cex, col="red") }


