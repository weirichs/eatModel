## usually called after defineModel()

getResults <- function(runModelObj, overwrite = FALSE, Q3 = TRUE, q3theta = c("pv", "wle", "eap"),
                       q3MinObs = 0, q3MinType = c("singleObs", "marginalSum"),
                       omitFit = FALSE, omitRegr = FALSE, omitWle = FALSE, omitPV = FALSE,
                       abs.dif.bound = 0.6, sig.dif.bound = 0.3, p.value = 0.9,
                       nplausible = NULL, ntheta = 2000, normal.approx = FALSE,
                       samp.regr = FALSE, theta.model=FALSE, np.adj=8, group = NULL,
                       beta_groups = TRUE, level = .95, n.iter = 1000, n.burnin = 500,
                       adj_MH = .5, adj_change_MH = .05, refresh_MH = 50,
                       accrate_bound_MH = c(.45, .55),	sample_integers=FALSE,
                       theta_init=NULL, print_iter = 20, verbose = TRUE, calc_ic=TRUE,
                       omitUntil=1, seed=NA) {

### checks ---------------------------------------------------------------------
  checkmate::assert_numeric(q3MinObs, lower = 0, len=1)
  
  q3MinType<- match.arg(q3MinType)
  q3theta  <- match.arg(q3theta)

  # logical arguments
  lapply(c(overwrite, Q3, omitFit, omitRegr, omitWle, omitPV),
         checkmate::assert_logical, len = 1)
  #lapply(c(), checkmate::assert_logical, len = 1)

  # if DIF was applied
  lapply(c(abs.dif.bound, sig.dif.bound, p.value), checkmate::assert_numeric, len = 1)

  # if software = "tam"
  lapply(c(normal.approx, samp.regr, theta.model, beta_groups, sample_integers,
           verbose, calc_ic), checkmate::assert_logical, len = 1)
  checkmate::assert_numeric(nplausible, len = 1, null.ok = TRUE)
  checkmate::assert_vector(group, null.ok = TRUE)
  checkmate::assert_matrix(theta_init, null.ok = TRUE)

  lapply(c(ntheta, np.adj, level, n.iter, n.burnin, adj_MH, adj_change_MH,
           refresh_MH, print_iter), checkmate::assert_numeric, len = 1)
  checkmate::assert_numeric(accrate_bound_MH)

  # omitUntil (given to plotDevianceConquest)
  checkmate::assert_numeric(omitUntil, len = 1)

### function -------------------------------------------------------------------

     ### Mehrmodellfall: rekursiver Aufruf
            if(inherits(runModelObj, "runMultiple")) {
                if(is.null ( attr(runModelObj, "split")[["nCores"]] ) || attr(runModelObj, "split")[["nCores"]] == 1 ) {
                   res <- lapply( runModelObj, FUN = function ( r ) {           ### erstmal single core auswertung
                          do  <- paste ( "getResults ( ", paste(names(formals(getResults)), car::recode(names(formals(getResults)), "'runModelObj'='r'"), sep =" = ", collapse = ", "), ")",sep="")
                          ret <- eval(parse(text=do))
                          return(ret)})
                   }  else  {
                          beg <- Sys.time()
                          if ( attr(runModelObj, "split")[["mcPackage"]] == "parallel") {
                               cl  <- makeCluster(attr(runModelObj, "split")[["nCores"]], type = "SOCK")
                          }  else  {
                               cl  <- future::makeClusterPSOCK(attr(runModelObj, "split")[["nCores"]], verbose=FALSE)
                          }
                          doMulticoreResults <- function (laufnummer, ... ) {
                              if(!"eatModel" %in% .packages()) {library(eatModel)}
                              if(!"TAM" %in% .packages()) {library(TAM)}
                              if(!"mirt" %in% .packages()) {library(mirt)}
                              do  <- paste ( "getResults ( ", paste(names(formals(getResults)), car::recode(names(formals(getResults)), "'runModelObj'='runModelObj[[laufnummer]]'"), sep =" = ", collapse = ", "), ")",sep="")
                              ret <- eval(parse(text=do))
                              return(ret)}
                          do  <- paste("clusterApply(cl = cl, x = 1:length(runModelObj), fun = doMulticoreResults,",paste(names(formals(getResults))[-1], names(formals(getResults))[-1], sep=" = ", collapse=", "),")", sep="")
                          res <- eval(parse(text=do))
                          stopCluster(cl)
                          cat(paste ( "Results of ",length(runModelObj), " analyses processed: ", sep="")); print( Sys.time() - beg, digits=3)
                   }
               res <- do.call("rbind", res )
               class(res) <- c("multipleResults", "data.frame")
               rownames(res) <- NULL
               return(res)
     ### hier ist der rekursive Aufruf beendet: das folgende geschieht fuer jedes Modell einzeln, technisch auf zweierlei Weisen, je nachdem ob Conquest oder TAM gerechnet wurde
            }  else {                                                           ### Einmodellfall
               if ( is.null(runModelObj)) {return(NULL)}
               isTa  <- FALSE
     ### Software = Conquest
               if(isTRUE(attr(runModelObj, "software") == "conquest") || is_true(obj = runModelObj[["software"]], crit = "conquest")) {
                    if ( isTRUE(Q3) ) {
                        if ( ncol ( runModelObj[["qMatrix"]]) !=2 ) {
                            cat("Q3 is only available for unidimensional models. Estimation will be skipped.\n")
                            Q3 <- FALSE
                        }
                    }
                    do  <- paste ( "res <- getConquestResults ( ", paste(names(formals(getConquestResults)), car::recode(names(formals(getConquestResults)), "'path'='runModelObj$dir'; 'analysis.name'='runModelObj$analysis.name'; 'model.name'='runModelObj$model.name'; 'qMatrix'='runModelObj$qMatrix'; 'all.Names'='runModelObj$all.Names'; 'renam'='runModelObj$renam'; 'deskRes'='runModelObj$deskRes'; 'discrim'='runModelObj$discrim'; 'daten'='runModelObj$daten'"), sep =" = ", collapse = ", "), ")",sep="")
                    eval(parse(text=do))                                        ### obere Zeile: baue Aufruf zusammen; rufe 'getConquestResults' mit seinen eigenen Argumenten auf
                    dir <- runModelObj[["dir"]]                                 ### wo Argumente neu vergeben werden, geschieht das in dem 'recode'-Befehl; so wird als 'path'-
                    name<- runModelObj[["analysis.name"]]                       ### Argument 'runModelObj$dir' uebergeben
                    allN<- runModelObj[["all.Names"]]
               }
     ### Software = TAM
               if( isTRUE(attr(runModelObj, "software") == "tam")  || is_true(obj = runModelObj[["software"]], crit = "tam")) {
                    isTa<- TRUE                                                 ### logisches Argument: wurde mit Tam gerechnet?
                    if ( Q3 ) {
                        if ( ncol ( attr(runModelObj, "defineModelObj")[["qMatrix"]]) !=2 ) {
                            cat("Q3 is only available for unidimensional models. Estimation will be skipped.\n")
                            Q3 <- FALSE
                        }
                    }                                                           ### untere zeile: hier wird ggf. die Anzahl der zu ziehenden PVs ueberschrieben
                    if(!is.null(nplausible)) { attr(runModelObj, "defineModelObj")[["n.plausible"]] <- nplausible }  else  { nplausible <- attr(runModelObj, "defineModelObj")[["n.plausible"]] }
                    do  <- paste ( "res <- getTamResults ( ", paste(names(formals(getTamResults)), car::recode(names(formals(getTamResults)),"'pvMethod'='attr(runModelObj, \"defineModelObj\")[[\"pvMethod\"]]'"),  sep =" = ", collapse = ", "), ")",sep="")
                    eval(parse(text=do))
                    dir <- attr(runModelObj, "defineModelObj")[["dir"]]         ### untere zeilen(n): tam summary ergaenzen, wenn mit tam gerechnet wurde
                    name<- attr(runModelObj, "defineModelObj")[["analysis.name"]]## wird als Attribut in Ergebnisstruktur angehangen (nicht so superclever;
                    allN<- attr(runModelObj, "defineModelObj")[["all.Names"]]   ### schoener waers, man wuerde das direkt in die Ergebnisstruktur einbauen)
               }
     ### Software = MIRT
               if(isTRUE(attr(runModelObj, "software") == "mirt")  || is_true(obj = runModelObj[["software"]], crit = "mirt") ) {
                    do  <- paste ( "res <- getMirtResults ( ", paste(names(formals(getMirtResults)),names(formals(getMirtResults)),  sep =" = ", collapse = ", "), ")",sep="")
                    eval(parse(text=do))
                    dir <- attr(runModelObj, "defineModelObj")[["dir"]]
                    allN<- attr(runModelObj, "defineModelObj")[["allNam"]]
                    name<- attr(runModelObj, "defineModelObj")[["analysis.name"]]
               }                                                                ### wenn es ein Rueckgabeobjekt gibt, wird das jetzt um einige technische Eintraege erweitert (geht noch nicht fuer mirt)
               if(!is.null(res) && isFALSE(attr(runModelObj, "software") == "mirt") ) {
                    stopifnot ( length(unique(res[,"model"])) == 1)             ### (das was frueher ueber Attribute gemacht wurde, aber das ist zu krass schlimm beschissen scheiss untransparent, das muss weg!!)
     ### Rueckgabeobjekt mit technischen Parametern 'anreichern', first: 'all.Names' ... auch hier muss die Itemnamenrueckbenennung im falle von conquest mit mehr als 11 Zeichen im Variablennamen stattfinden
                    if(inherits(runModelObj, "runMultiple")) {
                       renam <- do.call("rbind", lapply(runModelObj, FUN = function(y) {y[["renam"]]}))
                    } else {
                       renam <- runModelObj[["renam"]]
                    }
                    alln<- do.call("rbind", lapply(names(allN), FUN = function ( x ) {
                           if(length(allN[[x]]) > 0) {
                              if(x == "variablen" && !is.null(renam)) {
                                 dp <- recodeLookup(allN[[x]], renam[,c("new", "old")])
                              } else {
                                 dp <- allN[[x]]
                              }
                              res <- data.frame ( type = "tech", par = x, derived.par = dp)
                           } else {
                              res <- NULL
                           }
                           return(res)}))
                    res <- plyr::rbind.fill ( res, data.frame ( res[1,c("model", "source")], alln, stringsAsFactors = FALSE) ) |> suppressWarnings()
     ### Rueckgabeobjekt mit technischen Parametern 'anreichern', second: 'dif.setting'
                    difS<- list (abs.dif.bound = abs.dif.bound, sig.dif.bound = sig.dif.bound, p.value = p.value)
                    resD<- data.frame ( res[1,c("model", "source")], type = "tech", par = "dif", derived.par = names(difS), value = unlist(difS), stringsAsFactors = FALSE) |> suppressWarnings()
     ### Rueckgabeobjekt mit technischen Parametern 'anreichern', third: nodes, etc.
                    if ( inherits(runModelObj, "runConquest")) {                ### untere Zeile: bei method muesste 'gauss' etc. eigentlich in die value-Spalte, die wuerde aber dadurch zu character, obwohl sie numerisch sein soll. Deshalb steht da jetzt 1, und 'gauss' in 'derived.par'
                         resN<- data.frame ( res[1,c("model", "source")], type = "tech", par = c("method",rep("nodes", 3)), derived.par = c(runModelObj[["method"]],"nodes", "p.nodes", "f.nodes"), value = c(1,runModelObj[["nodes"]], runModelObj[["p.nodes"]], runModelObj[["f.nodes"]]), stringsAsFactors = FALSE) |> suppressWarnings()
                    }  else  {                                                  ### untere Zeile: Anzahl der nodes in TAM
                         nNod<- length(attr(runModelObj, "defineModelObj")[["control"]][["nodes"]])
                         resN<- data.frame ( res[1,c("model", "source")], type = "tech", par = c("QMC",rep("nodes", 1+nNod)), derived.par = c(NA,rep("discrete.theta",nNod), "snodes"), value = c(attr(runModelObj, "defineModelObj")[["control"]][["QMC"]],attr(runModelObj, "defineModelObj")[["control"]][["nodes"]], attr(runModelObj, "defineModelObj")[["control"]][["snodes"]]), stringsAsFactors = FALSE) |> suppressWarnings()
                    }
                    res <- plyr::rbind.fill ( res, resD, resN )
     ### jetzt wird der ganze scheiss bei Bedarf (wenn der user das will) noch auf der Festplatte gespeichert
                    id  <- unique(res[intersect(which(res[,"type"] == "tech"), which(res[,"par"] == "ID")),"derived.par"])
                    if(!is.null(dir)) {
                        stopifnot(length(id)==1)
                        item<-itemFromRes ( res )
                        if ( file.exists(file.path(dir, paste(name, "_items.csv",sep=""))) & overwrite == FALSE) {
                             cat(paste("Item results cannot be saved, file '",  file.path(dir, paste(name, "_items.csv",sep="")),"' already exists.\n    Please remove/rename existing file or use 'overwrite=TRUE'.\n",sep=""))
                        }  else  {
                             write.csv2(item, file.path(dir, paste(name, "_items.csv",sep="")), na="", row.names = FALSE)
                        }                                                       ### untere Zeilen: speichere wunschgemaess alle Personenparameter in einer Tabelle im Wideformat
                        txt <- capture.output ( wle <- wleFromRes(res) )        ### 'capture.output' wird benutzt um Warnungen in wleFromRes() zu unterdruecken
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
                        }                                                       ### Hier wird geprueft, welche Personenparameter vorliegen
                        alls<- list ( wle, pv, eap )                            ### wenn es Personenparameter gibt, werden sie eingelesen
                        allP<- NULL                                             ### alle vorhandenen Personenparameter werden zum Speichern in einen gemeinsamen Dataframe gemergt
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
                              q3m <- q3FromRes ( res )                          ### q3FromRes liest so viele q3-tabellen aus, wie es gibt, hier will ich nur die erste
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
               attr(res, "runModelAttributes") <- attributes(runModelObj)
               return(res)
               }}
