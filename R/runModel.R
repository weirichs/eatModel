### runModel is usually used after defineModel().

runModel <- function(defineModelObj, show.output.on.console = FALSE, show.dos.console = TRUE, wait = TRUE, onlySkeleton = FALSE) {
            argString <-  paste(names(formals(runModel))[-1], names(formals(runModel))[-1], sep=" = ", collapse=", ")
            if (inherits(defineModelObj, "defineMultiple") ) {                  ### erstmal fuer den Multimodellfall: nur dafuer wird single core und multicore unterschieden
                if(is.null ( attr(defineModelObj, "split")[["nCores"]] ) || attr(defineModelObj, "split")[["nCores"]] == 1 ) {
                   res <- lapply(defineModelObj, FUN = function ( r ) {         ### erstmal: single core
                          do  <- paste("runModel ( defineModelObj = r, ",argString, ")",sep="")
                          ret <- eval(parse(text=do))
                          return(ret)})
                }  else  {                                                      ### multicore
                   doIt<- function (laufnummer,  ... ) {
                          if(!"eatModel" %in% .packages()) {library(eatModel)}
                          if(!"TAM" %in% .packages()) {library(TAM)}
                          if(!"mirt" %in% .packages()) {library(mirt)}
                          do  <- paste("runModel ( defineModelObj = defineModelObj[[laufnummer]], ",argString, ")",sep="") 
                          ret <- eval(parse(text=do))
                          return(ret) }
                   beg <- Sys.time()
                   if ( attr(defineModelObj, "split")[["mcPackage"]] == "parallel") {
                        cl  <- makeCluster(attr(defineModelObj, "split")[["nCores"]], type = "SOCK")
                   }  else  {
                        cl  <- future::makeClusterPSOCK(attr(defineModelObj, "split")[["nCores"]], verbose=FALSE)
                   }
                   do  <- paste("clusterApply(cl = cl, x = 1:length(defineModelObj), fun = doIt , ",argString, ")",sep="") 
                   res <- eval(parse(text=do))
                   stopCluster(cl)
                   cat(paste ( length(defineModelObj), " analyses finished: ", sep="")); print( Sys.time() - beg, digits=3)
                }
                class(res) <- c("runMultiple", "list")
                attr(res, "split") <- attr(defineModelObj, "split")
                return(res)
            } else {                                                            ### ab hier fuer den single model Fall
    ### runModel for conquest
                if(inherits(defineModelObj, "defineConquest")) {
                   oldPfad <- getwd()
                   setwd(defineModelObj$dir)
                   suppressWarnings(system(paste(defineModelObj$conquest.folder," ",defineModelObj$input,sep=""),invisible=!show.dos.console,show.output.on.console=show.output.on.console, wait=wait) )
                   if(wait == FALSE) { Sys.sleep(0.2) }
                   setwd(oldPfad)                                               ### untere Zeile: Rueckgabeobjekt definieren: Conquest
                   class(defineModelObj) <- c("runConquest", "list")
                   attr(defineModelObj, "software") <- "conquest"
                   return ( defineModelObj )
                }
                if(inherits(defineModelObj, "defineMirt")) {
    ### runModel for mirt 
                   if(ncol(defineModelObj[["qMatrix"]]) == 2) {                 ### eindimensionales Modell
                      mirtMod <- 1
                   } else {                                                     ### mehrdimensionales Modell
                      mirtMod<- mirt.model(as.matrix(defineModelObj[["qMatrix"]][,-1]), COV= matrix(rep(TRUE, times = 2*(ncol(defineModelObj[["qMatrix"]])-1)), ncol(defineModelObj[["qMatrix"]])-1))
                   }                                                            
                   if(!is.null(defineModelObj[["allNam"]][["HG.var"]])) {
                      covdata<- defineModelObj[["daten"]][,defineModelObj[["allNam"]][["HG.var"]], drop=FALSE]
                      formula<- as.formula(paste0("~ ", paste(defineModelObj[["allNam"]][["HG.var"]], collapse = " + ")))
                   } else {
                      covdata<- NULL; formula <- NULL
                   }                                                            ### wenn untere Zeile TRUE, dann wird skeleton gebraucht
                   if(!is.null(defineModelObj[["anchor"]][["ank"]]) || !is.null(defineModelObj[["fixSlopeMat"]][["slopMat"]])) {pars <- "values"} else {pars <- NULL}  
                   if("Rasch" %in%  defineModelObj[["irtmodel"]][,"irtmod"]) {pars <- "values"}
                   skel <- mirt(data = defineModelObj[["daten"]][,defineModelObj[["allNam"]][["variablen"]]], model = mirtMod,  itemtype = defineModelObj[["irtmodel"]][,"irtmod"], SE = TRUE,  covdata=covdata, formula=formula, verbose =defineModelObj[["progress"]], pars="values")
                   # skel[grep("^d", skel[,"name"]),"est"] <- TRUE              ### standardmaessig alle Itemschwierigkeiten frei schaetzen lassen: problematisch
                   if(!is.null(pars)) {                                         ### constraints werden in adaptSkelForAnchor() umgesetzt 
                      if(isFALSE(onlySkeleton)) {message("Modify skeleton ... ")}## skeleton anpassen 
                      skelN<- skel <- adaptSkelForAnchor(allNam = defineModelObj[["allNam"]], skel = skel, anch = defineModelObj[["anchor"]], qmat = defineModelObj[["qMatrix"]], slope = defineModelObj[["fixSlopeMat"]], irtmodel =  defineModelObj[["irtmodel"]], est.slopegroups =defineModelObj[["est.slopegroups"]][["esg"]])
                   } else {                                                     ### wenn skeleton angepasst wurde, soll der angepasste skeleton als attribut gespeichert werden
                      skelN<- NULL                                              ### wenn er NICHT angepasst wurde, soll der originale (d.h., der von mirt erzeugte)
                   }                                                            ### skeleton als attribut gespeichert werden 
                   if(isTRUE(onlySkeleton)) {
                      return(skel)
                   } else {
                      mod  <- mirt(data = defineModelObj[["daten"]][,defineModelObj[["allNam"]][["variablen"]]], model = mirtMod,  itemtype = defineModelObj[["irtmodel"]][,"irtmod"], SE = TRUE,  covdata=covdata, formula=formula, verbose =defineModelObj[["progress"]], pars=skelN)
                      attr(mod, "defineModelObj") <- defineModelObj[-match("daten", names(defineModelObj))]
                      attr(mod, "personID") <- defineModelObj[["daten"]][,"ID"]
                      attr(mod, "software") <- "mirt"
                      attr(mod, "skeleton") <- skel
                      return(mod)
                   }
                }
    ### runModel for TAM
                if(inherits(defineModelObj, "defineTam")) {
                   if ( show.output.on.console == TRUE ) { defineModelObj[["control"]][["progress"]] <- TRUE }
                   if(length( defineModelObj[["all.Names"]][["HG.var"]])>0)     { Y <- defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["HG.var"]], drop=FALSE] } else { Y <- NULL }
                   if(length( defineModelObj[["all.Names"]][["weight.var"]])>0) { wgt <- as.vector(defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["weight.var"]]])} else {wgt <- NULL}
                   if(length( defineModelObj[["all.Names"]][["group.var"]])>0)  { group <- as.vector(defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["group.var"]]])} else {group <- NULL}
                   stopifnot(all(defineModelObj[["qMatrix"]][,1] == defineModelObj[["all.Names"]][["variablen"]]))
    ### Achtung! in alter Paketversion wurde der anchor parameter frame noch in 'defineModel' fuer TAM aufbereitet, neuerdings in 'runModel'.
    ### Grund: fuer partial credit muss sicherheitshalber erst ein 'skeleton' erzeugt werden, damit die richtigkeit der reihenfolge der
    ### Verankerungsparameter sichergestellt ist! Da der anchor partameter frame fuer mehrere Modelle gebraucht wird, aber nur fuer partial credit
    ### mittels skeleton erzeugt werden muss, geschieht das hier zweimal, erstmal allgemein (untere Zeile); fuer partial credit wird das dann nochmal ueberschrieben
                   anchor <- prepAnchorTAM(dfm = defineModelObj)
                   if(length(defineModelObj[["all.Names"]][["DIF.var"]]) == 0 ) {
                      if( defineModelObj[["irtmodel"]] %in% c("1PL", "PCM", "PCM2", "RSM", "GPCM", "GPCM.groups")) {
                          if ( isTRUE(defineModelObj[["fitTamMmlForBayesian"]]) ) {
    ### ueberschreibe 'originales' anchor objekt!
                               if(!is.null(anchor) && grepl("pcm", defineModelObj[["irtmodel"]], ignore.case=TRUE)) {
                                  beg    <- Sys.time()
                                  control<- defineModelObj[["control"]]
                                  control[["maxiter"]] <- 50
                                  skelet <- tam.mml(resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], constraint = defineModelObj[["constraint"]], pid = defineModelObj[["daten"]][,"ID"], Y = Y, Q = defineModelObj[["qMatrix"]][,-1,drop=FALSE], irtmodel = defineModelObj[["irtmodel"]], pweights = wgt, control = control, group=group)
                                  diffe  <- Sys.time() - beg
                                  if(as.numeric(diffe) > 0.2) {message(paste0("Generate skeleton for partial credit anchoring: ", timeFormat(diffe)))}
                                  anchor <- prepAnchorTAM(dfm = defineModelObj, skeleton = skelet[["xsi.fixed.estimated"]])
                               }
                               mod  <- tam.mml(resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], constraint = defineModelObj[["constraint"]], pid = defineModelObj[["daten"]][,"ID"], Y = Y, Q = defineModelObj[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = anchor, irtmodel = defineModelObj[["irtmodel"]], pweights = wgt, control = defineModelObj[["control"]], group=group)
                          }  else  {
                               mod  <- tamObjForBayesianPV (anchor = anchor, qMatrix = defineModelObj[["qMatrix"]], resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y=Y)
                          }
                      }
                      if( defineModelObj[["irtmodel"]] %in% c("2PL", "GPCM", "GPCM.groups", "2PL.groups", "GPCM.design", "3PL") )  {
                          if( defineModelObj[["irtmodel"]] == "3PL") {
                              if ( isTRUE(defineModelObj[["fitTamMmlForBayesian"]]) ) {
                                   mod  <- tam.mml.3pl(resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y = Y, Q = defineModelObj[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = anchor, pweights = wgt, est.guess =defineModelObj[["guessMat"]],  est.variance = defineModelObj[["estVar"]], control = defineModelObj[["control"]], group=group)
                              }  else  {
                                   mod  <- tamObjForBayesianPV (anchor = anchor, qMatrix = defineModelObj[["qMatrix"]], resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y=Y, slopeMatrix = defineModelObj[["fixSlopeMat"]])
                              }
                          }  else {
                              if ( defineModelObj[["fitTamMmlForBayesian"]] == TRUE ) {
                                   mod  <- tam.mml.2pl(resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y = Y, Q = defineModelObj[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = anchor, irtmodel = defineModelObj[["irtmodel"]], est.slopegroups=defineModelObj[["est.slopegroups"]],pweights = wgt, B.fixed = defineModelObj[["fixSlopeMat"]], est.variance = defineModelObj[["estVar"]], control = defineModelObj[["control"]], group=group)
                              }  else  {
                                   mod  <- tamObjForBayesianPV (anchor = anchor, qMatrix = defineModelObj[["qMatrix"]], resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y=Y, slopeMatrix = defineModelObj[["fixSlopeMat"]])
                              }
                          }
                      }
                   } else {
                     assign(paste("DIF_",defineModelObj[["all.Names"]][["DIF.var"]],sep="") , as.data.frame (defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["DIF.var"]]]) )
                     formel   <- as.formula(paste("~item - ",paste("DIF_",defineModelObj[["all.Names"]][["DIF.var"]],sep="")," + item * ",paste("DIF_",defineModelObj[["all.Names"]][["DIF.var"]],sep=""),sep=""))
                     facetten <- as.data.frame (defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["DIF.var"]]])
                     colnames(facetten) <- paste("DIF_",defineModelObj[["all.Names"]][["DIF.var"]],sep="")
                     if ( isTRUE(defineModelObj[["fitTamMmlForBayesian"]]) ) {
                          mod  <- tam.mml.mfr(resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], facets = facetten, constraint = defineModelObj[["constraint"]], formulaA = formel, pid = defineModelObj[["daten"]][,"ID"], Y = Y, Q = defineModelObj[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = anchor, irtmodel = defineModelObj[["irtmodel"]], pweights = wgt, control = defineModelObj[["control"]], group=group)
                     }  else  {
                          mod  <- tamObjForBayesianPV (anchor = anchor, qMatrix = defineModelObj[["qMatrix"]], resp = defineModelObj[["daten"]][,defineModelObj[["all.Names"]][["variablen"]]], pid = defineModelObj[["daten"]][,"ID"], Y=Y, slopeMatrix = defineModelObj[["fixSlopeMat"]])
                     }                                                          ### hier werden fuer 'tam' zusaetzliche Objekte als Attribute an das Rueckgabeobjekt angehangen
                   }                                                            ### Grund: Rueckgabeobjekt soll weitgehend beibehalten werden, damit alle 'tam'-Funktionen, die darauf aufsetzen, lauffaehig sind
                   attr(mod, "defineModelObj") <- defineModelObj[-match("daten", names(defineModelObj))]
                   attr(mod, "Y")              <- Y
                   attr(mod, "software")       <- "tam"
                   return(mod)  }  }   }


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
  class(tamObj) <- c("tamBayes", "list")
  return(tamObj)}





