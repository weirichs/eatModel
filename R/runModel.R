### runModel is usually used after defineModel().

# dient dazu, conquest ueber wine aus linux aufzurufen 
run_conquest_linux <- function(conquest_exe, syntax_file, wait) {
  script <- tempfile(fileext = ".sh")
  writeLines(c("#!/bin/bash",  paste("wine", shQuote(conquest_exe), shQuote(basename(syntax_file))) ), script)
  Sys.chmod(script, "0755")
  system2("xfce4-terminal",
          args = c("--disable-server",  paste0("--working-directory=", shQuote(dirname(syntax_file))),
                    paste0("--command=", shQuote(script))),  wait = wait)}

### ueberschreibe 'originales' anchor objekt! (hilfsfunktion fuer runmodel mit tam)
overwriteAnchorGpcmTAM <- function(anchor, dmo, Y, group, wgt) {
   if(!is.null(anchor) && grepl("pcm", dmo[["irtmodel"]], ignore.case=TRUE)) {
      beg    <- Sys.time()
      control<- dmo[["control"]]
      control[["maxiter"]] <- 50
      if(dmo[["irtmodel"]] %in% c("1PL", "PCM", "PCM2", "RSM")) {
         if(dmo[["irtmodel"]] == "PCM") {
            skelet <- generatePCM_skeleton_TAM(dmo=dmo, anchor=anchor)
            message(paste0("Model '",dmo[["analysis.name"]],"': Reconstruct the 'xsi.fixed.estimated' structure to prepare partial credit anchoring in TAM."))   
         } else {    
            skelet <- tam.mml(resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], constraint = dmo[["constraint"]], pid = dmo[["daten"]][,"ID"], Y = Y, Q = dmo[["qMatrix"]][,-1,drop=FALSE], irtmodel = dmo[["irtmodel"]], pweights = wgt, control = control, group=group)
         }
      } else {
         skelet <- tam.mml.2pl(resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], pid = dmo[["daten"]][,"ID"], Y = Y, Q = dmo[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = anchor, irtmodel = dmo[["irtmodel"]], est.slopegroups=dmo[["est.slopegroups"]],pweights = wgt, B.fixed = dmo[["fixSlopeMat"]], est.variance = dmo[["estVar"]], control = control, group=group)
      }         
      diffe  <- Sys.time() - beg
      if(as.numeric(diffe) > 0.2 && dmo[["irtmodel"]] != "PCM") {message(paste0("Model '",dmo[["analysis.name"]],"': Generate skeleton for partial credit anchoring: ", timeFormat(diffe)))}
      anchor <- prepAnchorTAM(dfm = dmo, skeleton = skelet[["xsi.fixed.estimated"]])
   }
   return(anchor)}
   
### die Funktion erzeugt die Liste, die normalerweise das gefittete tam-objekt mit mod[["xsi.fixed.estimated"]] zurueck gibt, nur mit leeren Werten
### in prepAnchorTAM werden hierin dann die VORHANDENEN anker-Werte eingetragen. Items, die nicht verankert werden, werden aus "xsi.fixed.estimated" entfernt
### diese Funktion erzeugt auch die NUmmerierung, die wichtig ist und nicht geaendert werden darf 
generatePCM_skeleton_TAM <- function(dmo, anchor){
   itemSeq <- do.call("rbind", lapply(dmo[["all.Names"]][["variablen"]], FUN = function(i) {
              cats <- paste0("Cat", 1:max(dmo[["daten"]][,i], na.rm=TRUE))
              rows <- paste(i, cats, sep="_")
              dfr  <- data.frame(empty=rep(NA,length(rows)),  xsi = 1)
              rownames(dfr) <- rows
              return(dfr)})) 
   itemSeq[,"empty"] <- 1:nrow(itemSeq)
   colnames(itemSeq) <- c("", "xsi")
   itemSeq <- as.matrix(itemSeq)
   ret     <- list(xsi.fixed.estimated = itemSeq)
   return(ret)}   

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
    ### runModel for conquest, mirt, tam
                if(inherits(defineModelObj, "defineConquest")) {ret <- runModelConquest(dmo=defineModelObj, show.dos.console=show.dos.console,show.output.on.console=show.output.on.console, wait=wait)}   
                if(inherits(defineModelObj, "defineMirt"))     {ret <- runModelMirt(dmo=defineModelObj)}
                if(inherits(defineModelObj, "defineTam"))      {ret <- runModelTam(dmo=defineModelObj, show.output.on.console=show.output.on.console)}
            }
            return(ret) }

runModelConquest <- function(dmo, show.dos.console,show.output.on.console, wait) {
    sysInfo  <- Sys.info()
    if(sysInfo[["sysname"]] == "Linux") {
       foo     <- run_conquest_linux(conquest_exe = dmo$conquest.folder, syntax_file = dmo$input, wait = wait) 
    } else {
       oldPfad <- getwd()
       setwd(dmo$dir)
       suppressWarnings(system(paste(dmo$conquest.folder," ",dmo$input,sep=""),invisible=!show.dos.console,show.output.on.console=show.output.on.console, wait=wait) )
       if(wait == FALSE) { Sys.sleep(0.2) }
       setwd(oldPfad)                                                           ### untere Zeile: Rueckgabeobjekt definieren: Conquest
    }
    class(dmo) <- c("runConquest", "list")
    attr(dmo, "software") <- "conquest"
    return(dmo) }      

runModelMirt <- function(dmo) {
    if(ncol(dmo[["qMatrix"]]) == 2) {                                           ### eindimensionales Modell
       mirtMod <- 1
    } else {                                                                    ### mehrdimensionales Modell
       mirtMod<- mirt.model(as.matrix(dmo[["qMatrix"]][,-1]), COV= matrix(rep(TRUE, times = 2*(ncol(dmo[["qMatrix"]])-1)), ncol(dmo[["qMatrix"]])-1))
    }                                                            
    if(!is.null(dmo[["allNam"]][["HG.var"]])) {
       covdata<- dmo[["daten"]][,dmo[["allNam"]][["HG.var"]], drop=FALSE]
       formula<- as.formula(paste0("~ ", paste(dmo[["allNam"]][["HG.var"]], collapse = " + ")))
    } else {
       covdata<- NULL; formula <- NULL
    }                                                                           ### wenn untere Zeile TRUE, dann wird skeleton gebraucht
    if(!is.null(dmo[["anchor"]][["ank"]]) || !is.null(dmo[["fixSlopeMat"]][["slopMat"]])) {pars <- "values"} else {pars <- NULL}  
    if("Rasch" %in%  dmo[["irtmodel"]][,2]) {pars <- "values"}
    if(!is.null(dmo[["allNam"]][["weight.var"]])) {
       wgtvec <- dmo[["daten"]][,dmo[["allNam"]][["weight.var"]]]
    } else {
       wgtvec <- NULL
    }
    tech <- dmo[["technical"]]
    skel <- mirt(data = dmo[["daten"]][,dmo[["allNam"]][["variablen"]]], model = mirtMod,  itemtype = dmo[["irtmodel"]][,2], SE = TRUE,  covdata=covdata, formula=formula, verbose =dmo[["progress"]], pars="values", method = dmo[["met"]][["method"]], quadpts = dmo[["met"]][["nodes"]], survey.weights = wgtvec, technical=tech)
    if(!is.null(pars)) {                                                        ### constraints werden in adaptSkelForAnchor() umgesetzt 
       if(isFALSE(onlySkeleton)) {message("Modify skeleton ... ")}              ### skeleton anpassen
       skelN<- skel <- adaptSkelForAnchor(allNam = dmo[["allNam"]], skel = skel, anch = dmo[["anchor"]], qmat = dmo[["qMatrix"]], slope = dmo[["fixSlopeMat"]], irtmodel =  dmo[["irtmodel"]], est.slopegroups =dmo[["est.slopegroups"]][["esg"]])
    } else {                                                                    ### wenn skeleton angepasst wurde, soll der angepasste skeleton als attribut gespeichert werden
       skelN<- NULL                                                             ### wenn er NICHT angepasst wurde, soll der originale (d.h., der von mirt erzeugte)
    }                                                                           ### skeleton als attribut gespeichert werden 
    if(isTRUE(onlySkeleton)) {
       return(skel)
    } else {
       mod  <- mirt(data = dmo[["daten"]][,dmo[["allNam"]][["variablen"]]], model = mirtMod,  itemtype = dmo[["irtmodel"]][,2], SE = TRUE,  covdata=covdata, formula=formula, verbose =dmo[["progress"]], pars=skelN, method = dmo[["met"]][["method"]], quadpts = dmo[["met"]][["nodes"]], survey.weights = wgtvec, technical=tech)
       attr(mod, "defineModelObj") <- dmo[-match("daten", names(dmo))]
       attr(mod, "personID") <- dmo[["daten"]][,"ID"]
       attr(mod, "software") <- "mirt"
       attr(mod, "skeleton") <- skel
       return(mod)
    }}   

runModelTam <- function(dmo, show.output.on.console) {
    if(show.output.on.console == TRUE) {dmo[["control"]][["progress"]] <- TRUE }
    if(length( dmo[["all.Names"]][["HG.var"]])>0)     { Y <- dmo[["daten"]][,dmo[["all.Names"]][["HG.var"]], drop=FALSE] } else { Y <- NULL }
    if(length( dmo[["all.Names"]][["weight.var"]])>0) { wgt <- as.vector(dmo[["daten"]][,dmo[["all.Names"]][["weight.var"]]])} else {wgt <- NULL}
    if(length( dmo[["all.Names"]][["group.var"]])>0)  { group <- as.vector(dmo[["daten"]][,dmo[["all.Names"]][["group.var"]]])} else {group <- NULL}
    stopifnot(all(dmo[["qMatrix"]][,1] == dmo[["all.Names"]][["variablen"]]))
    ### Achtung! in alter Paketversion wurde der anchor parameter frame noch in 'defineModel' fuer TAM aufbereitet, neuerdings in 'runModel'.
    ### Grund: fuer partial credit muss sicherheitshalber erst ein 'skeleton' erzeugt werden, damit die richtigkeit der reihenfolge der
    ### Verankerungsparameter sichergestellt ist! Da der anchor parameter frame fuer mehrere Modelle gebraucht wird, aber nur fuer partial credit
    ### mittels skeleton erzeugt werden muss, geschieht das hier zweimal, erstmal allgemein (untere Zeile); fuer partial credit wird das dann nochmal ueberschrieben
    anchor <- prepAnchorTAM(dfm = dmo)
    if(length(dmo[["all.Names"]][["DIF.var"]]) == 0) {
       anchor <- overwriteAnchorGpcmTAM(anchor=anchor, dmo=dmo, Y=Y, group=group, wgt=wgt)
       if(dmo[["irtmodel"]] %in% c("1PL", "PCM", "PCM2", "RSM")) {
          if(isTRUE(dmo[["fitTamMmlForBayesian"]])) {
             mod  <- tam.mml(resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], constraint = dmo[["constraint"]], pid = dmo[["daten"]][,"ID"], Y = Y, Q = dmo[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = anchor, irtmodel = dmo[["irtmodel"]], pweights = wgt, control = dmo[["control"]], group=group)
          } else {
             mod  <- tamObjForBayesianPV (anchor = anchor, qMatrix = dmo[["qMatrix"]], resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], pid = dmo[["daten"]][,"ID"], Y=Y)
          }
       }
       if(dmo[["irtmodel"]] %in% c("2PL", "GPCM", "GPCM.groups", "2PL.groups", "GPCM.design", "3PL"))  {
          if(dmo[["irtmodel"]] == "3PL") {
             if(isTRUE(dmo[["fitTamMmlForBayesian"]]) ) {
                mod  <- tam.mml.3pl(resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], pid = dmo[["daten"]][,"ID"], Y = Y, Q = dmo[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = anchor, pweights = wgt, est.guess =dmo[["guessMat"]],  est.variance = dmo[["estVar"]], control = dmo[["control"]], group=group)
             } else {
                mod  <- tamObjForBayesianPV (anchor = anchor, qMatrix = dmo[["qMatrix"]], resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], pid = dmo[["daten"]][,"ID"], Y=Y, slopeMatrix = dmo[["fixSlopeMat"]])
             }
          } else{
             if(dmo[["fitTamMmlForBayesian"]] == TRUE) {
                mod  <- tam.mml.2pl(resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], pid = dmo[["daten"]][,"ID"], Y = Y, Q = dmo[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = anchor, irtmodel = dmo[["irtmodel"]], est.slopegroups=dmo[["est.slopegroups"]],pweights = wgt, B.fixed = dmo[["fixSlopeMat"]], est.variance = dmo[["estVar"]], control = dmo[["control"]], group=group)
             } else {
                mod  <- tamObjForBayesianPV (anchor = anchor, qMatrix = dmo[["qMatrix"]], resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], pid = dmo[["daten"]][,"ID"], Y=Y, slopeMatrix = dmo[["fixSlopeMat"]])
             }
          }
       }
    } else {
       assign(paste("DIF_",dmo[["all.Names"]][["DIF.var"]],sep="") , as.data.frame (dmo[["daten"]][,dmo[["all.Names"]][["DIF.var"]]]) )
       formel   <- as.formula(paste("~item - ",paste("DIF_",dmo[["all.Names"]][["DIF.var"]],sep="")," + item * ",paste("DIF_",dmo[["all.Names"]][["DIF.var"]],sep=""),sep=""))
       if(grepl("PCM", dmo[["irtmodel"]])) {
          formel <- as.formula(paste0("~item+item:step + ",paste("DIF_",dmo[["all.Names"]][["DIF.var"]],sep="")," * item*step"))
       }
       facetten <- as.data.frame (dmo[["daten"]][,dmo[["all.Names"]][["DIF.var"]]])
       colnames(facetten) <- paste("DIF_",dmo[["all.Names"]][["DIF.var"]],sep="")
       if(isTRUE(dmo[["fitTamMmlForBayesian"]]) ) {
          if(grepl("PCM", dmo[["irtmodel"]])) {
             mod  <- tam.mml.mfr(resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], facets = facetten, formulaA = formel, pid = dmo[["daten"]][,"ID"], control = dmo[["control"]], group=group)
          } else {
             mod  <- tam.mml.mfr(resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], facets = facetten, constraint = dmo[["constraint"]], formulaA = formel, pid = dmo[["daten"]][,"ID"], Y = Y, Q = dmo[["qMatrix"]][,-1,drop=FALSE], xsi.fixed = anchor, irtmodel = dmo[["irtmodel"]], pweights = wgt, control = dmo[["control"]], group=group)
          }
       } else {
          mod  <- tamObjForBayesianPV (anchor = anchor, qMatrix = dmo[["qMatrix"]], resp = dmo[["daten"]][,dmo[["all.Names"]][["variablen"]]], pid = dmo[["daten"]][,"ID"], Y=Y, slopeMatrix = dmo[["fixSlopeMat"]])
       }                                                                        ### hier werden fuer 'tam' zusaetzliche Objekte als Attribute an das Rueckgabeobjekt angehangen
    }                                                                           ### Grund: Rueckgabeobjekt soll weitgehend beibehalten werden, damit alle 'tam'-Funktionen, die darauf aufsetzen, lauffaehig sind
    attr(mod, "defineModelObj") <- dmo[-match("daten", names(dmo))]
    attr(mod, "Y")              <- Y
    attr(mod, "software")       <- "tam"
    return(mod)}

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





