### Hilfsfunktion: transformiert die Q matrix (ggf. auch fuer 2pl) in das von TAM benoetigte Format
qMatToB <- function(qma, slp) {                                                 ### es muesste immer die erste Spalte aus 'slp' sein, egal wie diese nun genau benannt ist
      zei <- match( qma[,"item"], slp[,1])                                      ### reihenfolge aus Q Matrix ist die Refrenz (die stimmt!)
      for ( i in 1:length(zei) ) {
           ind <- which(qma[i,] ==1 )
           stopifnot(length(ind)==1, qma[i,"item"] == slp[zei[i],1])            ### zwei bedingungen simultan pruefen
           qma[i,ind] <- slp[zei[i],2] }
      return(qma)}

tamObjForBayesianPV <- function(anchor, qMatrix, slopeMatrix = NULL, resp, pid, Y) {
      warning("To date, bayesian plausible values imputation only works for binary between-item dimensionality models.")
      if ( !is.null(slopeMatrix)) {                                             ### Ladungen in Q-Matrix schreiben
           qMatrix <- qMatToB ( qma = qMatrix, slp = slopeMatrix)
      }
      xsi.obj<- as.matrix(data.frame ( V1 = 0, V2 = anchor[,"parameter"] * (-1)))## zweite Spalte in 'anchor' heisst immer 'parameter' - wird in Funktion 'anker()' so definiert
      B.obj  <- array(unlist(lapply(2:ncol(qMatrix),FUN = function ( col) {data.frame ( Cat0 = 0, Cat1 = qMatrix[,col])})), dim = c(nrow(qMatrix), 2, ncol(qMatrix)-1), dimnames = list(qMatrix[,"item"], c("Cat0", "Cat1"), paste0("Dim0",1:(ncol(qMatrix)-1)) ))
      tamObj <- list ( AXsi = xsi.obj, B = B.obj, resp = resp, Y=Y, pid = pid)
      class(tamObj) <- c("list", "tamBayes")
      return(tamObj)}


### prueft, ob Design verlinkt ist: checkLinking(design[1:15,c(1:5,ncol(design))], bookletColumn = "TH")
checkLinking <- function(design, blocks=NULL, bookletColumn=NULL, verbose=FALSE) {
  if(!all(sapply(design, inherits, what="character"))) design <- data.frame(lapply(design, as.character), stringsAsFactors=FALSE)
  if(!is.null(bookletColumn)) {
    if(length(bookletColumn) != 1) stop("Argument 'bookletColumn' must be of length 1.")
    book  <- eatTools::existsBackgroundVariables(dat = design, variable=bookletColumn)
  } else {
    design$bl <- paste0("T", 1:nrow(design))
    book <- "bl"
  }
  items <- setdiff(colnames(design), book)
  if(!is.null(blocks)) {
    stopifnot(is.vector(blocks))
    stopifnot(length(intersect(unlist(design), blocks)) > 0)
    if(any(grepl("[[:punct:]]", blocks))) if(verbose) message("Special characters and spaces will be removed from names in 'blocks'")
    blocks <- paste0("B", gsubAll(as.character(blocks), old=c(" ", "[[:punct:]]"), new=c("", "")))
  }
  design <-	removeNAd(design, items)
  items <- setdiff(colnames(design), book)
  if(any(grepl("[[:punct:]]", unlist(design[,items])))) if(verbose) message("Special characters and spaces will be removed from block identifiers in 'design'")
  design[,items] <- data.frame(lapply(design[,items], function(g) paste0("B", gsubAll(g, old=c(" ", "[[:punct:]]"), new=c("", "")))))
  if(!is.null(blocks)) {
    desd <- apply(design[,items], 2, function(k) ifelse(k %in% blocks, k, NA))
    design <- cbind(design[,book,drop=FALSE], desd)
    items <- setdiff(colnames(design), book)
    design <-	removeNAd(design, items)
    items <- setdiff(colnames(design), book)
  }
  dat   <- do.call(plyr::rbind.fill, apply(design, MARGIN = 1, FUN = simDat, booklet = book))
  link  <- checkLink(dataFrame = dat[,-1, drop = FALSE], remove.non.responser = TRUE, verbose = TRUE)
  return(link)
}


removeNAd <- function(design, items) {
  weg   <- which(rowSums(do.call("rbind", plyr::alply(design[,items], .margins = 1, .fun = is.na))) == ncol(design[,items]))
  if(length(weg)>0) design <- design[-weg,]
  weg2   <- which(colSums(do.call("rbind", plyr::alply(design, .margins = 1, .fun = is.na))) == nrow(design))
  if(length(weg2)>0) design <- design[,-weg2]
  return(design)
}

### Hilfsfunktion fuer 'checkLinking'
simDat <- function ( z, booklet ) {                                             ### erzeugt Datensatz aus einer Zeile des Designs
          if ( !length(na.omit(z[-match(booklet, names(z))])) == length(unique(na.omit(z[-match(booklet, names(z))] ))) ) { stop("Blocks are not unique in each line.\n")}
          items<- as.vector(sapply(na.omit(z[-match(booklet, names(z))]), FUN= function ( i ) { paste(i, 1:3, sep="_")}))
          pers <- paste(z[[booklet]], 11:22, sep="_")                           ### Funktion muss also ueber "apply" aufgerufen werden!
          mat  <- data.frame ( id = pers, matrix ( sample ( 0:1, size = length(pers) * length(items), replace = TRUE), ncol = length(items), nrow = length(pers)))
          colnames(mat)[-1] <- items
          return(mat)}

### Hilfsfunktion zur Bestimmung der Anzahl der Beobachtungen je Itempaar
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

### anchor            ... data frame with anchor parameters. first column is item identifier, second column parameter value
### mRef              ... numeric value: mean of the reference population
### sdRef             ... numeric value: standard deviation of the reference population
### addConst          ... additive constant for transformation
### multConst         ... multiplicative constant for transformation
### cutScores         ... list of two elements. "values" is a numeric vector of cut scores (increasing),
###                       "labels" is an optional character vector of cut score labels
###                       "labels" has be be of length(values)+1
#simEquiTable <- function ( anchor, mRef, sdRef, addConst = 500, multConst = 100, cutScores , dir , n = 2000, conquest.folder ) {
#                if ( length(which ( duplicated(anchor[,1])))>0) {
#                     cat(paste("Warning: Remove ",length(duplicated(anchor[,1]))," entries in the anchor parameter frame.\n",sep=""))
#                     anchor <- anchor[!duplicated(anchor[,1]),]
#                }
#                it  <- matrix ( data = anchor[,2], ncol = length(anchor[,2]), nrow = n, byrow = TRUE)
#                colnames(it) <- anchor[,1]
#                pop <- melt ( data.frame ( idstud = paste("P",1:n,sep=""), theta = rnorm (n = n, mean = mRef, sd = sdRef), it), id.vars = c("idstud", "theta"), value.name = "itemPar", variable.name = "item")
#                pop[,"resp"] <- item.logit(z = pop[,"theta"]-pop[,"itemPar"], thr = c(0.0), slope = 1)$x
#                popW<- dcast(pop,idstud~item, value.var = "resp")
#                mDef<- defineModel ( dat = popW, items = -1, id = "idstud", anchor = anchor, dir = dir, conquest.folder = conquest.folder, compute.fit = FALSE, analysis.name = "equSimTest")
#                mRun<- runModel(mDef)
#                equ <- get.equ ( file.path ( dir, "equSimTest.equ"))[[1]]
#                equ[,"estBista"] <- (equ[,"Estimate"] - mRef) / sdRef * multConst + addConst
#                equ[,"ks"]       <- num.to.cat ( x = equ[,"estBista"], cut.points = cutScores[["values"]], cat.values = cutScores[["labels"]])
#    ### jetzt noch die shortversion der Aequivalenztabelle erzeugen
#                shrt<- do.call("rbind", by ( data = equ, INDICES = equ[,"ks"], FUN = function ( sks ) {
#                       sks1<- data.frame ( do.call("cbind", lapply ( setdiff ( colnames(sks), "std.error"), FUN = function ( col ) {
#                              if ( length( unique ( sks[,col] )) > 1) {
#    ### hier werden spaltenspezifisch die Nachkommastellen bestimmt, auf die gerundet werden soll
#                                   dig <- as.numeric(recode ( col, "'Score'=0; 'Estimate'=2; 'estBista'=0"))
#                                   ret <- paste ( round(min(sks[,col]), digits = dig), round(max(sks[,col]), digits = dig), sep=" bis ")
#                              }  else  {
#                                   ret <- unique ( sks[,col] )
#                              }
#                              return(ret)})) )
#                       colnames(sks1) <- setdiff ( colnames(sks), "std.error")
#                       return(sks1)}))
#                return(list ( complete = equ, short = shrt))}
### Test:
### ret <- simEquiTable( anchor = data.frame ( item = paste("i",1:20,sep=""), par = rnorm(20, mean = -.1, sd = 1.5)), mRef = -0.05, sdRef = 0.9, cutScores = list ( values = 330+0:4*75, labels = c("1a", "1b", 2:5) ), dir = "c:/users/weirichs/test", conquest.folder = "N:/console_Feb2007.exe")

### neue Version derselben Funktion
simEquiTable <- function ( anchor, mRef, sdRef, addConst = 500, multConst = 100, cutScores) {
                anchor<- eatTools::makeDataFrame(anchor)
    ### various checks ...
                if ( ncol(anchor) != 2) {
                     warning(paste0("'anchor' has ",ncol(anchor)," columns. First column is used as item ID, second column is used as item parameter."))
                }
                if(!inherits(anchor[,2], c("integer", "numeric"))) {stop("Item parameter column must be numeric.")}
                if(length(unique(anchor[,1])) != nrow(anchor)) {stop("Item ID column has duplicated entries.")}
    ### temporaeren Datensatz mit allen moeglichen Summenscores erzeugen
                dtmp  <- data.frame(rbind(1*(lower.tri(matrix(1, nrow = nrow(anchor), ncol = nrow(anchor)))),1))
                dtmp  <- data.frame(dtmp, score = rowSums(dtmp) , irtoys::wle(dtmp, cbind(1, anchor[,2], 0)), stringsAsFactors = FALSE)
                dtmp[,"bista"] <- (dtmp[,"est"] - mRef) / sdRef * multConst + addConst
                dtmp[,"ks"]    <- eatTools::num.to.cat ( x = dtmp[,"bista"], cut.points = cutScores[["values"]], cat.values = cutScores[["labels"]])
    ### jetzt noch die shortversion der Aequivalenztabelle erzeugen
                shrt  <- do.call("rbind", by ( data = dtmp, INDICES = dtmp[,"ks"], FUN = function ( sks ) { data.frame ( score = paste(c(min(sks[,"score"]), max(sks[,"score"])), collapse=" bis "), estimate = paste(round(c(min(sks[,"est"]), max(sks[,"est"])),digits=2), collapse=" bis "), bista = paste(round(c(min(sks[,"bista"]), max(sks[,"bista"])),digits=0), collapse=" bis "), ks=unique(sks[,"ks"]), stringsAsFactors=FALSE)}))
                return(list ( complete = dtmp[,c("score", "est", "bista", "ks")], short = shrt))}


getResults <- function ( runModelObj, overwrite = FALSE, Q3 = TRUE, q3theta = c("pv", "wle", "eap"), q3MinObs = 0, q3MinType = c("singleObs", "marginalSum"), omitFit = FALSE, omitRegr = FALSE, omitWle = FALSE, omitPV = FALSE, abs.dif.bound = 0.6, sig.dif.bound = 0.3, p.value = 0.9,
              nplausible = NULL, ntheta = 2000, normal.approx = FALSE, samp.regr = FALSE, theta.model=FALSE, np.adj=8, group = NULL, beta_groups = TRUE, level = .95, n.iter = 1000, n.burnin = 500, adj_MH = .5, adj_change_MH = .05, refresh_MH = 50, accrate_bound_MH = c(.45, .55),	sample_integers=FALSE, theta_init=NULL, print_iter = 20, verbose = TRUE, calc_ic=TRUE, omitUntil=1) {
            q3MinType<- match.arg(q3MinType)
            q3theta  <- match.arg(q3theta )
            if(inherits(runModelObj, "runMultiple")) {                          ### Mehrmodellfall
                if(is.null ( attr(runModelObj, "split")[["nCores"]] ) || attr(runModelObj, "split")[["nCores"]] == 1 ) {
                   res <- lapply( runModelObj, FUN = function ( r ) {           ### erstmal single core auswertung
                          do  <- paste ( "getResults ( ", paste(names(formals(getResults)), car::recode(names(formals(getResults)), "'runModelObj'='r'"), sep =" = ", collapse = ", "), ")",sep="")
                          ret <- eval(parse(text=do))
                          return(ret)})
                   }  else  {
                          # if(!exists("detectCores"))   {library(parallel)}    ### jetzt multicore: muss dasselbe Objekt zurueckgeben!
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
                          cat(paste ( "Results of ",length(runModelObj), " analyses processed: ", sep="")); print( Sys.time() - beg)
                   }
               res <- do.call("rbind", res )
               class(res) <- c("data.frame", "multipleResults")
               rownames(res) <- NULL
               return(res)
     ### hier ist der rekursive Aufruf beendet: das folgende geschieht fuer jedes Modell einzeln, technisch auf zweierlei Weisen, je nachdem ob Conquest oder TAM gerechnet wurde
            }  else {                                                           ### Einmodellfall
               if ( is.null(runModelObj)) {return(NULL)}
               isTa  <- FALSE
               if(inherits(runModelObj, "runConquest")) {                       ### wurde mit Conquest gerechnet?
                    if ( isTRUE(Q3) ) {
                        if ( ncol ( runModelObj[["qMatrix"]]) !=2 ) {
                            cat("Q3 is only available for unidimensional models. Estimation will be skipped.\n")
                            Q3 <- FALSE
                        }
                    }
                    do    <- paste ( "res <- getConquestResults ( ", paste(names(formals(getConquestResults)), car::recode(names(formals(getConquestResults)), "'path'='runModelObj$dir'; 'analysis.name'='runModelObj$analysis.name'; 'model.name'='runModelObj$model.name'; 'qMatrix'='runModelObj$qMatrix'; 'all.Names'='runModelObj$all.Names'; 'deskRes'='runModelObj$deskRes'; 'discrim'='runModelObj$discrim'; 'daten'='runModelObj$daten'"), sep =" = ", collapse = ", "), ")",sep="")
                    eval(parse(text=do))                                        ### obere Zeile: baue Aufruf zusammen; rufe 'getConquestResults' mit seinen eigenen Argumenten auf
                    dir <- runModelObj[["dir"]]                                 ### wo Argumente neu vergeben werden, geschieht das in dem 'recode'-Befehl; so wird als 'path'-
                    name<- runModelObj[["analysis.name"]]                       ### Argument 'runModelObj$dir' uebergeben
                    allN<- runModelObj[["all.Names"]]                           ### Alternativ: es wurde mit TAM gerechnet
               }  else  {                                                       ### logisches Argument: wurde mit Tam gerechnet?
                    isTa<- TRUE                                                 ### hier wird ggf. die Anzahl der zu ziehenden PVs ueberschrieben
                    if ( isTRUE(Q3) ) {
                        if ( ncol ( attr(runModelObj, "qMatrix")) !=2 ) {
                            cat("Q3 is only available for unidimensional models. Estimation will be skipped.\n")
                            Q3 <- FALSE
                        }
                    }
                    if(!is.null(nplausible)) { attr(runModelObj, "n.plausible") <- nplausible }  else  { nplausible <- attr(runModelObj, "n.plausible") }
                    do    <- paste ( "res <- getTamResults ( ", paste(names(formals(getTamResults)), car::recode(names(formals(getTamResults)),"'pvMethod'='attr(runModelObj, \"pvMethod\")'"),  sep =" = ", collapse = ", "), ")",sep="")
                    eval(parse(text=do))
                    dir <- attr(runModelObj, "dir")                             ### untere zeilen(n): tam summary ergaenzen, wenn mit tam gerechnet wurde
                    name<- attr(runModelObj, "analysis.name")                   ### wird als Attribut in Ergebnisstruktur angehangen (nicht so superclever;
                    allN<- attr(runModelObj, "all.Names")                       ### schoener waers, man wuerde das direkt in die Ergebnisstruktur einbauen)
               }
               if(!is.null(res)) {                                              ### wenn es ein Rueckgabeobjekt gibt, wird das jetzt um einige technische Eintraege erweitert
                    stopifnot ( length(unique(res[,"model"])) == 1)             ### (das was frueher ueber Attribute gemacht wurde, aber das ist zu krass schlimm beschissen scheiss untransparent, das muss weg!!)
     ### Rueckgabeobjekt mit technischen Parametern 'anreichern', first: 'all.Names'
                    alln<- do.call("rbind", lapply(names(allN), FUN = function ( x ) {
                           if ( length( allN[[x]] ) > 0 ) {
                                res <- data.frame ( type = "tech", par = x, derived.par = allN[[x]])
                           }  else  {
                                res <- NULL
                           }
                           return(res)}))
                    res <- plyr::rbind.fill ( res, data.frame ( res[1,c("model", "source")], alln, stringsAsFactors = FALSE) )
     ### Rueckgabeobjekt mit technischen Parametern 'anreichern', second: 'dif.setting'
                    difS<- list (abs.dif.bound = abs.dif.bound, sig.dif.bound = sig.dif.bound, p.value = p.value)
                    resD<- data.frame ( res[1,c("model", "source")], type = "tech", par = "dif", derived.par = names(difS), value = unlist(difS), stringsAsFactors = FALSE)
                    res <- plyr::rbind.fill ( res, resD )
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
               return(res)
               }}

### Hilfsfunktion fuer equat1pl: konsistenzpruefungen fuer den itemparameter-Dataframe
checkItemParLists <- function (prmNorm, item, domain, testlet, value, dims = NULL) {
    ### wenn data.frame zwei spalten hat, muessen die Item- und value-Spalten nicht explizit benannt werden
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
    ### keine Missings in Spalten
           nomis<- sapply(prmNorm[,unlist(allF)], FUN = function ( i ) { length(which(is.na(i)))})
           if ( any(nomis>0)) {
                warning("Found ", length(which(nomis>0)), " column(s) in 'prmNorm' with missing values: '", paste(names(nomis[which(nomis>0)]), collapse= "', '"), "'")
           }
    ### items muessen unique sein
           tab  <- table(prmNorm[,c(allF[["item"]], allF[["domain"]]), drop=FALSE])
           if (!all(tab %in% 0:1)) {stop("Items must be unique for each domain in reference parameter frame 'prmNorm'.")}
    ### value-Spalte muss numerisch sein
           if(!inherits(prmNorm[,allF[["value"]]], "numeric")) {stop("Parameter value column in 'prmNorm' must be numeric.")}
    ### check: match domain names
           if (!is.null ( allF[["domain"]]) && !is.null(dims) ) {
                mis <- setdiff ( dims,  names(table(prmNorm[, allF[["domain"]] ])) )
                if ( length( mis ) > 0 ) { stop ( paste ( "Domain '",mis,"' is missing in 'prmNorm'.\n",sep="")) }
                uni <- by ( data = prmNorm, INDICES = prmNorm[, allF[["domain"]] ], FUN = function ( g ) {
                       if (!length(g[,allF[["item"]]]) == length(unique(g[,allF[["item"]]]))) { stop(paste ( "Item identifiers are not unique in 'prmNorm' for domain '",g[1,allF[["domain"]]],"'.\n",sep=""))}
                       }, simplify = FALSE)                                     ### check: items unique within domains?
           }
           return(allF)}

### hilfsfunktion fuer equat1pl: transformiert itemparameterliste in results-objekt, wenn nicht das Rueckgabeobjekt von getResults()
### die Fokusparameter enthaelt, sondern ein einfacher data.frame
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
           if ( length ( weg ) > 0 )  {                                         ### damit keine Spalten durch 'recode' doppelt benannt werden,
                results <- results[, -match(weg, colnames(results))]            ### muessen spalten, die sich durch die Recodierung aendern
           }                                                                    ### und zugleich schon im datensatz 'results' vergeben sind, raus
           allF <- allF[which(!sapply(allF, is.null))]
           toRec<- lapply(names(allF), FUN = function ( ff ) { paste ( "'",allF[[ff]],"'='",car::recode(ff, "'item'='item'; 'domain'='dimension'; 'value'='est'"),"'",sep="")})
           toRec<- paste(toRec, collapse = "; ")
           colnames(results) <- car::recode (colnames(results), toRec)
           return(list(results=results, dims=dims))}

### hilfsfunktion fuer equat1pl: baut leeres results objekt, falls equating nur durchgeschleift werden soll
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
           }  else  {                                                           ### jetzt fuer den Fall dass PVs bayesianische gezogen wurden und keine Itemparameter im results objekt stehen
                 resX <- results[which(!is.na(results[,"group"])),]
                 dimN <- by ( data = results, INDICES = results[,"group"], FUN = function ( prmDim ) {
                         eq <- list(B.est = c(Mean.Mean=0 , Haebara =0, Stocking.Lord=0), descriptives = c(N.Items =0, SD=NA,  Var=NA, linkerror=NA))
                         return ( list ( eq = eq, items = prmDim, method = method ) ) }, simplify = FALSE)
           }
           return(dimN)}

### hilfsfunktion fuer equat1pl: schreibt Informationen auf die Konsole
printToConsole <- function(d, nMods, it, prmDim, eq, allN, method, estimation, eqh, eqr) {
           cat(paste("\n",paste(rep("=",100),collapse=""),"\n \nModel No. ",match(d[1,"model"], names(nMods)),"\n    Model name:                ",d[1,"model"],"\n    Number of dimension(s):    ",length(unique(it[,"dimension"])),"\n    Name(s) of dimension(s):   ", paste( names(table(as.character(it[,"dimension"]))), collapse = ", "),"\n",sep=""))
           if  ( length(names(table(as.character(it[,"dimension"])))) > 1) {  cat(paste("    Name of current dimension: ",names(table(prmDim[,"dimension"]))," \n",sep=""))}
           cat(paste("    Number of linking items:   " , eq[["descriptives"]][["N.Items"]],"\n",sep=""))
           if ( !is.null(allN[["testlet"]]) ) { cat(paste( "    Number of testlets:        ",  eq[["ntl"]],"\n",sep="")) }
           cat(paste("    Linking method:            " , method,"\n",sep=""))
           if (method == "robust") { cat(paste("    Optimal trimming param.:   " , eqr[["kopt"]],"\n",sep="")) }
           if (method == "Haberman") {
               cat(paste("    Estimation method:         " , car::recode(estimation,"'OLS'='ordinary least squares'; 'BSQ'='bisquare weighted regression'; 'HUB'='regression using Huber weights'; 'MED'='median regression'; 'LTS'='trimmed least squares'; 'L1'='median polish'; 'L0'='minimizing number of interactions'"), "\n",sep=""))
               tf <- capture.output(summary(eqh))
               i1 <- grep("Used trimming factor", tf)
               i2 <- grep("Estimation information item intercepts", tf)
               i3 <- min(i1[which(i1>i2)])
               i4 <- unlist(strsplit(tf[i3], "="))
               cat(paste("    Used trimming factor:      " , round(as.numeric(eatTools::crop(i4[length(i4)])), digits = 3), "\n",sep=""))   }}

### hilfsfunktion fuer equat1pl: behandelt linking DIF
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
                          eq  <- eq1                                            ### Achtung! das muss stehen bleiben, da in obere Zeile noch Elemente aus 'eq' abgerufen werden!
    ### hier beginnt iterativer Ausschluss von Linking-DIF
                    }  else  {
                          info<- data.frame ( method = "iterativ", iter = 0 , itemExcluded = "" , DIF.excluded="", linking.constant = eq[["B.est"]][[method]], linkerror = eq[["descriptives"]][["linkerror"]] )
                          qp1 <- prmM
                          iter<- 1
                          while  ( length ( prbl ) > 0 ) {                      ### untere Zeile: finde maximalen dif-wert
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

### hilfsfunktion fuer equat1pl: wenn es keine Linking-Dif Items gibt bzw. wenn methode 'robust' ist
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

### hilfsfunktion fuer equat1pl
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

### fuer mehr als zwei Messzeitpunkte muss die Funktion rekursiv werden (ruft sich selber auf)
equat1pl<- function ( results , prmNorm , item = NULL, domain = NULL, testlet = NULL, value = NULL, excludeLinkingDif = TRUE, difBound = 1, iterativ = FALSE, method = c("Mean.Mean", "Haebara", "Stocking.Lord", "robust", "Haberman"),
           itemF = NULL, domainF = NULL, testletF = NULL, valueF = NULL, estimation=c("OLS", "BSQ", "HUB", "MED", "LTS", "L1", "L0"), b_trim=Inf, lts_prop=.5) {
           estimation <- match.arg(estimation)
           method     <- match.arg(method)
    ### Achtung! Funktion kann (was der default ist) sowohl Output aus 'runModel' weiterverarbeiten oder einfach nur zum Equaten zweiter verschiedener Itemparameterlisten genutzt werden
    ### Also muss die Funktion erstmal rauskriegen, was von beidem gemacht werden soll. Wenn 'isRunM' gleich TRUE, dann default
           isRunM<- all(c("model" , "source" , "var1" , "var2" , "type" , "indicator.group", "group", "par", "derived.par", "value") %in% names(results))
           if ( isRunM) {
                nMods <- table(results[,"model"])
                cat(paste("Found ", length(nMods), " model(s).\n   Equating is executed for each dimension in each model separately.\n",sep=""))
                dims  <- unique(unlist(by ( data = results, INDICES = results[,"model"], FUN = function ( x ) { names(table(as.character(itemFromRes(x)[,"dimension"])))})))
                if ( is.null(dims)) {
                     dims <- unique(na.omit(results[,"group"]))
                     warning(paste0("Cannot extract dimensions from 'results' object. This should only occur for bayesian plausible values imputation. Assume following dimensions: \n    '",paste(dims, collapse = "', '"),"'."))
                }
    ### jetzt beginnt der Fall, dass einfach nur zwei Itemparameterlisten equatet werden sollen. Dazu transformiert die Funktion den Input 'results' in das
    ### Format, das die Funktion 'itemFromRes()' erzeugt. Zusaetzlich wird eie weitere Spalte 'model' ergaenzt, damit die by-Funktion darueber schleifen kann
           }  else  {
                resList <- transformItemParListIntoResults (results = results, itemF = itemF, domainF = domainF, testletF = testletF, valueF = valueF)
                results <- resList[["results"]]
                dims    <- nMods <- resList[["dims"]]
           }
           if ( missing ( prmNorm) ) {                                          ### kein Equating: Rueckgabe NULL, 'results' wird durchgeschleift
                if ( isFALSE(isRunM) ) { stop("No norm parameter defined ('prmNorm' is missing).\n")}
                cat("No norm parameter defined ('prmNorm' is missing). Treat current sample as drawn from the reference population.\n")
    ### Kein equating: baue 'leeres' Rueckgabeobjekt
                items <- by ( data = results, INDICES = results[,"model"], FUN = function ( d ) {
                         dimN <- buildEmptyResultsObject(d=d, method = method, results=results)
                         return(dimN)}, simplify = FALSE)
                ret   <- list(items = items, results = results)                 ### die Klasse des Rueckgabeobjekts heisst hier "eq2tom", das ist das Equatingobjekt
                class(ret) <- c("eq2tom", class(ret))                           ### fuer 2 Messzeitpunkte (time of measurement)
                return(ret)
    ### Equating: baue 'befuelltes' Rueckgabeobjekt
           }  else {                                                            ### plausibility checks
                prmNorm<- eatTools::makeDataFrame(prmNorm)
                allN   <- checkItemParLists(prmNorm =prmNorm, item = item, domain = domain, testlet = testlet, value = value, dims=dims)
    ### Fuer jedes Modell und jede Dimension innerhalb jedes Modells findet Equating separat statt
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
    ### items muessen unique sein
                               if ( length(prmDim[, "item"]) != length(unique(prmDim[, "item"])) ) {  stop(paste("Items are not unique for model '",as.character(d[1,"model"]),"'.\n",sep="")) }
                               eq  <- equAux ( x = prmDim[ ,c("item", "est")], y = prmM[,c(allN[["item"]], allN[["value"]], allN[["testlet"]])] )
                               if ( eq[["descriptives"]][["N.Items"]] > 0) {
                                     if ( method == "robust") {
                                         prm<- merge(prmM[,c(allN[["item"]], allN[["value"]], allN[["testlet"]])], prmDim[ ,c("item", "est")], by.y="item", by.x = allN[["item"]], all=FALSE)
                                         eqr<- sirt::linking.robust(prm)        ### Haberman: der waehlt das Vorzeichen nach alphabetischer Reihenfolge in "study"
                                     }                                          ### Achtung: fuer Haberman spielt es eine Rolle, wieviele nicht-link-items zusaetzlich mit drin sind
                                     if ( method == "Haberman") {               ### also die Linkingkonstante ist eine andere, wenn man nur die gemeinsamen items im Objekt 'prm' drinhat, oder wenn alle drin sind
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
    ### Konsolenoutput (Teil 1) erstellen
                               foo <- printToConsole(d=d, nMods=nMods, it=it, prmDim=prmDim, eq=eq, allN=allN, method=method, estimation=estimation, eqh=eqh, eqr=eqr)
    ### Gibt es Items mit linking dif?
                               if ( method != "robust" && method != "Haberman" && length( prbl ) > 0 ) {
                                    eld  <- handleLinkingDif(prmDim=prmDim,prbl=prbl, eq=eq, difBound=difBound, dif=dif, method=method, excludeLinkingDif=excludeLinkingDif, iterativ=iterativ,prmM=prmM, allN=allN)
                               }  else  {                                       ### hier folgt der Abschnitt, wenn es keine Linking-Dif Items gibt bzw. wenn methode 'robust' ist
                                    eld  <- noLinkingDif(method=method, eq=eq, eqr=eqr, eqh=eqh)
                               }
    ### Konsolenoutput (Teil 2) ausgeben lassen
                               if( isFALSE(iterativ) && excludeLinkingDif && !is.null(eld[["info2"]])) {
                                   cat("\nItems with DIF:\n")
                                   print(eatTools::roundDF(eld[["info2"]][,1:2], digits = 3)); flush.console()
                               }
                               cat("\n")
                               print(eatTools::roundDF(eld[["info"]])); flush.console()
                               cat("\n")
    ### Output fuer die Methoden robust und Haberman in die einheitliche Struktur der anderen Methoden bringen
                               if ( method %in% c("robust", "Haberman")) {
                                    eld <- createOutput(method=method, eqr=eqr, prm=prm, eqh=eqh, info=eld[["info"]])
                               }
                               ret <- list ( eq = eld[["eq"]], items = prmDim, info = eld[["info"]], method = method )
                               return ( ret ) }, simplify =FALSE)
                       return(dimN) }, simplify = FALSE)
              ret  <- list ( items = items, results = results)                  ### die Klasse des Rueckgabeobjekts heisst hier "eq2tom", das ist das Equatingobjekt
              class(ret) <- c("eq2tom", class(ret))                             ### fuer 2 Messzeitpunkte (time of measurement)
              return(ret)                                                       ### "results"-Objekt wird durchgeschleift
              }  }

### Hilfsfunktion fuer equat1pl
equAux  <- function ( x, y ) {
           eq  <- sirt::equating.rasch(x = x, y = y[,1:2])                      ### kein Jackknife
           if ( ncol(y)==3) {                                                   ### jackknife
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

### Hilfsfuntion fuer "transformToBista", solange eatrep im Uebergangsstadium ist
adaptEatRepVersion <- function ( x ) {
     if ( inherits(x, "data.frame"))  {
           return ( x )
     }  else  {
           x <- x[[1]][[1]]
           stopifnot ( inherits(x, "data.frame") )
           return(x)
     } }

### Hilfsfuntion fuer "transformToBista"
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

### Hilfsfuntion fuer "transformToBista"
createItemVeraObj <- function(itempars, roman){
       pCols      <- colnames(itempars)[grep("^itemP", colnames(itempars))]
       allCols    <- na.omit(match ( c("dimension","item", pCols, "itemDiscrim", "estTransf", "infit", "estTransfBista", "traitLevel"), colnames(itempars)))
       itemVera   <- itempars[,allCols]
       colnames(itemVera) <- car::recode ( colnames(itemVera), "'dimension'='domain'; 'item'='iqbitem_id'; 'itemDiscrim'='trennschaerfe'; 'estTransf'='logit'; 'estTransfBista'='bista'; 'traitLevel'='kstufe'")
       colnames(itemVera)[match(pCols, colnames(itemVera))] <- paste0("lh", eatTools::removePattern ( string = pCols, pattern = "^itemP"))
       if ( roman == TRUE ) {                                                   ### sollen roemische Zahlen fuer Kompetenzstufen verwendet werden?
            if (!all(itemVera[,"kstufe"] %in% c("1a", "1b", 1:5))) {stop(paste("Competence levels do not match allowed values. '1a', '1b', '1', '2', '3', '4', '5' is allowed. '",paste(names(table(itemVera[,"kstufe"])), collapse = "', '"),"' was found.\n",sep=""))}
            itemVera[,"kstufe"] <- car::recode (itemVera[,"kstufe"], "'1a'='Ia'; '1b'='Ib'; '1'='I'; '2'='II'; '3'='III'; '4'='IV'; '5'='V'")
       }
    ### jetzt den scheiss reshapen fuer vera-3 mathe, wenn es separate werte fuer global- und domaenenspezifische modelle gibt
       if ( length ( unique ( itemVera[,"iqbitem_id"])) != length ( itemVera[,"iqbitem_id"]) ) {
            cat("Found duplicated entries in 'item-ID' column. This should only occur for subject 'math' in grade 3.\n")
    ### vorher rauskriegen, ob jedes item nur zu einer domaene und einem globalmodell gehoert, dann 'domain' umbenennen, um NAs im Ergebnis zu vermeiden
            tab  <- table(itemVera[,c("domain", "iqbitem_id")])
            if ( !"GL" %in% rownames(tab)) {
                 cat("Cannot find 'global' entry in the 'domain' column. Cancel reshaping.\n")
            }  else  {
                 if ( !sum(tab[which(rownames(tab) == "GL"),]) == ncol(tab)) {  ### hat jedes Item einen Wert auf 'global'?
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
       return(itemVera) }

transformToBista <- function ( equatingList, refPop, cuts, weights = NULL, defaultM = 500, defaultSD = 100, roman = FALSE, vera = TRUE, idVarName = NULL, years = NULL ) {
    ### wenn equatet wurde, sollte auch 'refPop' definiert sein (es sei denn, es wurde verankert skaliert)
    ### wenn 'refPop' fehlt, wird es fuer alle gegebenen Dimensionen anhand der Gesamtstichprobe berechnet
       mr  <- FALSE                                                             ### default: 'refPop' fehlt nicht. Wenn doch, wird es aus Daten generiert und spaeter
       if(missing(refPop)) {                                                    ### (nachdem ggf. transformiert wurde!) auf Konsole angezeigt
          mr  <- TRUE
          cat("'refPop' was not defined. Treat current sample as drawn from the reference population.\n")
          flush.console()
       }
       if( missing(cuts)) { cutsMis <- TRUE }  else  { cutsMis <- FALSE }
       nam1<- names(equatingList[["items"]])                                    ### hier stehen in der regel die beiden Modellnamen, also quasi names(table(equatingList[["results"]][,"model"]))
    ### Fuer jedes Modell und jede Dimension findet Transformation separat statt
    ### Schritt 1: 'refPop' bestimmen, wenn es sie noch nicht gibt ... 'refPop' ist ggf. verschieden ueber Itemgruppen, aber immer gleich ueber Personengruppen!
    ### fuer die Bestimmung von 'refPop' wird nur ueber Itemgruppen gesplittet!
       it     <- itemFromRes(equatingList[["results"]])
       dims   <- unique(it[,"dimension"])
       if ( is.null(dims)) {
            dims <- unique(na.omit(equatingList[["results"]][,"group"]))
            warning(paste0("Cannot extract dimensions from 'results' object. This should only occur for bayesian plausible values imputation. Assume following dimensions: \n    '",paste(dims, collapse = "', '"),"'."))
            if(vera==TRUE) {
               warning("'vera' must be FALSE for bayesian plausible values imputation. Set 'vera' to FALSE.")
               vera <- FALSE
            }
       }
    ### Hotfix: 'id'-Variable identifizieren ... um Kompatibilitaet mit aelteren Paketversionen zu wahren, kann man hier die ID auch von Hand eingeben
       id     <- unique(equatingList[["results"]][intersect(which(equatingList[["results"]][,"type"] == "tech"), which(equatingList[["results"]][,"par"] == "ID")),"derived.par"])
       if(length(id)!=1) {
           id   <- getIdVarName(id=NULL, idVarName, verbose=TRUE)
       }
       refList<- lapply ( dims, FUN = function (dimname) {
                 rex  <- pvFromRes(equatingList[["results"]][unique(c(which(equatingList[["results"]][,"group"] == dimname),which(equatingList[["results"]][,"type"] == "tech"))), ], toWideFormat = FALSE, idVarName=idVarName, verbose=FALSE)
                 if (is.null(rex)) {return(NULL)}                               ### NULL wird zurueckgegeben, wenn keine PVs in der Ergebnisstrauktur vorhanden waren
                 if ( is.null(weights) ) {
                      txt <- capture.output ( msd <- eatRep::repMean ( datL = rex, ID = id, imp = "imp", dependent = "value", na.rm = TRUE))
    ### Achtung!! ggf. anpassen fuer neue eatRep-Version
                      msd <- adaptEatRepVersion(msd)                            ### in transformToBista() muessen die messages von mergeAttr nicht uber capture.output abgefangen werden
                 }  else  {                                                     ### wie in anker(), da transformToBista() immer nur single core aufgerufen werden kann
                    rex <- eatTools::mergeAttr ( rex, weights , by.x = id, by.y = colnames(weights)[1], all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "cases", xName = paste0("plausible values for dimension ",dimname), yName = "weights", verbose = c("match", "dataframe"))
                    mis <- which(is.na(rex[,colnames(weights)[2]]))
                    if ( length(mis) > 0 ) {                                    ### missings in the weights frame are not allowed
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
    ### wenn 'refPop' nicht definiert wurde, wird es hier mit Werten gesetzt, die direkt aus der Stichprobe (= Normpopulation) berechnet wurden
       ref    <- do.call("rbind", lapply(refList, FUN = function ( u ) { u[["rp"]] }))
       if ( isTRUE(mr) ) {
          refPop <- ref
       }   else  {
    ### wenn 'refPop' NUR FUER EINE DIMENSION nicht definiert wurde, werden hier die nicht definierten ('NA') Werte durch Werte aus der Stichprobe (= Normpopulation fuer genau diese Dimension) ersetzt
          mis <- which(is.na(refPop))
          if ( length(mis) >0) {
               stopifnot ( nrow(ref ) == nrow(refPop))
               mat <- merge( 1:nrow(refPop), 1:ncol(refPop), by = NULL)         ### rauskriegen, in welchen Zeilen und Spalten die werte fehlen
               refPop[unique(mat[mis,"x"]), unique(mat[mis,"y"])] <- ref[unique(mat[mis,"x"]), unique(mat[mis,"y"])]
          }
       }
       if(ncol ( refPop ) == 3) {
          cat ( paste("The 'refPop' data.frame does not include information about reference population mean/SD on Bista metric. Values will be defaulted to ",defaultM,"/",defaultSD,".\n",sep=""))
          refPop[,4] <- defaultM; refPop[,5] <- defaultSD
       }  else  {
          if ( ncol ( refPop) != 5 ) { stop ( "Invalid 'refPop'.\n") }
       }
    ### fuer die Transformation selbst wird nach Modellen und Dimensionen getrennt
       modN<- lapply(nam1, FUN = function ( mod ) {                             ### aeussere Schleife: geht ueber modelle
              nam2 <- names(equatingList[["items"]][[mod]])
              dimN <- lapply(nam2, FUN = function ( dims ) {                    ### innere Schleife: geht ueber Dimensionen (innerhalb von modellen)
    ### check: sind Personen innerhalb jeder Dimension (und jeder Imputation) unique? ... 'redMD' ist ein reduziertes Results-Objekt: nur die interessierende Dimension des interessierenden Modells + saemtliche "tech"-Variablen
                      resMD<- equatingList[["results"]][unique(c(intersect(which(equatingList[["results"]][,"model"] == mod), which(equatingList[["results"]][,"group"] == dims)),  which(equatingList[["results"]][,"type"] == "tech"))),]
                      rex  <- pvFromRes(resMD, toWideFormat = TRUE, idVarName = idVarName, verbose=FALSE)
                      if (!is.null(rex)) {
                          if ( length ( rex[,id]) != unique(length ( rex[,id])) ) {
                               stop(paste( "Model '",mod,"', Dimension '",dims,"': cases according to '", id,"' variable are not unique.\n",sep=""))
                          }
                      }
    ### check: keine verankerten parameter?
                      offSet  <- grep("offset", as.character(resMD[,"par"]))
                      if(length(offSet)>0) {  resMD[,"par"] <- car::recode ( resMD[,"par"], "'offset'='est'") }
                      itFrame <- itemFromRes(resMD)
                      if ( !is.null(itFrame) && !itFrame[1,"dimension"] %in% refPop[,1] ) {
                            cat(paste("Cannot found dimension '",itFrame[1,"dimension"],"' in the first column of the 'refPop' argument. Skip transformation ... \n",sep=""))
                            return ( list ( itempars = NULL, personpars = NULL, rp = NULL))
                      }  else  {
                            if ( is.null ( itFrame ) ) {
                                cat(paste0("Model '",mod,"', dimension '",dims,"': No item parameters found. This should only occur for bayesian plausible values imputation. Transformation of item parameters will be skipped.\n"))
                            }  else  {
    ### wenn cuts vom user definiert, wird hier auf plausibilitaet geprueft (kuenftig ggf. auslagern in separate funktion)
                                if ( isFALSE(cutsMis) ) {
                                     if ( !itFrame[1,"dimension"] %in% names(cuts) ) { stop(paste("Cannot found dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
                                     mat1<- match( itFrame[1,"dimension"], names(cuts))
                                     if ( !"values" %in% names(cuts[[mat1]]) ) { stop(paste("'cuts' must be a named list. Cannot found 'values' element for dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
                                     if ( length(cuts[[mat1]])>1) {
                                         if ( !"labels" %in% names(cuts[[mat1]]) ) { stop(paste("'cuts' must be a named list. Cannot found 'labels' element for dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
                                     }
                                }
                                mat <- match( itFrame[1,"dimension"], refPop[,1])
    ### 1. Transformation fuer Itemparameter
                                itFrame[,"estTransf"] <- itFrame[,"est"] + equatingList[["items"]][[mod]][[dims]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dims]][["method"]] ]]
    ### Achtung, heikel: wenn equatet wurde, aber der Datensatz aus der Normpopulation kommt, werden hier die empirischen Mittelwerte,
    ### die oben (mit oder ohne Gewichte) berechnet wurden, nochmal transformiert ... sollte praktisch nie der Fall sein.
                                if ( isTRUE(mr) ) {
                                     if ( equatingList[["items"]][[mod]][[dims]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dims]][["method"]] ]] != 0) {
                                          cat("W A R N I N G: Preceding Equating without 'refPop' definition. Sure you want to use current sample as drawn from the reference population?\n")
                                          refPop[mat,2] <- refPop[mat,2]+ equatingList[["items"]][[mod]][[dims]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dims]][["method"]] ]]
                                     }
                                }
                                itFrame[,"estTransf625"]   <- itFrame[,"estTransf"] + log(0.625/(1-0.625))
                                itFrame[,"estTransfBista"] <- (itFrame[,"estTransf625"] - refPop[mat,2]) / refPop[mat,3] * refPop[mat,5] + refPop[mat,4]
                                if ( isFALSE(cutsMis) ) {
    ### Achtung: dieser Umweg ist notwendig, weil 'num.to.cat' Attribute ausgibt die unten wieder gebraucht werden!
                                     traitLevel            <- eatTools::num.to.cat(x = itFrame[,"estTransfBista"], cut.points = cuts[[mat1]][["values"]], cat.values = cuts[[mat1]][["labels"]])
                                     itFrame[,"traitLevel"]<- traitLevel
                                }
    ### Achtung!! Linkingfehler sollte eigentlich nur ausgegeben werden, wenn das vorherige equating NICHT durchgeschleift wurde!
                                itFrame[,"linkingConstant"]<- equatingList[["items"]][[mod]][[dims]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dims]][["method"]] ]]
                                itFrame[,"linkingMethod"]  <- equatingList[["items"]][[mod]][[dims]][["method"]]
                                itFrame[,"nLinkitems"]     <- equatingList[["items"]][[mod]][[dims]][["eq"]][["descriptives"]][["N.Items"]]
                                itFrame[,"linkingError"]   <- equatingList[["items"]][[mod]][[dims]][["eq"]][["descriptives"]][["linkerror"]]
    ### Transformation des Linkingfehlers entsprechend der Rechenregeln fuer Varianzen. ist geprueft, dass dasselbe rauskommt, wie wenn man Parameter transformiert und dann Linkingfehler bestimmt
                                itFrame[,"linkingErrorTransfBista"] <- ( (itFrame[,"linkingError"]^2) * (refPop[mat,5]^2) / (refPop[mat,3]^2) )^0.5
    ### Deltamethode, wie in eatTrend (Funktion 'seKompstuf'). Dazu wird MW und SD der Fokuspopulation benoetigt! (wurde oben als 'msd' berechnet)
    ### das ganze findet nur statt, wenn sowohl cut scores bereits definiert sind und wenn equatet wurde (denn nur dann gibt es einen Linkingfehler, den man transformieren kann)
                            }
                            pv  <- pvFromRes(resMD, toWideFormat = FALSE, idVarName=idVarName, verbose=FALSE)
                            equ <- equatingList[["items"]][[mod]][[dims]][["eq"]][["B.est"]][[ equatingList[["items"]][[mod]][[dims]][["method"]] ]]
    ### Hotfix fuer bayesianisch
                            if (!exists("mat")) { mat <- match(dims,  refPop[,1]) }
                            pv[,"valueTransfBista"] <- (pv[,"value"] + equ - refPop[mat,2]) / refPop[mat,3] * refPop[mat,5] + refPop[mat,4]
    ### Dazu muss zuerst Mittelwert und SD der Fokuspopulation bestimmt werden.
                            if (!is.null(pv)) {
                                if ( is.null(weights) ) {
                                     txt <- capture.output ( msdF <- eatRep::repMean ( datL = pv, ID = id, imp = "imp", dependent = "valueTransfBista", na.rm = TRUE))
                                     msdF<- adaptEatRepVersion(msdF)
                                }  else  {
                                     pvF <- eatTools::mergeAttr ( pv, weights , by.x = id, by.y = colnames(weights)[1], all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "cases", xName = paste0("plausible values for dimension ",dims), yName = "weights", verbose = c("match", "dataframe"))
                                     mis <- which(is.na(pvF[,colnames(weights)[2]]))
                                     if ( length(mis) > 0 ) {                       ### missings in the weights frame are not allowed
                                          cat(paste ( "Found ",length(mis)," missing values in the 'weights' frame.\n    Cases with missing values on weighting variable will be ignored for transformation.\n",sep=""))
                                          pvF <- pvF[-mis,]
                                     }
                                     txt <- capture.output ( msdF <- eatRep::repMean ( datL = pvF, ID = id, imp = "imp", wgt = colnames(weights)[2], dependent = "valueTransfBista", na.rm = TRUE) )
                                     msdF<- adaptEatRepVersion(msdF)
                                }
                                msdFok <- c(msdF[intersect(which(msdF[,"parameter"] == "mean"), which(msdF[,"coefficient"] == "est")),"value"], msdF[intersect(which(msdF[,"parameter"] == "sd"), which(msdF[,"coefficient"] == "est")),"value"])
                            }  else  {
                                cat("Results object does not contain any plausible values. Skip transformation of linking error for competence levels.\n")
                            }
                            if ( !is.null(pv) && !is.null ( itFrame )) {        ### Cuts mit Schwelle nach unten und nach oben offen
                                if ( cutsMis == FALSE & !is.null ( equatingList[["items"]] )) {
                                     cts <- c( -10^6, cuts[[mat1]][["values"]], 10^6)
                                     le  <- do.call("rbind", lapply ( (length(cts)-1):1 , FUN = function ( l ) {
                                            kmp<- c(cts[l], cts[l+1])           ### Linkingfehler fuer einzelnen Kompetenzintervalle; absteigend wie bei karoline
                                            a1 <- sum ( dnorm ( ( kmp - refPop[mat,4]) / refPop[mat,5] ) * c(-1,1) / refPop[mat,5] )
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
                                     itFrame <- eatTools::mergeAttr ( itFrame, le, by = "traitLevel", sort = FALSE, all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "trait levels", xName = "item parameter list", yName = "linking error list", verbose = c("match"))
                                     itFrame <- itFrame[,c(ori, "linkingErrorTraitLevel")]
                                }
                                itFrame[,"refMean"]        <- refPop[mat,2]
                                itFrame["refSD"]           <- refPop[mat,3]
                                itFrame[,"refTransfMean"]  <- refPop[mat,4]
                                itFrame[,"refTransfSD"]    <- refPop[mat,5]
                            }
    ### 2. Transformation der Personenparameter: kann auch dann stattfinden, wenn PVs bayesianisch gezogen wurden
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
    ### 'refPop' Informationstabelle bauen
                            rp  <- refPop[mat,]
                            colnames(rp) <- c("domain", "refMean", "refSD", "bistaMean", "bistaSD")
                            rp  <- cbind ( model = mod, rp, focusMean = msdFok[1], focusSD = msdFok[2])
                            return(list ( itempars = itFrame, personpars = pv, rp = rp))
                      }  })                                                     ### untere Zeile: hier muss man rbind.fill nehmen, das gibt sonst im LV2021 einen Fehler, wenn bei der PV-Ziehung manche Items verankert sind, andere nicht
              itempars<- do.call(plyr::rbind.fill, lapply ( dimN, FUN = function ( x ) { x[["itempars"]]}))
              perspar <- do.call("rbind", lapply ( dimN, FUN = function ( x ) { x[["personpars"]]}))
              rp      <- do.call("rbind", lapply ( dimN, FUN = function ( x ) { x[["rp"]]}))
              return( list ( itempars = itempars, personpars = perspar, rp=rp)) } )
       personpars <- do.call("rbind", lapply ( modN, FUN = function ( x ) { x[["personpars"]]}))
       itempars   <- do.call(plyr::rbind.fill, lapply ( modN, FUN = function ( x ) { x[["itempars"]]}))
       rp         <- do.call("rbind", lapply ( modN, FUN = function ( x ) { x[["rp"]]}))
    ### jetzt noch die Itemparameterliste fuer die Vergleichsarbeiten reduzieren und aufbereiten
       if ( vera == FALSE ) {
           itemVera <- NULL
       }  else  {
           itemVera <- createItemVeraObj(itempars=itempars, roman=roman)
       }
    ### optional: separate linkingfehlerobjekte erzeugen
       if (!is.null(years)) {
           stopifnot(length(years) == 2 && length(unique(years)) == 2)
           leo  <- createLinkingErrorObject(itempars=itempars, years=years)
       }  else  {
           leo  <- NULL
       }
       context    <- equatingList[["results"]][which(equatingList[["results"]][,"type"]=="tech"),]
       ret        <- list ( itempars = itempars, personpars = personpars, refPop = refPop, means = rp, all.Names = context, itemparsVera = itemVera, linkingErrors = leo)
       class(ret) <- c("list", "transfBista")
       return( ret ) }


runModel <- function(defineModelObj, show.output.on.console = FALSE, show.dos.console = TRUE, wait = TRUE) {
            if (inherits(defineModelObj, "defineMultiple") ) {                  ### erstmal fuer den Multimodellfall: nur dafuer wird single core und multicore unterschieden
                if(is.null ( attr(defineModelObj, "split")[["nCores"]] ) || attr(defineModelObj, "split")[["nCores"]] == 1 ) {
                   res <- lapply(defineModelObj, FUN = function ( r ) {         ### erstmal: single core
                          ret <- runModel ( defineModelObj = r, show.output.on.console = show.output.on.console, show.dos.console = show.dos.console, wait = wait)
                          return(ret)})
                }  else  {                                                      ### multicore
                   # if(!exists("detectCores"))   {library(parallel)}
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
                   cat(paste ( length(defineModelObj), " analyses finished: ", sep="")); print( Sys.time() - beg)
                }
                class(res) <- c("runMultiple", "list")
                attr(res, "split") <- attr(defineModelObj, "split")
                return(res)
            } else {                                                            ### ab hier fuer den single model Fall
                if(inherits(defineModelObj, "defineConquest")) {                 ### hier fuer conquest
                   oldPfad <- getwd()
                   setwd(defineModelObj$dir)
                   suppressWarnings(system(paste(defineModelObj$conquest.folder," ",defineModelObj$input,sep=""),invisible=!show.dos.console,show.output.on.console=show.output.on.console, wait=wait) )
                   if(wait == FALSE) { Sys.sleep(0.2) }
                   setwd(oldPfad)                                               ### untere Zeile: Rueckgabeobjekt definieren: Conquest
                   class(defineModelObj) <- c("runConquest", "list")
                   return ( defineModelObj )
                }
                if(inherits(defineModelObj, "defineTam")) {
                   if ( show.output.on.console == TRUE ) { control$progress <- TRUE }
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
                   attr(mod, "qMatrix")      <- defineModelObj[["qMatrix"]]     ### hier werden fuer 'tam' zusaetzliche Objekte als Attribute an das Rueckgabeobjekt angehangen
                   attr(mod, "n.plausible")  <- defineModelObj[["n.plausible"]] ### Grund: Rueckgabeobjekt soll weitgehend beibehalten werden, damit alle 'tam'-Funktionen, die darauf aufsetzen, lauffaehig sind
                   attr(mod, "dir")          <- defineModelObj[["dir"]]
                   attr(mod, "analysis.name")<- defineModelObj[["analysis.name"]]
                   attr(mod, "all.Names")    <- defineModelObj[["all.Names"]]
                   attr(mod, "deskRes")      <- defineModelObj[["deskRes"]]
                   attr(mod, "discrim")      <- defineModelObj[["discrim"]]
                   attr(mod, "irtmodel")     <- defineModelObj[["irtmodel"]]
                   attr(mod, "pvMethod")     <- defineModelObj[["pvMethod"]]
                   attr(mod, "Y")            <- Y
                   return(mod)  }  }   }

### Hilfsfunktiuon fuer 'defineModel': Jetzt wird die aufbereitete Liste aus 'splitModels' abgearbeitet. ACHTUNG: Argumente in 'splittedModels' ueberschreiben default- und vom Nutzer gesetzte Argumente in 'defineModel'!
### Der Funktionsaufruf von 'doAufb' variiert je nach single- oder multicore handling. Die Funktion muss wissen, ob sie mit single- oder multicore aufgerufen wird,
### weil davon abhaengt, ob messages einfach so ausgegeben, oder ueber capture() gespeichert und erst nach dem multicore auf die konsole geprintet werden sollen
doAufb <- function ( m, matchCall, anf, verbose, dir, multicore ) {
          matchL <- match(m, unlist(lapply(matchCall[["splittedModels"]][["models.splitted"]], FUN = function ( l ) { l[["model.no"]] } )))
          mess1  <- NULL                                                        ### Nachrichtenobjekt initialisieren
          if(!is.null(matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["qMatrix"]])) {
     ### check: wenn superSplitter BERUHEND AUF ITEM GROUPING genutzt wird, wird 'items'-Argument von 'defineModel' ignoriert; wenn das also im 'matchCall' NICHT NULL ist, wird es ignoriert
             if ( !is.null(matchCall[["items"]]) )  {                           ### Warnung nur beim ersten Schleifendurchlauf anzeigen!
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
     ### Personen im Datensatz selektieren: Achtung: wenn keine Personen in "person.grouping", nimm alle!
          if(!is.null(matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["person.grouping"]])) {
             persMis<- setdiff ( matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["person.grouping"]][,1], matchCall[["dat"]][,matchCall[["id"]]])
             if( length ( persMis ) > 0) {
                 mess1 <- c(mess1, paste0( "Warning: ",length(persMis) ," from ",nrow(matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["person.grouping"]])," persons not found in data.\n"))
             }
             persons<- intersect ( matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["person.grouping"]][,1], matchCall[["dat"]][,matchCall[["id"]]])
             datSel <- matchCall[["dat"]][match(persons, matchCall[["dat"]][,matchCall[["id"]]]),]
          }  else  { datSel <- matchCall[["dat"]] }
     ### Unterverzeichnisse definieren
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
     ### Aufruf von 'defineModel' generieren, Teil 1. Achtung: wenn der Nutzer eigenhaendig neue Argumente in <models>[["models"]] einfuegt, muessen die hier in <models>[["models.splitted"]] uebernommen werden!
          overwr1<- list( dat=datSel, items = itemSel, qMatrix = qMatrix, analysis.name = nameI, dir = dirI)
          overwrF<- setdiff ( colnames(matchCall[["splittedModels"]][["models"]]), c("model.no", "model.name", "model.subpath", "dim", "Ndim", "group", "Ngroup"))
          if(length(overwrF)>0) {                                               ### wenn der Nutzer zusaetzliche Spalten in <models>$models spezifiziert, muessen
             notAllow <- setdiff ( overwrF, names(formals(defineModel)))        ### die Spaltennamen zu Argumenten von 'defineModel' passen, sonst werden die ignoriert
             if ( length ( notAllow ) > 0 ) {
                  if ( m == anf ) {                                             ### folgende Warnung soll nur einmal erscheinen, obwohl es fuer jedes Modell geschieht (Konsole nicht mit Meldungen zumuellen)
                       mess1 <- c(mess1, paste("Column(s) '",paste(notAllow, collapse = "', '"),"' of 'splittedModels' definition frame do not match arguments of 'defineModel()'. Columns will be ignored.\n", sep=""))
                  }
                  overwrF <- setdiff (overwrF, notAllow)
             }                                                                  ### wenn der Nutzer zusaetzliche Spalten in <models>$models spezifiziert, duerfen bestimmte
             notAllow2<- intersect ( overwrF, names(overwr1))                   ### Argumente von 'defineModel' NICHT benutzt werden, die werden dann auch ignoriert
             if ( length ( notAllow2 ) > 0 ) {
                  if ( m == anf ) {                                             ### folgende Warnung soll nur einmal erscheinen, obwohl es fuer jedes Modell geschieht (Konsole nicht mit Meldungen zumuellen)
                       mess1 <- c(mess1, paste("Column(s) '",paste(notAllow2, collapse = "', '"),"' of 'splittedModels' definition frame are not allowed to be modified by user. Columns will be ignored.\n", sep=""))
                  }
                  overwrF <- setdiff (overwrF, notAllow2)
             }
             notAllow3<- intersect ( overwrF, names(matchCall))                 ### wenn der Nutzer bspw. im Splitter die nodes modellspezifisch setzt, duerfen die im Aufruf von 'defineModel' nicht nochmal gesetzt werden ...
             if ( length ( notAllow3 ) > 0 ) {                                  ### gibt nur eine Warnung, das Ignorieren geschieht automatisch
                  if ( m == anf ) {
                       mess1 <- c(mess1, paste("Column(s) '",paste(notAllow3, collapse = "', '"),"' were defined twice, in <models>$models and 'defineModel'. The latter one will be ignored.\n", sep=""))
                  }
             }                                                                  ### wenn nach den ganzen checks immer noch zusaetzliche Argumente uebrig sind, werden die jetzt
             if ( length ( overwrF ) > 0 ) {                                    ### in 'overwr1' ergaenzt ... und in 'splittedModels' fuer die 'sprechenden Ausgaben'
                  for ( hh in overwrF ) {
                        overwr1[[hh]] <- matchCall[["splittedModels"]][["models"]][which(matchCall[["splittedModels"]][["models"]][,"model.no"] == m),hh]
                        matchCall[["splittedModels"]][["models.splitted"]][[matchL]][[hh]] <- matchCall[["splittedModels"]][["models"]][which(matchCall[["splittedModels"]][["models"]][,"model.no"] == m),hh]
                  }
             }
          }  else  { overwrF <-  NULL }
     ### Aufruf von 'defineModel' generieren, Teil 2. Fuege Argumente aus 'matchCall' in 'overwr1' hinzu ... aber nur, die es nicht schon gibt
          zusatz <- setdiff ( setdiff ( names(matchCall), "splittedModels"), names( overwr1))
          if ( length ( zusatz ) > 0 ) { overwr1 <- c(overwr1, matchCall[zusatz]) }
     ### sprechende Ausgaben, wenn verbose == TRUE
          if(!is.null(matchCall[["items"]])) {allVars<- list(variablen=matchCall[["items"]])}
          if(exists("itemSel"))              {allVars<- list(variablen=itemSel)}### Hotfix: anzahl der Items bestimmen
          allNams<- lapply(allVars, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = matchCall[["dat"]], variable=ii)})
          overwr3<- data.frame ( arg = c("Model name", "Number of items", "Number of persons", "Number of dimensions"), eval = as.character(c(matchCall[["splittedModels"]][["models.splitted"]][[matchL]][["model.name"]],length(allNams[["variablen"]]), nrow(datSel) , nDim)), stringsAsFactors = FALSE)
           if ( length ( overwrF) > 0 )  {
                zusatz <- lapply ( overwrF, FUN = function ( y ) { matchCall[["splittedModels"]][["models"]][which(matchCall[["splittedModels"]][["models"]][,"model.no"] == m),y] })
                zusatz <- data.frame ( arg = overwrF, eval = as.character ( zusatz ) )
                overwr3<- rbind ( overwr3, zusatz)
           }
           overwr3[,"leerz"] <- max (nchar(overwr3[,"arg"])) - nchar(overwr3[,"arg"]) + 1
           txt    <- apply(overwr3, MARGIN = 1, FUN = function ( j ) { paste("\n    ", j[["arg"]], ":", paste(rep(" ", times = j[["leerz"]]), sep="", collapse=""), j[["eval"]], sep="")})
           nDots  <- max(nchar(overwr3[,"arg"])) + max(nchar(overwr3[,"eval"])) + 6
           if(verbose == TRUE ) {
              cat(paste("\n\n",paste(rep("=",times = nDots), sep="", collapse=""),"\nModel No. ",m, paste(txt,sep="", collapse=""), "\n",paste(rep("=",times = nDots), sep="", collapse=""),"\n\n", sep=""))
              if(!is.null(mess1)) { cat(mess1)}
           }
     ### Achtung! Rueckgabe haengt davon ab, ob multicore Handling stattfinden soll! zuerst single core
          if(is.null ( matchCall[["splittedModels"]][["nCores"]] ) | matchCall[["splittedModels"]][["nCores"]] == 1 ) {
             ret    <- do.call("defineModel", args = overwr1)                   ### single core handling: die verschiedenen Modelle werden
          }  else  {                                                            ### bereits jetzt an "defineModel" zurueckgegeben und seriell verarbeitet
             attr(overwr1[["dat"]], "multicore") <- TRUE                        ### hier sollte die Eigenschaft 'multicore' als Attribut des Datensatzes durchgeschleift werden
             ret    <- overwr1                                                  ### multicore: die verschiedenen Modelle werden noch nicht weiter verarbeitet,
          }                                                                     ### es wird lediglich der Modellaufruf generiert, der dann spaeter an die einzelnen
          return(ret) }                                                         ### cores weitergegeben wird

### dat               ... Datensatz als R-Dataframe
### items             ... wo stehen Items im datensatz, z.B. 5:120 oder -c(1:5)
### id                ... wo steht ID-variable, entweder Spaltennummer oder Variablenname als String
### qMatrix           ... falls nicht definiert: eindimensional, ansonsten muss hier qMatrix als R-Dataframe uebergeben werden;
###                       Dabei enthaelt die erste Spalte stets den Item-Bezeichner, die folgenden die Zugehoerigkeit zu Dimensionen
### DIF.var           ... eine DIF-Variable, entweder als Spaltennummer oder Variablenname als String
### HG.var            ... ein oder mehrere Hintergrundvariablen, entweder als Spaltennummern oder Variablennamen, z.B. c(4,6), c("alter","geschlecht")
### weight.var        ... Gewichtungsvariable, entweder als Vektor/Skalar mit Spaltennummern oder Variablennamen
### anchor            ... optional: data.frame mit Ankerparametern, 1. Spalte muss Itembezeichner sein, zweite Spalte Parameter
### dir               ... wo sollen Conquest-Inputdateien hingeschrieben werden?
### analysis.name     ... Dateiname fuer Conquest-Input (nur Praefix, Suffixe werden automatisch vergeben)
### conquest.folder   ... optional: vollstaendiger Conquest-Pfad, z.B. "N:/iqb/console.exe" Wenn spezifiziert, wird Batchdatei ausgegeben
### set.constraints   ... "none" , "cases" (default) , "items" ; bei Anchor wird automatisch auf "none" gesetzt
### boundary          ... Es wird eine Nachricht ausgegeben, wenn Personen im Datensatz sind, die weniger als [boundary] Items bearbeitet haben
### remove.boundary   ... Logical: Sollen Personen, die weniger als [boundary] Items beantwortet haben, aus den Daten entfernt werden?
### nodes             ... (optional): positive ganzzahlige zahl (integer)
### f.nodes  			    ... positive ganzahlige Zahl (integer); Conquest-Handbuch S. 225
### p.nodes    		    ... positive ganzahlige Zahl (integer); Conquest-Handbuch S. 225
### method		  	    ... optional: "gauss" (default), "quadrature", "montecarlo"; Conquest-Handbuch S. 225
### std.err           ... optional: "full", "quick" (default), "none"; Conquest-Handbuch S. 167ff
### n.interations	    ... positive ganzahlige Zahl (integer); maximale Anzahl von Iterationen (geht fuer tam und conquest)
### converge		      ... Gleitkommazahl; Conquest-Handbuch S. 225
### distribution      ... "normal" (default), "discrete"; Conquest-Handbuch S. 167ff
### equivalence.table ... Gibt ggf. Tabelle mit Umrechnungen Rohwert-Normwert; moegliche Werte sind "wle" (default); "mle" oder NULL (keine Tabelle wird ausgegeben)
###                       (Conquest-handbuch, S.166)
### use.letters       ... optional: sollen Daten in Buchstaben umcodiert werden und fuer Conquest-Analyse Buchstaben statt Ziffern benutzt werden? Ist moeglicherweise bei Partial Credit mit mehr als zehn Kategorien
###                       bedeutsam, da Conquest Probleme mit zweistelligen Codes hat. Buchstaben erlauben also bis maximal 26 Kategorien und bleiben einstellig.
### model.statement   ... model statement, wie es in Conquest erscheinen soll
### suppress.logfile  ... logical: should creation of logfile in Conquest be suppressed? Maybe important in Simulations studies (logfile is asked to overwrite in Conquest)
### check.for.linking ... logical: check whether all items are connected by design
### software          ... specifies software for analysis. If software == "lme4", only the long-format data.frame is returned. if software == "tam", dataset is prepared and tam is evaluated
### est.slopegroups   ... Matrix; erste Spalte: Item-ID, zweite Spalte: numerisch, definiert Gruppen fuer die eine gleiche Trennschaerfe angenommen werden soll ("restringierte 2pl-Modelle"); siehe auch Hilfe zu TAM
### guessMat          ... Matrix; erste Spalte: Item-ID, zweite Spalte: numerisch, definiert Gruppen mit gleichem Guessingparameter, unterschiedliche Zahlen bedeuten, dass jeweils ein separater
###                       Rateparameter geschaetzt wird. 0 oder fehlender Eintrag: kein Rateparameter wird geschaetzt
defineModel <- function(dat, items, id, splittedModels = NULL, irtmodel = c("1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", "2PL.groups", "GPCM.design", "3PL"),
               qMatrix=NULL, DIF.var=NULL, HG.var=NULL, group.var=NULL, weight.var=NULL, anchor = NULL, domainCol=NULL, itemCol=NULL, valueCol=NULL,check.for.linking = TRUE,
               minNperItem = 50, removeMinNperItem = FALSE, boundary = 6, remove.boundary = FALSE, remove.no.answers = TRUE, remove.no.answersHG = TRUE, remove.missing.items = TRUE, remove.constant.items = TRUE,
               remove.failures = FALSE, remove.vars.DIF.missing = TRUE, remove.vars.DIF.constant = TRUE, verbose=TRUE, software = c("conquest","tam"), dir = NULL,
               analysis.name, schooltype.var = NULL, model.statement = "item",  compute.fit = TRUE, pvMethod = c("regular", "bayesian"), fitTamMmlForBayesian = TRUE, n.plausible=5, seed = NULL, conquest.folder=system.file("exec", "console_Feb2007.exe", package = "eatModel"),
               constraints=c("cases","none","items"),std.err=c("quick","full","none"), distribution=c("normal","discrete"), method=c("gauss", "quadrature", "montecarlo", "quasiMontecarlo"),
               n.iterations=2000,nodes=NULL, p.nodes=2000, f.nodes=2000,converge=0.001,deviancechange=0.0001, equivalence.table=c("wle","mle","NULL"), use.letters=FALSE,
               allowAllScoresEverywhere = TRUE, guessMat = NULL, est.slopegroups = NULL, fixSlopeMat = NULL, slopeMatDomainCol=NULL, slopeMatItemCol=NULL, slopeMatValueCol=NULL,
               progress = FALSE, Msteps = NULL, increment.factor=1 , fac.oldxsi=0, export = list(logfile = TRUE, systemfile = FALSE, history = TRUE, covariance = TRUE, reg_coefficients = TRUE, designmatrix = FALSE) )   {
                  ismc <- attr(dat, "multicore")                                ### findet der Aufruf innerhalb einer multicore session statt?
                  dat  <- eatTools::makeDataFrame(dat, name = "dat")
     ### Sektion 'multiple models handling': jedes Modell einzeln von 'defineModel' aufbereiten lassen
     ### Hier wird jetzt erstmal nur die bescheuerte Liste aus 'splitModels' aufbereitet (wenn der Nutzer sie verhunzt hat)
     ### das findet natuerlich nur statt, wenn es 'splitModels' gibt, wenn also MEHRERE Modelle simultan verarbeitet werden sollen
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
     ### ok, die bescheuerte Liste aus 'splitModels' ist aufbereitet. Jetzt wird, sofern mehrere Modelle spezifiziert wurden, 'defineModel' mehrmals
     ### hintereinander aufgerufen (= rekursiv). Das geschieht ueber die Funktion 'doAufb' und kann prinzipiell ueber single- oder multicore erfolgen
     ### WICHTIG: 'doAufb' ruft irgendwann 'defineModel' auf, braucht also alle Argumente dieser Funktion ... sonst werden defaults genommen, und das will man ja hier nicht mehr
     ### in einem ersten Schritt muessen also alle Argumente von 'default Models' gesammelt werden
                     cl1 <- as.list(match.call(definition = defineModel))       ### sammle hier alle Argumente, die der Nutzer bei Aufruf von 'defineModel' selbststaendig definiert hat!
                     cl1[["analysis.name"]] <- NULL                             ### aber 'analysis.name' muss weg, selbst wenn der Nutzer das definiert hat, denn es wird ggf. neu konstruiert bzw. erweitert durch den splitter
                     cl1[["dat"]] <- as.name("dat")                             ### Hotfix, keine Ahnung, weswegen das frueher auch ohne ging, jetzt aber eine Fehlermeldung gibt
                     cl1[["splittedModels"]] <- splittedModels                  ### ersetze nun alle Argumente, die die Funktion bis hierher geaendert hat!
                     cll <- list()                                              ### boah, wieso geht das hier nicht mehr?!? cll <- lapply ( cl1[2:length(cl1)], eval )
                     for ( u in 2:length(cl1)) {cll[[u-1]] <- eval(cl1[[u]])}
                     names(cll) <- names(cl1)[-1]
                     anf <- mods[1]
     ### single core handling: Funktion "doAufb" wird seriell fuer alle "mods" aufgerufen
                     if(is.null ( splittedModels[["nCores"]] ) | splittedModels[["nCores"]] == 1 ) {
                        models <- lapply ( mods, FUN = doAufb, matchCall = cll, anf=anf, verbose=verbose, dir=dir, multicore=FALSE)
     ### wenn multicore handling, dann wird das Objekt "model" an cores verteilt und dort weiter verarbeitet. Ausserdem werden Konsolenausgaben in das stringobjekt "txt" weitergeleitet
                     }  else  {
                        txt <- capture.output ( models <- lapply ( mods, FUN = doAufb, matchCall = cll, anf=anf, verbose=verbose, dir=dir, multicore=TRUE) )
                        # if(!exists("detectCores"))   {library(parallel)}
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
                        cat(paste ( length(models), " models were prepared for estimation: ", sep="")); print( Sys.time() - beg)
     ### Trenne Aufbereitungsergebnisse von Konsolennachrichten
                        models <- lapply(mods, FUN = function ( m ) { m[["res"]] } )
     ### multicore gibt keine Ausgaben auf die Konsole, die muessen ueber "capture.output" eingefangen und separat ausgegeben werden
                        txts<- lapply(mods, FUN = function ( m ) { m[["txt"]] } )
                        luec<- which(txt == "")
                        pos <- luec[which ( diff(luec) == 1L )]
                        dif2<- which(diff(pos) == 1L)                           ### Hotfix!
                        if(length(dif2)>0) { pos <- pos [ -dif2 ] }
                        pos <- c(pos, length(txt)+1)
                        txtP<- lapply ( 1:(length(pos)-1), FUN = function ( u ) { txt[ pos[u] : (pos[u+1]-1) ] })
                        txtG<- NULL
                        stopifnot(length(txtP) == length(txts))
                        for ( j in 1:length(txtP) ) {
                              txtG <- c(txtG, txtP[[j]], txts[[j]])
                        }                                                       ### Hotfix 2: Zusaetzliche warnungen aus 'doAufb' abfangen
                        fl  <- min(grep( pattern = "====", x = txt))            ### 'fl' = first line
                        if ( fl > 1) {
                             if (!all(txt[1:(fl-1)] == "" )) {
                                 txtG <- c("", txt[1:(fl-3)], txtG)
                             }
                        }
                        cat(txtG, sep="\n")
                     }
                  attr(models, "split") <- splittedModels
                  class(models)    <- c("defineMultiple", "list")
                  return(models)                                                ### Das ist die Rueckgabe fuer den Mehrmodellfall
                  }  else  {
     ### ACHTUNG: hier beginnt jetzt der 'single model Fall' von 'defineModel' ###
                     irtmodel <- match.arg(irtmodel)
                     if ( is.null(Msteps) ) {                                   ### den Default fuer Msteps so setzen wie in TAM
                          if ( irtmodel == "3PL" ) { Msteps <- 10 } else { Msteps <- 4 }
                     }
                     software <- match.arg(software)
                     method   <- match.arg(method)
                     pvMethod <- match.arg(pvMethod)
                     if(software == "conquest") {
                        original.options <- options("scipen")                   ### lese Option fuer Anzahl der Nachkommastellen
                        options(scipen = 20)                                    ### setze Option fuer Anzahl der Nachkommastellen
                        if(missing(analysis.name)) {stop("Please specify 'analysis.name' or use 'software = \"tam\"'\n")}
                     }  else  {
                        if(missing(analysis.name)) {analysis.name <- "not_specified"}
                     }
                     if(length(model.statement)!=1)                {stop("'model.statement' has to be of length 1.\n")}
                     if(!inherits(model.statement, "character"))   {stop("'model.statement' has to be of class 'character'.\n")}
                     if(missing(dat))   {stop("No dataset specified.\n") }      ### 11.04.2014: nutzt Hilfsfunktionen von repMean etc.
                     if(is.null(items)) {stop("Argument 'items' must not be NULL.\n",sep="")}
                     if(length(items) == 0 ) {stop("Argument 'items' has no elements.\n",sep="")}
                     if ( length(items) != length(unique(items)) ) {
                          warning(paste0("Warning: Item identifier is not unique. Only ",length(unique(items))," unique item identifiers will be used."))
                          items <- unique(items)
                     }
                     if(length(id) != 1 ) {stop("Argument 'id' must be of length 1.\n",sep="")}
                     allVars     <- list(ID = id, variablen=items, DIF.var=DIF.var, HG.var=HG.var, group.var=group.var, weight.var=weight.var, schooltype.var = schooltype.var)
                     all.Names   <- lapply(allVars, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = dat, variable=ii)})
     ### wenn software = conquest, duerfen variablennamen nicht mehr als 11 Zeichen haben!
                     if(software == "conquest") {
                         if(max(nchar(all.Names[["variablen"]]))>11) {stop("In Conquest, maximum length of variable names must not exceed 11 characters. Please shorten variables names.\n")}
                     }
     ### ID-Variable pruefen und ggf. aendern
                     dat <- checkID_consistency(dat=dat, allNam=all.Names, software=software)
     ### Verzeichnis ('dir') pruefen oder erzeugen
                     dir <- checkDir(dir=dir, software=software)
     ### pruefen, ob es Personen gibt, die weniger als <boundary> items gesehen haben (muss VOR den Konsistenzpruefungen geschehen)
                     dat <- checkBoundary(dat=dat, allNam=all.Names, boundary=boundary, remove.boundary=remove.boundary)
     ### Sektion 'explizite Variablennamen ggf. aendern' ###
                     subsNam <- .substituteSigns(dat=dat, variable=unlist(all.Names[-c(1:2)]), all.Names = all.Names)
                     if(software == "conquest" || !is.null(all.Names[["DIF.var"]])) {
                        if(!all(subsNam$old == subsNam$new)) {                  ### Conquest erlaubt keine gross geschriebenen und expliziten Variablennamen, die ein "." oder "_" enthalten
                           sn     <- subsNam[which( subsNam$old != subsNam$new),]
                           cat("Conquest neither allows '.', '-', and '_' nor upper case letters in explicit variable names. Delete signs from variables names for explicit variables.\n"); flush.console()
                           recStr <- paste("'",sn[,"old"] , "' = '" , sn[,"new"], "'" ,sep = "", collapse="; ")
                           colnames(dat) <- car::recode(colnames(dat), recStr)
                           all.Names     <- lapply(all.Names, FUN = function ( y ) { car::recode(y, recStr) })
                           if(model.statement != "item") {
                              cat("    Remove deleted signs from variables names for explicit variables also in the model statement. Please check afterwards for consistency!\n")
                              for ( uu in 1:nrow(sn))  {model.statement <- gsub(sn[uu,"old"], sn[uu,"new"], model.statement)}
                           }
                        }
                        if("item" %in% unlist(all.Names[-c(1:2)])) { stop("Conquest does not allow labelling explicit variable(s) with 'Item' or 'item'.\n") }
                     }                                                          ### untere Zeilen: Dif-Variablen und Testitems duerfen sich nicht ueberschneiden
                     if(length(intersect(all.Names$DIF.var, all.Names$variablen))>0)    {stop("Test items and DIF variable have to be mutually exclusive.\n")}
                     if(length(intersect(all.Names$weight.var, all.Names$variablen))>0) {stop("Test items and weighting variable have to be mutually exclusive.\n")}
                     if(length(intersect(all.Names$HG.var, all.Names$variablen))>0)     {stop("Test items and HG variable have to be mutually exclusive.\n")}
                     if(length(intersect(all.Names$group.var, all.Names$variablen))>0)  {stop("Test items and group variable have to be mutually exclusive.\n")}
     ### Sektion 'Q matrix ggf. erstellen und auf Konsistenz zu sich selbst und zu den Daten pruefen' ###
                     if(is.null(qMatrix)) { qMatrix <- data.frame ( item = all.Names$variablen, Dim1 = 1, stringsAsFactors = FALSE) } else {
                         qMatrix <- checkQmatrixConsistency(qMatrix)            ### pruefe Konsistenz der q-matrix
                         notInDat<- setdiff(qMatrix[,1], all.Names$variablen)
                         notInQ  <- setdiff( all.Names$variablen , qMatrix[,1])
                         if(length(notInDat)>0) {
                            cat(paste("Following ", length(notInDat)," item(s) missed in data frame will be removed from Q matrix: \n    ",paste(notInDat,collapse=", "),"\n",sep=""))
                            qMatrix <- qMatrix[-match(notInDat, qMatrix[,1]),]
                            if(nrow(qMatrix) == 0) { stop("No common items in Q matrix and data.\n")}
                         }
                         if(length(notInQ)>0) {
                            cat(paste("Following ", length(notInQ)," item(s) missed in Q matrix will be removed from data: \n    ",paste(notInQ,collapse=", "),"\n",sep=""))
                         }                                                      ### Wichtig! Sicherstellen, dass Reihenfolge der Items in Q-Matrix mit Reihenfolge der Items im Data.frame uebereinstimmt!
                         all.Names[["variablen"]] <- qMatrix[,1]  } ;   flush.console()
     ### Sektion 'Alle Items auf einfache Konsistenz pruefen'
                      cic <- checkItemConsistency(dat=dat, allNam = all.Names, remove.missing.items=remove.missing.items, verbose=verbose, removeMinNperItem=removeMinNperItem, minNperItem=minNperItem, remove.constant.items=remove.constant.items, model.statement=model.statement)
     ### Sektion 'Hintergrundvariablen auf Konsistenz zu sich selbst und zu den Itemdaten pruefen'. Ausserdem Stelligkeit (Anzahl der benoetigten character) fuer jede Variable herausfinden
                      cbc <- checkBGV(allNam = cic[["allNam"]], dat=cic[["dat"]], software=software, remove.no.answersHG=remove.no.answersHG, remove.vars.DIF.missing=remove.vars.DIF.missing, namen.items.weg=cic[["namen.items.weg"]], remove.vars.DIF.constant=remove.vars.DIF.constant)
     ### Sektion 'Itemdatensatz zusammenbauen' (fuer Conquest ggf. mit Buchstaben statt Ziffern)
                      if(length(cbc[["namen.items.weg"]])>0)  {
                         cat(paste("Remove ",length(unique(cbc[["namen.items.weg"]]))," test item(s) overall.\n",sep=""))
                         cbc[["allNam"]]$variablen <- setdiff(cbc[["allNam"]]$variablen, unique(cbc[["namen.items.weg"]]) )
                         qMatrix             <- qMatrix[match(cbc[["allNam"]]$variablen, qMatrix[,1]),]
                      }
     ### Sektion 'Personen ohne gueltige Werte identifizieren und ggf. loeschen'. Gibt dat, perNA, datL zurueck
                      pwvv<- personWithoutValidValues(dat=cbc[["dat"]], allNam=cbc[["allNam"]], remove.no.answers=remove.no.answers)
     ### Sektion 'Summenscores fuer Personen pruefen'
                      cpsc<- checkPersonSumScores(datL = pwvv[["datL"]], allNam = cbc[["allNam"]], dat=pwvv[["dat"]], remove.failures=remove.failures)
     ### Sektion 'Verlinkung pruefen'
                      if(check.for.linking == TRUE) {                           ### Dies geschieht auf dem nutzerspezifisch reduzierten/selektierten Datensatz
                         linkNaKeep <- checkLink(dataFrame = cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], remove.non.responser = FALSE, verbose = FALSE )
                         linkNaOmit <- checkLink(dataFrame = cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], remove.non.responser = TRUE, verbose = FALSE )
                         if(linkNaKeep == FALSE & linkNaOmit == TRUE )  {cat("Note: Dataset is not completely linked. This is probably only due to missings on all cases.\n")}
                         if(linkNaKeep == TRUE )                        {cat("Dataset is completely linked.\n")}
                      }
     ### Sektion 'Anpassung der Methode (gauss, monte carlo) und der nodes'
                      met <- adaptMethod(method=method, software=software, nodes=nodes)
     ### Sektion 'Datensaetze softwarespezifisch aufbereiten: Conquest' ###
                      if(length(cbc[["namen.all.hg"]])>0) {all.hg.char <- sapply(cbc[["namen.all.hg"]], FUN=function(ii) {max(nchar(as.character(na.omit(cpsc[["dat"]][,ii]))))})} else {all.hg.char <- NULL}
                      if ( software == "conquest" )   {                         ### untere Zeile: wieviele character muss ich fuer jedes Item reservieren?
                          var.char  <- sapply(cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], FUN=function(ii) {max(nchar(as.character(na.omit(ii))))})
                          no.number <- setdiff(1:length(var.char), grep("[[:digit:]]",var.char))
                          if(length(no.number)>0) {var.char[no.number] <- 1}    ### -Inf steht dort, wo nur missings sind, hier soll die Characterbreite auf 1 gesetzt sein
                          if(use.letters == TRUE)   {                           ### sollen Buchstaben statt Ziffern benutzt werden? Dann erfolgt hier Recodierung.
                             rec.statement <- paste(0:25,"='",LETTERS,"'",sep="",collapse="; ")
                             for (i in cbc[["allNam"]][["variablen"]])  {             ### Warum erst hier? Weil Pruefungen (auf Dichotomitaet etc. vorher stattfinden sollen)
                                  cpsc[["dat"]][,i] <- car::recode(cpsc[["dat"]][,i], rec.statement)}
                             var.char <- rep(1,length(cbc[["allNam"]][["variablen"]]))}## var.char muss nun neu geschrieben werden, da nun alles wieder einstellig ist!
                      }
     ### Sektion 'deskriptive Ergebnisse berechnen und durchschleifen' ###
                      daten   <- data.frame(ID=as.character(cpsc[["dat"]][,cbc[["allNam"]][["ID"]]]), cpsc[["dat"]][,cbc[["namen.all.hg"]], drop = FALSE], cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], stringsAsFactors = FALSE)
                      deskRes <- desk.irt(daten = daten, itemspalten = match(cbc[["allNam"]][["variablen"]], colnames(daten)), percent = TRUE)
                      crit    <- which (deskRes[,"valid"] < minNperItem)
                      if ( length(crit)>0) {
                           cat ( paste ( "Following ",length(crit), " items with less than ",minNperItem," item responses:\n",sep=""))
                           options(width=1000)
                           print(deskRes[crit,-match(c("item.nr", "Label", "KB", "Codes", "Abs.Freq", "Rel.Freq"), colnames(deskRes))], digits = 3)
                      }
                      discrim <- item.diskrim(daten,match(cbc[["allNam"]][["variablen"]], colnames(daten)))
                      if ( length ( cbc[["allNam"]][["schooltype.var"]] ) > 0 ) { ### jetzt ggf. noch schulformspezifische p-Werte, falls gewuenscht
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
     ### erstmal testen, ob die Characterzahl wirklich einheitlich ist ... datensatz wird dazu nicht auf festplatte geschrieben
                          txt  <-  capture.output ( gdata::write.fwf(daten , colnames = FALSE,rownames = FALSE, sep="",quote = FALSE,na=".", width=fixed.width))
                          stopifnot(length(table(nchar(txt)))==1)               ### Check: hat der Resultdatensatz eine einheitliche Spaltenanzahl? Muss unbedingt sein!
                          rm(txt)                                               ### Speicher sparen
                          gdata::write.fwf(daten , file.path(dir,paste(analysis.name,".dat",sep="")), colnames = FALSE,rownames = FALSE, sep="",quote = FALSE,na=".", width=fixed.width)
                          colnames(lab) <- c("===>","item")                     ### schreibe Labels!
                          write.table(lab,file.path(dir,paste(analysis.name,".lab",sep="")),col.names = TRUE,row.names = FALSE, dec = ",", sep = " ", quote = FALSE)
                          if(!is.null(conquest.folder))     {
                             batch <- paste( normalize.path(conquest.folder),paste(analysis.name,".cqc",sep=""), sep=" ")
                             write(batch, file.path(dir,paste(analysis.name,".bat",sep="")))}
                          foo <- gen.syntax(Name=analysis.name, daten=daten, all.Names = cbc[["allNam"]], namen.all.hg = cbc[["namen.all.hg"]], all.hg.char = all.hg.char, var.char= max(var.char), model=qMatrix, anchored=anchor, pfad=dir, n.plausible=n.plausible, compute.fit = compute.fit,
                                            constraints=constraints, std.err=std.err, distribution=distribution, method=met[["method"]], n.iterations=n.iterations, nodes=met[["nodes"]], p.nodes=p.nodes, f.nodes=f.nodes, converge=converge,deviancechange=deviancechange, equivalence.table=equivalence.table, use.letters=use.letters, model.statement=model.statement, conquest.folder = conquest.folder, allowAllScoresEverywhere = allowAllScoresEverywhere, seed = seed, export = export)
                          if(!is.null(anchor))  {
                             write.table(ankFrame[["resConquest"]], file.path(dir,paste(analysis.name,".ank",sep="")) ,sep=" ", col.names = FALSE, row.names = FALSE, quote = FALSE)
                          }
     ### wenn Conquest gewaehlt, dann ggf. Logfile umbenennen, falls es bereits (unter demselben namen) existiert
                          if(file.exists( file.path ( dir,  paste(analysis.name,".log",sep=""))) )  {
                             cat(paste("Found existing log file '",paste(analysis.name,".log",sep=""), "' in folder '",dir,"'\nConquest analysis will overwrite log file. Original log file will be saved as '",paste(analysis.name,"_old.log'\n",sep=""),sep=""))
                             do <- file.rename(from = file.path(dir, paste(analysis.name,".log",sep="")), to = file.path(dir, paste(analysis.name,"_old.log",sep="")))
                          }
     ### Sektion 'Rueckgabeobjekt bauen', hier fuer Conquest                    ### setze Optionen wieder in Ausgangszustand
                          options(scipen = original.options); flush.console()   ### Achtung: setze Konsolenpfade in Hochkommas, da andernfalls keine Leerzeichen in den Ordner- bzw. Dateinamen erlaubt sind!
                          ret <- list ( software = software, input = paste("\"", file.path(dir, paste(analysis.name,"cqc",sep=".")), "\"", sep=""), conquest.folder = paste("\"", conquest.folder, "\"", sep=""), dir=dir, analysis.name=analysis.name, model.name = analysis.name, qMatrix=qMatrix, all.Names=cbc[["allNam"]], deskRes = deskRes, discrim = discrim, perNA=pwvv[["perNA"]], per0=cpsc[["per0"]], perA = cpsc[["perA"]], perExHG = cbc[["perExHG"]], itemsExcluded = cbc[["namen.items.weg"]], daten=daten)
                          class(ret) <-  c("defineConquest", "list")
                          return ( ret )  }
     ### Sektion 'Rueckgabeobjekt fuer tam'
                      if ( software == "tam" )   {
                          cat(paste("Q matrix specifies ",ncol(qMatrix)-1," dimension(s).\n",sep=""))
                          anchor          <- prepAnchorTAM(ank = ankFrame[["resTam"]], allNam = cbc[["allNam"]])
                          est.slopegroups <- prepEstSlopegroupsTAM(esg = est.slopegroups, allNam = cbc[["allNam"]])
                          fixSlopeMat     <- prepFixSlopeMatTAM(fsm = fixSlopeMat, allNam = cbc[["allNam"]], qma =  qMatrix, slopeMatDomainCol=slopeMatDomainCol, slopeMatItemCol=slopeMatItemCol, slopeMatValueCol=slopeMatValueCol, dat=daten)
                          guessMat        <- prepGuessMat(guessMat, allNam = fixSlopeMat[["allNam"]])
                          control <- list ( snodes = met[["snodes"]] , QMC=met[["QMC"]], convD = deviancechange ,conv = converge , convM = .0001 , Msteps = Msteps , maxiter = n.iterations, max.increment = 1 ,
                                     min.variance = .001 , progress = progress , ridge=0 , seed = seed , xsi.start0=FALSE,  increment.factor=increment.factor , fac.oldxsi= fac.oldxsi)
                          if ( !is.null(met[["nodes"]])) { control$nodes <- met[["nodes"]] }
                          ret     <- list ( software = software, constraint = match.arg(constraints) , qMatrix=qMatrix, anchor=anchor,  all.Names=fixSlopeMat[["allNam"]], daten=daten, irtmodel=irtmodel, est.slopegroups = est.slopegroups, guessMat=guessMat, control = control, n.plausible=n.plausible, dir = dir, analysis.name=analysis.name, deskRes = deskRes, discrim = discrim, perNA=pwvv[["perNA"]], per0=cpsc[["per0"]], perA = cpsc[["perA"]], perExHG = cbc[["perExHG"]], itemsExcluded = cbc[["namen.items.weg"]], fixSlopeMat = fixSlopeMat[["slopMat"]], estVar = fixSlopeMat[["estVar"]], pvMethod = pvMethod,  fitTamMmlForBayesian=fitTamMmlForBayesian)
                          class(ret) <-  c("defineTam", "list")
                          return ( ret )    }   }  }

### Hilfsfunktion fuer defineModel() ... scheint buggy zu sein
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

### Hilfsfunktion fuer defineModel()
prepFixSlopeMatTAM <- function (fsm, allNam, qma, slopeMatDomainCol, slopeMatItemCol, slopeMatValueCol, dat){
       if(!is.null(fsm))  {                                                     ### Achtung: wenn Items identifiers NICHT unique sind (z.B., Item gibt es global und domaenenspezifisch, dann wird jetzt 'fixSlopeMat' auf die Dimension in der Q Matrix angepasst ... das ist nur erlaubt, wenn es ein eindimensionales Modell ist!!
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
              cat(paste("Following ",length(weg2), " items in matrix for items with fixed slopes ('fixSlopeMat') which are not in dataset:\n",sep=""))
              cat("   "); cat(paste(weg2, collapse=", ")); cat("\n")
              cat("Remove these item(s) from 'fixSlopeMat' matrix.\n")
              fsm <- fsm[-match(weg2,fsm[,1]),]
           }
           weg3   <- setdiff(allNam[["variablen"]], fsm[,1])
           if(length(weg3)>0) {
              cat(paste("Following ",length(weg3), " items in dataset without fixed slopes in 'fixSlopeMat'. Slope(s) will be estimated freely.\n",sep=""))
              cat("   "); cat(paste(weg3, collapse=", ")); cat("\n")
           }
     ### Achtung, grosser Scheiss: wenn man nicht (wie oben) eine Reihenfolgespalte angibt, aendert die untere 'by'-Schleife die Sortierung!
           if ( nrow(fsm) != length(unique(fsm[,1])) ) { stop ( "Item identifiers in 'fixSlopeMat' are not unique.\n")}
           fsm[,"reihenfolge"] <- 1:nrow(fsm)
           dims  <- (1:ncol(qma))[-1]                                           ### Slopematrix muss itemweise zusammengebaut werden
           slopMa<- do.call("rbind", by ( data = fsm, INDICES = fsm[,"reihenfolge"], FUN = function (zeile ) {
                    zeile <- zeile[,-ncol(zeile)]
                    stopifnot ( nrow(zeile) == 1 )
                    qSel  <- qma[which( qma[,1] == zeile[[1]]),]
                    anzKat<- length(unique(na.omit(dat[,as.character(zeile[[1]])])))
                    zeilen<- anzKat * length(dims)                              ### fuer jedes Items gibt es [Anzahl Kategorien] * [Anzahl Dimensionen] Zeilen in der TAM matrix
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
       return(list(allNam=allNam, estVar=estVar, slopMat = slopMa))}

### Hilfsfunktion fuer defineModel()
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
           }                                                                    ### untere zeile: pruefen, ob keine fehlenden oder leeren Eintraege in der Liste sind
           weg3 <- c(which(is.na(esg[,2])), which(esg[,2] ==""))
           if(length(weg3)>0) {stop("Items in 'est.slopegroups' with missing or empty values.\n")}
           esg  <- as.numeric(as.factor(as.character(esg[match(allNam[["variablen"]], esg[,1]),2])))
       }
       return(esg)}

### Hilfsfunktionen fuer defineModel
checkContextVars <- function(x, varname, type = c("weight", "DIF", "group", "HG"), itemdata, suppressAbort = FALSE, internal = FALSE)   {
                     type <- match.arg(arg = type, choices = c("weight", "DIF", "group", "HG"))
                     stopifnot(length(x) == nrow(itemdata))
                     if(missing(varname))  {varname <- "ohne Namen"}
                     if(!inherits(x, "numeric") && isTRUE(internal))  {         ### ist Variable numerisch?
                        if (type == "weight") {stop(paste(type, " variable has to be 'numeric' necessarily. Automatic transformation is not recommended. Please transform by yourself.\n",sep=""))}
                        cat(paste(type, " variable has to be 'numeric'. Variable '",varname,"' of class '",class(x),"' will be transformed to 'numeric'.\n",sep=""))
                        x <- suppressWarnings(unlist(eatTools::asNumericIfPossible(x = data.frame(x, stringsAsFactors = FALSE), transform.factors = TRUE, maintain.factor.scores = FALSE, force.string = FALSE)))
                        if(!inherits(x, "numeric"))  {                          ### erst wenn asNumericIfPossible fehlschlaegt, wird mit Gewalt numerisch gemacht, denn fuer Conquest MUSS es numerisch sein
                           x <- as.numeric(as.factor(x))
                        }
                        cat(paste("    '", varname, "' was converted into numeric variable of ",length(table(x))," categories. Please check whether this was intended.\n",sep=""))
                        if(length(table(x)) < 12 ) { cat(paste("    Values of '", varname, "' are: ",paste(names(table(x)), collapse = ", "),"\n",sep=""))}
                     }
                     toRemove<- NULL
                     mis     <- length(table(x))
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
                                        y       <- paste0("V", x)               ### wenn x numerisch ist, sind die Spaltennamen in completeMissingGroupwise nicht mehr den levels von x zuweisbar, da haengt R dann ein X ran
                                        n.werte <- lapply(itemdata, FUN=function(iii){by(iii, INDICES=list(y), FUN=table)})
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


### Hilfsfunktion fuer defineModel
prepAnchorTAM <- function (ank, allNam) {
        if(!is.null(ank)) {
            stopifnot(ncol(ank) == 2 )                                          ### Untere Zeile: Wichtig! Sicherstellen, dass Reihenfolge der Items in Anker-Statement der Reihenfolge im datensatz entspricht
            notInData   <- setdiff(ank[,1], allNam[["variablen"]])              ### messages entfernt, denn die werden ja schon in anker() durch mergeAttr() ausgegeben
            if(length(notInData)>0)  {ank <- ank[-match(notInData, ank[,1]),]}
            ank[,1]    <- match(as.character(ank[,1]), allNam[["variablen"]])
        }
        return(ank)}

### Hilfsfunktion fuer defineModel
checkBGV <- function(allNam, dat, software, remove.no.answersHG, remove.vars.DIF.missing, namen.items.weg, remove.vars.DIF.constant){
            weg.dif <- NULL; weg.hg <- NULL; weg.weight <- NULL; weg.group <- NULL# initialisieren
     ### Gibt es ueberhaupt irgendwelche Kovariaten?
            if(length(allNam[["HG.var"]])>0 || length(allNam[["group.var"]])>0 || length(allNam[["DIF.var"]])>0 || length(allNam[["weight.var"]]) >0  ) {
               varClass<- sapply(c(allNam[["HG.var"]],allNam[["group.var"]],allNam[["DIF.var"]], allNam[["weight.var"]]),FUN = function(ii) {class(dat[,ii])})
               if ( isFALSE(all(sapply(varClass, length) == 1)) ) {
                    fehler <- which(sapply(varClass, length) != 1)
                    cat(paste0("Following ",length(fehler), " variables with more that one class:"))
                    print(varClass[names(fehler)]); stop()
               }
            }
     ### Hintergrundvariablen (conditioning model)
            if(length(allNam[["HG.var"]])>0)    {
               varClass<- sapply(allNam[["HG.var"]], FUN = function(ii) {class(dat[,ii])})
               notNum  <- which(varClass %in% c("factor", "character"))
               if(length(notNum)>0) {
                  cat(paste("Warning: Background variables '",paste(names(varClass)[notNum], collapse="', '"),"' of class \n    '",paste(varClass[notNum],collapse="', '"),"' will be converted to indicator variables.\n",sep=""))
                  ind <- do.call("cbind", lapply ( names(varClass)[notNum], FUN = function ( yy ) {
                         if ( length(which(is.na(dat[,yy])))>0) { stop(paste0("Found ",length(which(is.na(dat[,yy]))), " missings on background variable '",yy,"'."))}
                         newFr <- model.matrix( as.formula (paste("~",yy,sep="")), data = dat)[,-1,drop=FALSE]
                         cat(paste("    Variable '",yy,"' was converted to ",ncol(newFr)," indicator(s) with name(s) '",paste(colnames(newFr), collapse= "', '"), "'.\n",sep=""))
                         return(newFr) }))
                  if(software == "conquest") {                                  ### ggf. fuer Conquest Namen der HG-Variablen aendern
                      subNm <- .substituteSigns(dat=ind, variable=colnames(ind))
                      if(!all(subNm$old == subNm$new)) {
                          sn  <- subNm[which( subNm$old != subNm$new),]
                          reSt<- paste("'",sn[,"old"] , "' = '" , sn[,"new"], "'" ,sep = "", collapse="; ")
                          colnames(ind) <- car::recode(colnames(ind), reSt)     ### wenn background-variable ein faktor ist, werden daraus jetzt numerische dummies
                      }                                                         ### entferne Originalnamen aus allNam[["HG.var"]] und ergaenze neue Namen
                  }                                                             ### ergaenze neue Variablen im Datensatz
                  allNam[["HG.var"]] <- setdiff ( allNam[["HG.var"]], names(varClass)[notNum])
                  allNam[["HG.var"]] <- c(allNam[["HG.var"]], colnames(ind))
                  if ( length(allNam[["HG.var"]]) > 99 && software == "conquest" ) {
                       warning(paste0(length(allNam[["HG.var"]]), " background variables might be problematic in 'Conquest'. Recommend to use 'TAM' instead."))
                  }                                                             ### Warnung wenn mehr als 100 HG-Variablen und Conquest
               dat <- data.frame ( dat, ind, stringsAsFactors = FALSE )
               }
               hg.info <- lapply(allNam[["HG.var"]], FUN = function(ii) {checkContextVars(x = dat[,ii], varname=ii, type="HG", itemdata=dat[,allNam[["variablen"]], drop = FALSE], suppressAbort = TRUE, internal=TRUE )})
               for ( i in 1:length(hg.info)) { dat[, hg.info[[i]][["varname"]] ] <- hg.info[[i]]$x }
               wegVar  <- unlist(lapply(hg.info, FUN = function ( uu ) { uu[["toRemove"]] }))
               if(length(wegVar)>0) { allNam[["HG.var"]] <- setdiff ( allNam[["HG.var"]], wegVar) }
               weg.hg  <- unique(unlist(lapply(hg.info, FUN = function ( y ) {y$weg})))
               if(length(weg.hg)>0) {                                           ### untere Zeile: das removen geschieht erst etwas spaeter, wenn datensatz zusammengebaut ist
                   if ( remove.no.answersHG == TRUE ) {
                        cat(paste("Remove ",length(weg.hg)," cases with missings on at least one HG variable.\n",sep=""))
                   }  else  {
                        cat(paste(length(weg.hg)," cases with missings on at least one HG variable will be kept according to 'remove.no.answersHG = FALSE'.\n",sep=""))
                        weg.hg <- NULL
                   }
               }
            }
     ### Gruppenvariablen
            if(length(allNam$group.var)>0)  {
                group.info <- lapply(allNam$group.var, FUN = function(ii) {checkContextVars(x = dat[,ii], varname=ii, type="group", itemdata=dat[,allNam[["variablen"]], drop = FALSE], internal=TRUE)})
                for ( i in 1:length(group.info)) { dat[, group.info[[i]]$varname ] <- group.info[[i]]$x }
                weg.group  <- unique(unlist(lapply(group.info, FUN = function ( y ) {y$weg})))
                if(length(weg.group)>0)  {                                      ### untere Zeile: das removen geschieht erst etwas spaeter, wenn datensatz zusammengebaut ist
                    cat(paste("Remove ",length(weg.group)," cases with missings on group variable.\n",sep=""))
                }
            }
     ### DIF-Variablen
            if(length(allNam$DIF.var)>0)  {
                dif.info <- lapply(allNam$DIF.var, FUN = function(ii) {checkContextVars(x = dat[,ii], varname=ii, type="DIF", itemdata=dat[,allNam[["variablen"]], drop = FALSE], internal = TRUE)})
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
                if(length(weg.dif)>0)  {                                        ### untere Zeile: removen geschieht erst etwas spaeter, wenn datensatz zusammengebaut ist
                    cat(paste("Remove ",length(weg.dif)," cases with missings on DIF variable.\n",sep=""))
                }
            }
     ### Gewichtungsvariablen
            if(length(allNam$weight.var)>0)  {
                if(length(allNam$weight.var)!=1) {stop("Use only one weight variable.")}
                weight.info <- lapply(allNam$weight.var, FUN = function(ii) {checkContextVars(x = dat[,ii], varname=ii, type="weight", itemdata=dat[,allNam[["variablen"]], drop = FALSE], internal = TRUE)})
                for ( i in 1:length(weight.info)) { dat[, weight.info[[i]]$varname ] <- weight.info[[i]]$x }
                weg.weight  <- unique(unlist(lapply(weight.info, FUN = function ( y ) {y$weg})))
                if(length(weg.weight)>0) {                                      ### untere Zeile: remove geschieht erst etwas spaeter, wenn datensatz zusammengebaut ist
                    cat(paste("Remove ",length(weg.weight)," cases with missings on weight variable.\n",sep=""))
                }

            }                                                                   ### untere Zeile, Achtung: group- und DIF- bzw. group- und HG-Variablen duerfen sich ueberschneiden!
     ### jetzt alles rausschmeissen, was wegen irgendeines Grundes raus soll
            namen.all.hg <- unique(c(allNam$HG.var,allNam$group.var,allNam$DIF.var,allNam$weight.var))
            weg.all <- unique(c(weg.dif, weg.hg, weg.weight, weg.group))
            perExHG <- NULL
            if(length(weg.all)>0) {
               cat(paste("Remove",length(weg.all),"case(s) overall due to missings on at least one explicit variable.\n"))
               perExHG<- dat[weg.all, allNam[["ID"]] ]
               dat    <- dat[-weg.all,]
            }
            return(list(dat=dat, allNam=allNam, namen.items.weg=namen.items.weg,perExHG=perExHG, namen.all.hg=namen.all.hg))}

### Hilfsfunktion fuer defineModel
checkItemConsistency <- function(dat, allNam, remove.missing.items, verbose, removeMinNperItem, minNperItem, remove.constant.items, model.statement){
          namen.items.weg <- NULL                                               ### initialisieren
     ### Wandle NaN in NA, falls es welche gibt
          is.NaN <- do.call("cbind", lapply(dat[,allNam[["variablen"]], drop = FALSE], FUN = function (uu) { is.nan(uu) } ) )
          if(sum(is.NaN) > 0 ) {
             cat(paste("Found ",sum(is.NaN)," 'NaN' values in the data. Convert 'NaN' to 'NA'.\n",sep=""))
             for ( j in allNam[["variablen"]]) {
                   weg <- which ( is.nan(dat[,j] ))
                   if(length(weg)>0) {  dat[weg,j] <- NA }
             }
          }
     ### sind die responses numerisch bzw. stehen da Ziffern drin? (notfalls sowas wie as.character(1) )
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
     ### identifiziere Items ohne jegliche gueltige Werte
          if(length(n.mis) >0) {
             cat(paste("Serious warning: ",length(n.mis)," testitems(s) without any values.\n",sep=""))
             if(verbose == TRUE) {cat(paste("    ", paste(names(n.mis), collapse=", "), "\n", sep=""))}
             if(remove.missing.items == TRUE) {
                 cat(paste("Remove ",length(n.mis)," variable(s) due to solely missing values.\n",sep=""))
                 namen.items.weg <- c(namen.items.weg, names(n.mis))
             }
          }
     ### identifiziere Items mit Anzahl gueltiger Werte < minNperItem
          if ( removeMinNperItem == TRUE ) {                                    ### identifiziere Items mit weniger gueltigen Werte als in 'minNperItem' angegeben (nur wenn 'removeMinNperItem' = TRUE)
               nValid <- unlist(lapply(dat[,allNam[["variablen"]], drop = FALSE], FUN = function ( ii ) { length(na.omit ( ii )) }))
               below  <- which ( nValid < minNperItem )
               if ( length ( below ) > 0 ) {
                    cat (paste ( "Found ", length(below), " items with less than ", minNperItem, " valid responses. These items will be removed.\n", sep=""))
                    namen.items.weg <- unique ( c(namen.items.weg, names(below)))
               }
          }
     ### identifiziere konstante Items (Items ohne Varianz)
          constant <- which(n.werte == 1)
          if(length(constant) >0) {
             cat(paste("Warning: ",length(constant)," testitems(s) are constants.\n",sep=""))
             if(verbose == TRUE) {foo <- lapply(names(constant),FUN=function(ii) {cat(paste(ii,": ",names(table(dat[,ii])), " ... ",length(na.omit(dat[,ii]))," valid responses", sep="")); cat("\n")})}
             if(remove.constant.items == TRUE) {
                 cat(paste("Remove ",length(constant)," variable(s) due to solely constant values.\n",sep=""))
                 namen.items.weg <- c(namen.items.weg, names(constant))
             }
          }
     ### identifiziere alle Items, die nicht dichotom (="ND") sind
          n.rasch  <- which( !isDichot )                                        ### (aber nicht die, die bereits wegen konstanter Werte aussortiert wurden!)
          if(length(n.rasch) >0 )   {                                           ### also polytome Items oder Items, die mit 1/2 anstatt 0/1 kodiert sind
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

### Hilfsfunktion fuer defineModel
checkID_consistency <- function(dat, allNam, software){
          dat[,allNam[["ID"]] ] <- as.character(dat[,allNam[["ID"]] ])
          doppelt     <- which(duplicated(dat[,allNam[["ID"]]]))
          if(length(doppelt)>0)  {stop(paste( length(doppelt) , " duplicate IDs found!",sep=""))}
          if(software == "conquest") {
              notAllowed  <- grep("-|\\.", dat[,allNam[["ID"]] ])               ### fuer Conquest: unerlaubte Zeichen aus ID-Variable loeschen
              if ( length(notAllowed)>0) {
                   cat("Conquest neither allows '.' nor '-' in ID variable. Delete signs from ID variable.\n")
                   dat[,allNam[["ID"]] ] <- eatTools::removePattern(string = eatTools::removePattern(string=dat[,allNam[["ID"]] ], pattern="\\."), pattern = "-")
                   if ( length ( which(duplicated(dat[,allNam[["ID"]]])))>0) {
                         dat[,allNam[["ID"]] ] <- paste0(1:nrow(dat),dat[,allNam[["ID"]] ])
                   }                                                            ### wenn ID jetzt nicht mehr unique ist, wieder unique machen
              }
          }
          return(dat)}

### Hilfsfunktion fuer defineModel
checkDir <- function(dir, software) {
            if(!is.null(dir)) {                                                 ### Sofern ein verzeichnis angegeben wurde (nicht NULL),
                dir <- eatTools::crop(dir,"/")                                  ### das Verzeichnis aber nicht existiert, wird es jetzt erzeugt
                if(dir.exists(dir) == FALSE) {
                   cat(paste("Warning: Specified folder '",dir,"' does not exist. Create folder ... \n",sep=""))
                   dir.create(dir, recursive = TRUE)
                }
            }  else  {                                                          ### sicher ist sicher ...
                if (software == "conquest") {stop("Argument 'dir' must be specified if software = 'conquest'.\n")}
            }
            return(dir)}

### Hilfsfunktion fuer defineModel
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

### Hilfsfunktion fuer defineModel
personWithoutValidValues <- function (dat, allNam, remove.no.answers){
          if(inherits(try(datL  <- reshape2::melt(data = dat, id.vars = unique(unlist(allNam[-match("variablen", names(allNam))])), measure.vars = allNam[["variablen"]], na.rm=TRUE)  ),"try-error"))  {
             cat("W A R N I N G ! ! !   Error in melting for unknown reasons. Try workaround.\n"); flush.console()
             allHG <- setdiff(unique(unlist(allNam[-match("variablen", names(allNam))])), allNam[["ID"]] )
             stopifnot(length(allHG)>0)                                         ### dies ist ein Workaround, wenn "melt" fehltschlaegt (Fehler nicht reproduzierbar)
             datL  <- reshape2::melt(data = dat, id.vars = allNam[["ID"]], measure.vars = allNam[["variablen"]], na.rm=TRUE)
             datL  <- merge(datL, dat[,unique(unlist(allNam[-match("variablen", names(allNam))]))], by = allNam[["ID"]], all=TRUE)
          }
          wegNV <- setdiff(dat[,allNam[["ID"]]], unique(datL[,allNam[["ID"]]]))
          perNA <- NULL
          if(length(wegNV)>0)   {                                               ### identifiziere Faelle mit ausschliesslich missings
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

### Hilfsfunktion fuer defineModel
checkPersonSumScores <- function(datL, allNam, dat, remove.failures){
          minMax<- do.call("rbind", by ( data = datL, INDICES = datL[,"variable"], FUN = function ( v ) {
                   v[,"valueMin"] <- min(v[,"value"])                           ### obere Zeile: hier wird variablenweise der kleinstmoegliche Wert gesucht
                   v[,"valueMax"] <- max(v[,"value"])                           ### da der hier verwendete Longdatensatz 'datL' oben mit 'na.rm = TRUE' erzeugt wurde,
                   return(v)}))                                                 ### sind hier diejenigen Personen mit ausschliesslich Missings bereits eliminiert
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
             # alle<- na.omit(match(allTru, dat[,allNam[["ID"]]]))
             perA<- numT
          }
          return(list(dat=dat, per0=per0, perA=perA))}

### Hilfsfunktion fuer defineModel
adaptMethod <- function(method, software,nodes){
        snodes <- NULL; QMC <- NULL                                             ### initialisieren
        if(method == "quasiMontecarlo" && software == "conquest") {
           cat("Method 'quasiMontecarlo' is not available for software 'conquest'. Set method to 'montecarlo'.\n")
           method <- "montecarlo"
        }
        if(method %in% c("montecarlo", "quasiMontecarlo"))  {
           if(nodes < 500 ) {
              warning(paste0("Due to user specification, only ",nodes," nodes are used for '",method,"' estimation. Please note or re-specify your analysis."))
           }
           if(is.null(nodes) )   {
              cat(paste("'",method,"' has been chosen for estimation method. Number of nodes was not explicitly specified. Set nodes to 1000.\n",sep=""))
              if(software == "conquest") {nodes <- 1000}
              if(software == "tam" )     {nodes <- NULL; snodes <- 1000; QMC <- as.logical(car::recode ( method, "'montecarlo'=FALSE; 'quasiMontecarlo'=TRUE"))}
           }  else  {
              if(software == "tam" )     {
                 snodes <- nodes; nodes <- NULL; QMC <- as.logical(car::recode ( method, "'montecarlo'=FALSE; 'quasiMontecarlo'=TRUE"))
              }
           }
        }  else  {
           if ( is.null(nodes) )   {
                cat(paste("Number of nodes was not explicitly specified. Set nodes to 20 for method '",method,"'.\n",sep=""))
                if ( software == "conquest" ) { nodes <- 20; snodes <- NULL; QMC <- NULL }
  	            if ( software == "tam" )      { nodes <- seq(-6,6,len=20); snodes <- 0; QMC <- FALSE}
    			}  else {
                if ( software == "tam" )      { nodes <- seq(-6,6,len=nodes); snodes <- 0; QMC <- FALSE }
          }
    		}
    		return(list(method=method, nodes=nodes, snodes=snodes, QMC=QMC))}

.substituteSigns <- function(dat, variable, all.Names = NULL ) {
                    if(!is.null(variable)) {
           					   variableNew <- tolower(gsub("_|\\.|-", "", variable))
     ### DIF-variablen duerfen zusaetzlich keine Zahlen enthalten
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
    ### alle Spalten ausser der ersten muessen numerisch oder integer sein
             if ( !all(nClass)) {
                  warning(paste0("Found non-numeric indicator column(s) in the Q matrix. Transform column(s) '",paste(colnames(qmat)[ which(nClass==FALSE)+1], collapse = "', '") ,"' into numeric format."))
                  qmat <- data.frame ( qmat[,1,drop=FALSE], eatTools::asNumericIfPossible(qmat[,-1,drop=FALSE]), stringsAsFactors = FALSE)
             }
    ### es duerfen nur werte von 0 und 1 auftreten (keine missings)
             werte <- eatTools::tableUnlist(qmat[,-1,drop=FALSE], useNA="always")
             if(length(setdiff( names(werte) , c("0","1", "NA")))<0) {stop("Q matrix must not contain entries except '0' and '1'.\n")}
             if(werte[match("NA", names(werte))] > 0) {stop("Missing values in Q matrix.\n")}
    ### Indikatorspalten duerfen nicht konstant 0 sein (konstant 1 ginge, das waere dann within item multidimensionality)
             wertes<- lapply(qmat[,-1,drop=FALSE], FUN = function (col)  {all ( col == 0)})
             konst <- which(wertes == TRUE)
             if ( length(konst)>0) {
                  cat(paste0("Column(s) '",paste(names(konst), collapse = "', '"),"' in Q matrix are konstant with value 0. Delete column(s).\n"))
                  qmat <- qmat[,-match(names(konst), colnames(qmat)), drop=FALSE]
             }
    ### keine doppelten Eintraege in Itemspalte
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
    ### items loeschen, die auf keiner dimension laden
             zeilen<- apply(qmat, 1, FUN = function ( y ) { all ( names(table(y[-1])) == "0")  })
             weg   <- which(zeilen == TRUE)
             if(length(weg)>0) {
                cat(paste("Note: Following ",length(weg)," item(s) in Q matrix do not belong to any dimension. Delete these item(s) from Q matrix.\n",sep=""))
                cat("    "); cat(paste(qmat[weg,1],collapse=", ")); cat("\n")
                qmat  <- qmat[-weg,]
             }
             return(qmat)}


checkLink <- function(dataFrame, remove.non.responser = FALSE, sysmis = NA, verbose = TRUE)   {
             if(!is.na(sysmis))  {
               na <- which(is.na(dataFrame))
               if(length(na)>0)  {
                  warning(paste0("'",sysmis,"' was specified to denote 'sysmis' in the data. ",length(na)," 'NA'-values were found in the dataset anyway. \n         Hence, ",sysmis," and 'NA' will be handled as 'sysmis'."))
               }
               dataFrame <- as.data.frame(lapply(dataFrame, FUN=function(ii) {car::recode(ii, paste(sysmis,"= NA", collapse = "; ") ) } ) )
             }
             if ( remove.non.responser == TRUE ) {
                na <- which( rowSums(is.na(dataFrame)) == ncol ( dataFrame ) )
                if(length(na)>0) {
                   dataFrame <- dataFrame[-na,]
                   if(verbose == TRUE ) {cat(paste("Remove ",length(na)," cases with missing on all items.\n", sep = ""))}
                }
             }
             non.missing.cases <- lapply(dataFrame, FUN=function(ii) {which(!is.na(ii))})
             all.cases <- non.missing.cases[[1]]
             i <- 2
             total.abbruch     <- FALSE
             while( (i < length(non.missing.cases) + 1 ) & !total.abbruch )  {
                  if(length( intersect(all.cases,non.missing.cases[[i]])) > 0 )  {
                     all.cases <- unique(c(all.cases, non.missing.cases[[i]] ) )
                  }  else   {
                     overlap        <- FALSE
                     remain.columns <- length(non.missing.cases) + 1 - i
                     ii             <- 1
                     while (overlap == FALSE & ii < remain.columns )  {
                           non.missing.cases <- non.missing.cases[c(setdiff(1:length(non.missing.cases),i),i)]
                          if(length( intersect(all.cases,non.missing.cases[[i]])) > 0 ) {overlap <- TRUE}
                           ii <- ii + 1
                     }
                     if (overlap == FALSE) {total.abbruch <- TRUE}
                     if (overlap == TRUE)  {all.cases <- unique(c(all.cases, non.missing.cases[[i]] ) ) }
                  }
                  i <- i + 1
             }
             if (length(all.cases) != nrow(dataFrame))   {
                if (verbose == TRUE) {cat("WARNING! Dataset is not completely linked.\n") }
                if ( remove.non.responser == TRUE ) {
                     missed <- setdiff ( 1:nrow(dataFrame), all.cases)
                     misFra <- reshape2::melt ( data.frame ( id = 1:length(missed), dataFrame[missed,]), id.vars = "id", na.rm=TRUE)
                     cat ( paste ( "W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).\n                  ",length(missed)," cases unconnected. Following items are unconnected: \n",sep=""))
                     cat("                  "); cat ( paste ( unique(as.character(misFra[,"variable"])), collapse = ", ")); cat("\n")
                }
                return(FALSE)
             }
             if (length(all.cases) == nrow(dataFrame))   {
                if (verbose == TRUE) {cat("Dataset is completely linked.\n") }
                return(TRUE)
             }  }

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
    ### Achtung!! wenn das Modell polytom war, muss p-Wert aus 'itn'-File ausgelesen werden!
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
         shw3 <- shw31 <- NULL                                                  ### initialisieren
         if(is.null ( deskRes ) ) { return(NULL)}
         deskR<- merge(deskRes, qL[,-match("value", colnames(qL))], by.x = "item.name", by.y = colnames(qMatrix)[1], all=TRUE)
         if ( isPoly == FALSE ) {
               shw3 <- data.frame ( model = model.name, source = "conquest", var1 = deskR[,"item.name"], var2 = NA , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "itemP",  derived.par = NA, value = deskR[,"item.p"], stringsAsFactors = FALSE)
         }
         shw4 <- data.frame ( model = model.name, source = "conquest", var1 = deskR[,"item.name"], var2 = NA , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "Nvalid",  derived.par = NA, value = deskR[,"valid"], stringsAsFactors = FALSE)
         shw4 <- shw4[!duplicated(shw4[,"var1"]),]
    ### Achtung! wenn in dem 'deskRes'-Objekt noch mehr p-Werte (schulformspezifische p-Werte drinstehen, werden die jetzt auch in die Ergebnisstruktur eingetragen)
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
         if(length(shw) <= 4 )  {  return(NULL)}                                ### ggf. Parameter zusaetzlicher Conquest-Terme einlesen, wenn length(shw) <= 4, gibt es keinen zusaetzlichen Terme
         res   <- NULL                                                          ### initialisieren
         read  <- 2 : (length(shw) - 3)                                         ### Diese Terme muessen eingelesen werden
         for ( i in names(shw)[read] ) {
               cols <- unlist(isLetter(i))                                      ### versuche Spalte(n) zu identifizieren
               if( !all(cols %in% colnames(shw[[i]])) ) {
                   cat(paste("Cannot identify variable identifier for term '",i,"' in file '",shwFile,"'. Skip procedure.\n",sep=""))
               }  else  {
                   if(length(cols) == 1 ) {
                      var1 <- paste( cols, shw[[i]][,cols],sep="_")
                   } else {
                      var1 <- unlist(apply(shw[[i]][,cols], MARGIN=1, FUN = function ( y ) {
                              paste ( unlist(lapply ( 1:length(y), FUN = function ( yy ) { paste(names(y)[yy], y[yy],sep="_")})), sep="", collapse = "_X_")  }))
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
         if(ncol(qMatrix) == 2) {                                               ### eindimensionaler Fall
            res  <- data.frame ( model = model.name, source = "conquest", var1 = colnames(qMatrix)[2], var2 = NA , type = "distrpar", indicator.group = NA, group = "persons", par = "var",  derived.par = NA, value = shw$cov.structure, stringsAsFactors = FALSE)
         }  else  {                                                             ### mehrdimensional
            stopifnot(nrow(shw$cov.structure) == ncol(qMatrix))                 ### (Residual-)Varianzen und (Residual-)Korrelationen der lat. Dimensionen
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
         reg  <- shw$regression                                                 ### untere Zeile: Dimensionen analog zu Q matrix umbenennen
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
         wle  <- get.wle( file.path(path, wleFile) )                            ### unten: ins mittel-longformat, um rel zu bestimmen
         res  <- NULL                                                           ### initialisieren
         for ( i in 1:nrow(altN)) { colnames(wle) <- gsub(  paste(".",altN[i,"nr"],"$",sep=""), paste("_", altN[i,"to"],sep="") , colnames(wle))}
         wleL <- reshape2::melt(wle, id.vars = "ID", measure.vars = colnames(wle)[-c(1:2)], na.rm=TRUE)
         foo  <- data.frame ( eatTools::halveString( as.character(wleL[,"variable"]), pattern = "_"), stringsAsFactors=FALSE)
         colnames(foo) <- c("par", "group")                                     ### halveString statt strsplit nehmen, weil es sonst schiefgeht, wenn dimensionsname einen Unterstrich enthaelt; foo muss immer zwei spalten haben, das ist nicht so, wenn man strsplit nimmt und der dimensionsname einen Unterstrich enthaelt
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
         drinI <- match( shw[["item"]][,"item"], colnames(daten))               ### ggf.: welche Items im Datensatz stehen nicht im Showfile (*.shw)?
         drinP <- match(theta[,1], daten[,"ID"])                                ### ggf.: welche Personen im Datensatz stehen nicht im PV-File
         stopifnot(length(which(is.na(drinP))) == 0 , length(which(is.na(drinI))) == 0 )
         q3.res<- sirt::Q3(dat = daten[drinP,drinI], theta = theta[,2], b = shw[["item"]][,"ESTIMATE"], progress = FALSE)
         nObs  <- NULL
         if ( !is.null(q3MinObs) ) {                                            ### untere Zeile: paarweise Anzahl Beobachtungen je Itempaar
              if ( q3MinObs > 1 ) { nObs <- nObsItemPairs ( responseMatrix = daten[,all.Names[["variablen"]]], q3MinType = q3MinType ) }
         }
         matL  <- reshapeQ3 (mat = q3.res$q3.matrix, q3MinObs = q3MinObs, nObs = nObs)
         if( nrow(matL)== 0) { return(NULL)}
         res   <- data.frame ( model = model.name, source = "conquest", var1 = matL[,"Var1"],  var2 = matL[,"Var2"] , type = "fixed",indicator.group = "items",group = paste(names(table(shw1[,"group"])), collapse="_"), par = "q3", derived.par = NA, value = matL[,"value"] , stringsAsFactors = FALSE)
         return(res)}

getConquestDeviance <- function ( path, analysis.name, omitUntil = omitUntil) {
    ### erstmal zusaetzliche Informationen (Anzahl nodes etc.) gewinnen
         cqc  <- scan(file.path ( path, paste0(analysis.name, ".cqc")),what="character",sep="\n",quiet=TRUE)
         such <- c("method", "nodes")
         ret  <- lapply(such, FUN = function ( su ) {
                 indm <- grep(paste0(su, "="), cqc)
                 if ( length(indm)>1) {                                         ### schlechter Hotfix, f_nodes entfernen
                    hf   <- grep("f_nodes", cqc)
                    indm <- setdiff(indm, hf)
                 }
                 if(length(indm) != 1) {
                    cat(paste("Cannot identify '",su,"'from cqc file.\n",sep=""))
                    met <- NULL
                 }  else  {
                    pos1<- nchar(unlist(strsplit(cqc[indm], su))[1])            ### position finden, an der 'method' steht
                    pos2<- which(sapply(1:nchar(cqc[indm]), FUN = function(x){ substr(cqc[indm],x,x) == ","}))
                    pos2<- min(pos2[which(pos2>pos1)])
                    met <- eatTools::removePattern(substr(cqc[indm], pos1+1, pos2-1), paste0(su,"="))
                 }
                 return(met)})                                                  ### Zeit als Differenz von cqc und shw file
         tme  <- file.info ( file.path ( path, paste0(analysis.name, ".shw")))[["mtime"]] - file.info ( file.path ( path, paste0(analysis.name, ".cqc")))[["mtime"]]
         grDevices::pdf(file = file.path ( path, paste0(analysis.name, "_dev.pdf")), width = 10, height = 7.5)
         plotDevianceConquest ( logFile = list ( path=path, analysis.name=analysis.name, ret=ret, tme=tme), omitUntil = omitUntil)
         grDevices::dev.off() }


getConquestResults<- function(path, analysis.name, model.name, qMatrix, all.Names, abs.dif.bound , sig.dif.bound, p.value, deskRes, discrim, omitFit, omitRegr, omitWle, omitPV, daten, Q3=Q3, q3theta=q3theta, q3MinObs =  q3MinObs, q3MinType = q3MinType, omitUntil) {
         allFiles <- list.files(path=path, pattern = analysis.name, recursive = FALSE)
         qL       <- reshape2::melt(qMatrix, id.vars = colnames(qMatrix)[1], variable.name = "dimensionName", na.rm=TRUE)
         qL       <- qL[which(qL[,"value"] != 0 ) , ]
         varName  <- colnames(qMatrix)[1]
         ret      <- NULL                                                       ### Rueckgabeobjekt initialisieren
    ### Sektion 'Konvergenz pruefen' (log)
         logFile  <- paste(analysis.name, "log", sep=".")
         isConv   <- converged ( dir = path, logFile = logFile )
         isPoly   <- length(unique(deskRes[,"Codes"]))>1                        ### war modell polytom? damit das geht, muss es immer deskriptive Statistiken geben, muss also in 'defineModel' obligatorisch sein!
    ### Deviance als pdf plotten
         plotPdf  <- getConquestDeviance(path=path, analysis.name = analysis.name, omitUntil = omitUntil)
    ### Itemparameter auslesen (itn)
         ret      <- rbind(ret, getConquestItn (model.name=model.name, analysis.name=analysis.name, qMatrix=qMatrix, qL=qL, allFiles=allFiles, isPoly=isPoly, path=path))
    ### Descriptives auslesen
         ret      <- rbind(ret, getConquestDesc (model.name=model.name, deskRes = deskRes, qMatrix=qMatrix, qL = qL, isPoly=isPoly))
    ### Diskrimination auslesen
         ret      <- rbind(ret, getConquestDiscrim (model.name=model.name, discrim = discrim, qMatrix=qMatrix, qL = qL))
    ### Itemparameter auslesen (shw): alle folgenden Funktionen werden nur aufgerufen, wenn es ein showfile gibt
         shwFile  <- paste(analysis.name, "shw", sep=".")
         if (!shwFile %in% allFiles) {
             cat("Cannot find Conquest showfile.\n")
         } else {
             shw  <- get.shw( file = file.path(path, shwFile) )                 ### untere Zeile: 'reine' itemparameter auslesen
             if(is.null( dim(shw$cov.structure) )) {from <- NA} else { from <- shw$cov.structure[-ncol(shw$cov.structure),1]}
             altN <- data.frame ( nr = 1:(ncol(qMatrix)-1), pv = paste("dim", 1:(ncol(qMatrix)-1),sep="."), from = from ,  to = colnames(qMatrix)[-1], stringsAsFactors = FALSE)
             shw[["item"]]  <- merge(shw[["item"]], qL[,-match("value", colnames(qL))], by.x = "item", by.y = colnames(qMatrix)[1], all=TRUE)
             shw12<- getConquestShw (model.name=model.name, qMatrix=qMatrix, qL=qL, shw=shw, altN=altN)
             ret  <- rbind(ret, shw12[["shw1"]], shw12[["shw2"]])
             ret  <- rbind(ret, getConquestInfit (model.name=model.name, shw=shw))
             ret  <- rbind(ret, getConquestAdditionalTerms (model.name=model.name, qMatrix=qMatrix, shw=shw, shwFile = shwFile))
    ### reliabilitaeten ergaenzen
             ret  <- rbind(ret, data.frame ( model = model.name, source="conquest", var1=NA, var2=NA,type="tech", indicator.group="persons", group = colnames(qMatrix)[-1], par="eap", derived.par = "rel", value = shw[["reliability"]][,"eap.rel"], stringsAsFactors=FALSE))
    ### Populationsparameter und Regressionsparameter aus Showfile auslesen (shw)
             ret  <- rbind(ret, getConquestPopPar (model.name=model.name, qMatrix=qMatrix, shw=shw))
             ret  <- rbind(ret, getConquestRegPar (model.name=model.name, shw=shw, altN = altN))
    ### Sektion 'Modellindizes auslesen' (shw)
             ret  <- rbind(ret, data.frame ( model = model.name, source = "conquest", var1 = NA, var2 = NA , type = "model", indicator.group = NA, group = NA, par = c("deviance", "Npar"),  derived.par = NA, value = shw$final.deviance , stringsAsFactors = FALSE))
    ### Personenparameter auslesen (wle) ... da hierzu das Objekt 'altN' gebraucht wird, das aus dem shw-file erzeugt wird, geht das Auslesen von WLEs nur, wenn das Auslesen von shw geklappt hat
             wles <- getConquestWles (model.name=model.name, analysis.name=analysis.name, qMatrix=qMatrix, allFiles=allFiles, omitWle = omitWle, altN = altN, path=path)
             ret  <- rbind(ret, wles[["res"]])
             pvs  <- getConquestPVs (model.name=model.name, analysis.name=analysis.name, omitPV = omitPV, altN = altN, path=path, allFiles=allFiles)
             ret  <- rbind(ret, pvs[["res"]])
    ### Q3 erzeugen
             ret  <- rbind(ret, getConquestQ3 (model.name=model.name, shw=shw,Q3=Q3, q3theta=q3theta, omitWle=omitWle, omitPV=omitPV, pv=pvs[["pv"]],wle=wles[["wle"]],daten=daten,all.Names=all.Names, q3MinObs=q3MinObs, q3MinType=q3MinType, shw1 = shw12[["shw1"]]))
         }                                                                      ### schliesst die Bedingung 'shw file vorhanden'
         if(!is.null(ret)) {
             attr(ret, "isConverged") <- isConv
             attr(ret, "available")   <- list ( itn =  paste(analysis.name, "itn", sep=".") %in% allFiles, shw =  paste(analysis.name, "shw", sep=".") %in% allFiles, wle = ( paste(analysis.name, "wle", sep=".") %in% allFiles) & (omitWle == FALSE), pv = ( paste(analysis.name, "pvl", sep=".") %in% allFiles) & (omitPV == FALSE))
         }
         return(ret)}

### Teilfunktionen fuer 'getTamResults()' zum Auslesen der Itemparameter (Schwierigkeiten, p-Werte) auslesen
getTamItempars    <- function(runModelObj, qL, qMatrix, leseAlles) {
         if(leseAlles == FALSE) {return(NULL)}
         if ( is.null(attr(runModelObj, "all.Names")[["DIF.var"]])) {           ### wenn kein DIF: konventionell mergen
              xsis <- merge(data.frame ( item = rownames(runModelObj[["xsi"]]), runModelObj[["xsi"]], stringsAsFactors = FALSE), qL[,-match("value", colnames(qL))],  by.x = "item", by.y = colnames(qMatrix)[1], all = TRUE)
         }  else  {                                                             ### bei DIF anders mergen
              xsis <- mergeDimensionIfDIF (dat = data.frame ( item = rownames(runModelObj[["xsi"]]), runModelObj[["xsi"]], stringsAsFactors = FALSE), qmatLong = qL[,-match("value", colnames(qL))], datMergingVar="item", remove = "toMerge")
         }
         shw1 <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = xsis[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = xsis[,"dimensionName"], par = "est",  derived.par = NA, value = xsis[,"xsi"], stringsAsFactors = FALSE)
         shw2 <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = xsis[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = xsis[,"dimensionName"], par = "est",  derived.par = "se", value = xsis[,"se.xsi"], stringsAsFactors = FALSE)
         if ( !is.null(attr(runModelObj, "all.Names")[["DIF.var"]])) {          ### wenn DIF: items umbenennen
               shw1 <- renameDifParameters (dat=shw1, qmatLong = qL[,-match("value", colnames(qL))])
               shw2 <- renameDifParameters (dat=shw2, qmatLong = qL[,-match("value", colnames(qL))])
         }
         toOff<- shw2[ which(shw2[,"value"] == 0 ), "var1"]                     ### verankerte Parameter identifizieren
         if(length(toOff)>0) {
            shw1[match(toOff, shw1[,"var1"]), "par"] <- "offset"
            shw2  <- shw2[-which(shw2[,"value"] == 0 ),] }                      ### entferne Zeilen aus shw2, die in der "value"-Spalte NA haben, danach: p-Werte einfuegen
         return(list ( shw1=shw1, shw2=shw2))}

getTamDescriptives    <- function(runModelObj, qL, qMatrix, leseAlles) {
         if(leseAlles == FALSE || is.null ( attr(runModelObj, "deskRes") )) {return(NULL)}
         deskR<- merge(attr(runModelObj, "deskRes"), qL[,-match("value", colnames(qL))],  by.x = "item.name", by.y = colnames(qMatrix)[1], all = TRUE)
         shw3 <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = as.character(deskR[,"item.name"]), var2 = NA , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "itemP",  derived.par = NA, value = deskR[,"item.p"], stringsAsFactors = FALSE)
         shw4 <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = as.character(deskR[,"item.name"]), var2 = NA , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "Nvalid",  derived.par = NA, value = deskR[,"valid"], stringsAsFactors = FALSE)
    ### Achtung! wenn in dem 'deskRes'-Objekt noch mehr p-Werte (schulformspezifische p-Werte drinstehen, werden die jetzt auch in die Ergebnisstruktur eingetragen)
         cols <- setdiff ( colnames(deskR)[grep("^item.p", colnames(deskR))], "item.p")
         if ( length ( cols ) > 0 ) {
              colsR <- data.frame ( original = cols, reduziert = eatTools::removePattern ( string = cols, pattern = "item.p.") , stringsAsFactors = FALSE)
              shw31 <- do.call("rbind", apply ( colsR, MARGIN = 1, FUN = function ( zeile ) { data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = as.character(deskR[,"item.name"]), var2 = zeile[["reduziert"]] , type = "fixed", indicator.group = "items", group = deskR[,"dimensionName"], par = "itemP",  derived.par = NA, value = deskR[,zeile[["original"]]], stringsAsFactors = FALSE) }))
              return(rbind(shw3, shw4, shw31))
         }  else  {
              return(rbind(shw3, shw4))
         } }

getTamDiscrim    <- function(runModelObj, qL, qMatrix, leseAlles) {
         if(leseAlles == FALSE || is.null ( attr(runModelObj, "discrim") )) {return(NULL)}
         discR<- merge( attr(runModelObj, "discrim") , qL[,-match("value", colnames(qL))],  by.x = "item.name", by.y = colnames(qMatrix)[1], all = TRUE)
         shw5 <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = discR[,"item.name"], var2 = NA , type = "fixed", indicator.group = "items", group = discR[,"dimensionName"], par = "itemDiscrim",  derived.par = NA, value = discR[,"item.diskrim"], stringsAsFactors = FALSE)
         return(shw5)}

getTam2plDiscrim <- function(runModelObj, qMatrix, leseAlles, regr, omitRegr) {
         if(leseAlles == FALSE || !attr(runModelObj, "irtmodel") %in% c("2PL", "2PL.groups", "GPCM", "3PL") ) {return(NULL)}
         shw6 <- do.call("rbind", lapply (  1 : length ( colnames( qMatrix ) [-1] ) , FUN = function ( dims ) {
                 if ( isFALSE(omitRegr) ) {                                     ### wenn omitRegr == FALSE, werden die Diskriminationsparameter aus diesem Objekt,
                      obj <- regr[["B"]]                                        ### ansonsten aus dem direkten TAM-Rueckgabeobjekt ausgelesen
                 } else {                                                       ### wenn omitRegr == FALSE, kommen die Diskriminationen als data.frame,
                      obj <- as.data.frame ( runModelObj[["B"]])                ### andernfalls als array, muss also umgewandelt werden
                      colnames(obj) <- paste0("B.", gsub("Dim0", "Dim", colnames(obj)))
                      obj[,"item"]  <- rownames(obj)                            ### Spalten rauswerfen, in denen ausschliesslich nullen stehen
                      isNull        <- which(sapply(obj, FUN = function ( x ) { all(x==0)})==TRUE)
                      if (length (isNull)>0) {
                          obj <- obj[,-isNull]
                      }
                 }
                 cols  <- grep(paste0(".Dim",dims,"$" ), colnames(obj), value=TRUE)
                 tamMat<- obj[,c("item",cols)]
                 weg   <- which(tamMat[,2] == 0)
                 if(length(weg)>0) {tamMat <- tamMat[-weg,]}
                 shw6D <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = tamMat[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = colnames(qMatrix)[dims+1], par = "estSlope",  derived.par = NA, value = tamMat[,2], stringsAsFactors = FALSE)
                 if (ncol(tamMat) == 3 ) {
                     shw6se<- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = tamMat[,"item"], var2 = NA , type = "fixed", indicator.group = "items", group = colnames(qMatrix)[dims+1], par = "estSlope",  derived.par = "se", value = tamMat[,3], stringsAsFactors = FALSE)
                 }  else  {
                     shw6se<- NULL
                 }
                 return(rbind(shw6D, shw6se)) }))
         return(shw6)}

### Hilfsfunktion fuer getTamInfit() und andere
mergeDimensionIfDIF <- function(dat, qmatLong, datMergingVar, remove) {
         dat[,datMergingVar]    <- as.character(dat[,datMergingVar])
         dat[,"toMerge"] <- eatTools::halveString(dat[,datMergingVar], ":", first=TRUE)[,1]
         dat  <- merge(dat, qmatLong[,c("item", "dimensionName")], by.x = "toMerge", by.y = "item", all.x = TRUE)
         dat  <- dat[,-match(remove, colnames(dat))]
         return(dat)}

### Hilfsfunktion fuer getTamInfit() und andere
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


getTamInfit    <- function(runModelObj, qL, qMatrix, leseAlles, omitFit) {
         if(leseAlles == FALSE || omitFit == TRUE ) {return(NULL)}
         infit<- tam.fit(runModelObj, progress=FALSE)                           ### Achtung: wenn DIF-Analyse, dann misslingt untere Zeile: Workarond!
         fits <- merge(infit[["itemfit"]], qL[,-match("value", colnames(qL))],  by.x = "parameter", by.y = colnames(qMatrix)[1], all = TRUE)
         if ( is.null(attr(runModelObj, "all.Names")[["DIF.var"]])) {           ### wenn kein DIF: mergen
              ret  <- rbind(data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = fits[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = fits[,"dimensionName"], par = "est",  derived.par = "infit", value = fits[,"Infit"], stringsAsFactors = FALSE),
                            data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = fits[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = fits[,"dimensionName"], par = "est",  derived.par = "outfit", value = fits[,"Outfit"], stringsAsFactors = FALSE))
         }  else  {                                                             ### wenn DIF: workaround ... DIF-Parameter umbenennen, so dass es konsistent zu "getConquestResults" ist
              ret  <- rbind(data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = infit$itemfit[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = NA, par = "est",  derived.par = "infit", value = infit$itemfit[,"Infit"], stringsAsFactors = FALSE),
                            data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = infit$itemfit[,"parameter"], var2 = NA , type = "fixed", indicator.group = "items", group = NA, par = "est",  derived.par = "outfit", value = infit$itemfit[,"Outfit"], stringsAsFactors = FALSE) )
              ret  <- mergeDimensionIfDIF(dat=ret, qmatLong=qL[,-match("value", colnames(qL))], datMergingVar="var1", remove = c("group", "toMerge"))
              colnames(ret) <- car::recode(colnames(ret), "'dimensionName'='group'")
              ret  <- renameDifParameters(dat=ret, qmatLong=qL[,-match("value", colnames(qL))])
         }
         return(ret)}

getTamPopPar    <- function(runModelObj, qMatrix, leseAlles) {
         if(leseAlles == FALSE ) {return(NULL)}
         if(ncol(qMatrix) == 2) {                                               ### eindimensionaler Fall
            ret  <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = colnames(qMatrix)[2], var2 = NA , type = "distrpar", indicator.group = NA, group = "persons", par = "var",  derived.par = NA, value = runModelObj[["variance"]][1,1] , stringsAsFactors = FALSE)
         }  else  {                                                             ### mehrdimensional: (Residual-)Varianzen und (Residual-)Korrelationen der lat. Dimensionen
            cov1 <- runModelObj[["variance"]]
            colnames(cov1) <- colnames(qMatrix)[-1]
            rownames(cov1) <- colnames(qMatrix)[-1]
            cor1 <- cov2cor(cov1)
            for (ii in 1:nrow(cor1))   {                                        ### loesche alles oberhalb der Hauptdiagonalen
                 cor1[ii,ii:ncol(cor1)] <- NA}
            cor1 <- reshape2::melt(cor1, measure.vars = colnames(cor1), na.rm=TRUE)
            vars <- Matrix::diag(cov1)
            ret  <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = c(names(vars),as.character(cor1[,"Var1"])) , var2 = c(rep(NA, length(vars)), as.character(cor1[,"Var2"])) , type = "random", indicator.group = NA, group = "persons", par = c(rep("var",length(vars)), rep("correlation", nrow(cor1))) ,  derived.par = NA, value = c(unlist(vars), cor1[,"value"]), stringsAsFactors = FALSE)
         }
         return(ret)}

getTamRegPar    <- function(runModelObj, leseAlles, qMatrix, omitRegr, regr) {
         if(leseAlles == FALSE || omitRegr == TRUE ) {return(NULL)}
         if( !isTRUE(all.equal ( dim(runModelObj$beta) , c(1,1))))  {           ### wird nur gemacht, wenns auch Regressionsparameter gibt
             regr <- data.frame ( reg.var = rownames(regr$beta), regr$beta, stringsAsFactors = FALSE)
             regr <- reshape2::melt(regr, id.vars = "reg.var", na.rm=TRUE)
             regr2<- data.frame ( par = "est", derived.par = car::recode(unlist(lapply(strsplit(as.character(regr[,"variable"]),"\\."), FUN = function ( l ) {l[1]})), "'se'='se'; else=NA"), group = colnames(qMatrix)[as.numeric(eatTools::removePattern( string = unlist(lapply(strsplit(as.character(regr[,"variable"]),"\\."), FUN = function ( l ) {l[2]})), pattern = "Dim")) + 1], regr, stringsAsFactors = FALSE)
             regr3<- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = regr2[,"reg.var"], var2 = NA , type = "regcoef", indicator.group = NA, group = regr2[,"group"], par = regr2[,"par"],  derived.par = regr2[,"derived.par"], value = regr2[,"value"] , stringsAsFactors = FALSE)
         }  else {
             return(NULL)
         }
         return(regr3)}

getTamModInd    <- function(runModelObj, leseAlles) {
         if(leseAlles == FALSE ) {return(NULL)}
         return(data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = NA, var2 = NA , type = "model", indicator.group = NA, group = NA, par = c("deviance", "Npar", "AIC", "BIC"), derived.par = NA, value = unlist(runModelObj[["ic"]][c("deviance", "Npars", "AIC", "BIC")]), stringsAsFactors = FALSE))}

getTamWles    <- function(runModelObj, qMatrix, leseAlles, omitWle) {
         if(leseAlles == FALSE || omitWle == TRUE ) {return(NULL)}
         txt  <- capture.output(wle  <- tam.wle(runModelObj, progress = FALSE)) ### Achtung: im eindimensionalen Fall enthalten die Spaltennamen keine Benennung der Dimension
         eind1<- ncol(wle) == 7                                                 ### ist das eindimensional?
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
         rel  <- reshape2::melt(as.data.frame ( wle)[1,weg1], na.rm = TRUE)     ### das Auslesen der Dimensionsnamen fuer den ein- und zweidimensionalen Fall (Objekt 'trans')
         if ( "variable" %in% colnames(rel)) {                                  ### ist noch nicht wirklich schoen, muesste gegebenenfalls (wenn mal Zeit ist) elegantisiert werden
               rel[,"original"] <- unlist(lapply(strsplit(as.character(rel[,"variable"]),"\\."), FUN = function (l) {l[length(l)]}))
               rel  <- merge(rel, trans, by="original", all=TRUE)
         }  else  {
               rel[,"uebersetzt"] <- colnames(qMatrix)[-1]
         }
         res  <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = c(wleL[,"pid"],rep(NA,nrow(rel))), var2 = NA , type = c(rep("indicator",nrow(wleL)), rep("tech", nrow(rel))), indicator.group = "persons", group = c(wleL[,"group"],rel[,"uebersetzt"]), par = c(wleL[,"par"],rep("wle", nrow(rel))),  derived.par = c(wleL[,"derived.par"],rep("rel",nrow(rel))), value = c(wleL[,"value"] ,rel[,"value"]), stringsAsFactors = FALSE)
         return(res)}

getTamPVs <- function ( runModelObj, qMatrix, leseAlles, omitPV, pvMethod, tam.pv.arguments) {
         if(omitPV == TRUE ) {return(NULL)}
         for ( i in names( tam.pv.arguments )) { assign(i, tam.pv.arguments[[i]]) }
         if(leseAlles == TRUE ) {
            if ( pvMethod == "regular" ) {                                      ### konventionell
                 do   <- paste ( "tam.pv ( ", paste(names(formals(tam.pv)), car::recode ( names(formals(tam.pv)), "'tamobj'='runModelObj'"), sep =" = ", collapse = ", "), ")",sep="")
            } else {                                                            ### bayesianisch, aber mit vorhandedem 'tam.mml'-Objekt
                 if ( is.null ( attr(runModelObj, "Y") ) ) {
                      warning("Conditioning model was not defined ('Y' is NULL).")
                      Y1 <- NULL
                 } else {
                      Y1 <- data.frame ( intercpt = 1, attr(runModelObj, "Y"))
                 }
                 do   <- paste ( "tam.pv.mcmc ( ", paste(names(formals(tam.pv.mcmc)), car::recode ( names(formals(tam.pv.mcmc)), "'tamobj'='runModelObj'; 'Y'='Y1'"), sep =" = ", collapse = ", "), ")",sep="")
            }
         }  else  {                                                             ### PV-Ziehung ohne vorhandenes 'tam.mml'-Objekt
            class(runModelObj) <- "list"                                        ### gibt sonst eine komische Warnung bei TAM
            stopifnot ( pvMethod == "bayesian")
            do   <- paste ( "tam.pv.mcmc ( ", paste(names(formals(tam.pv.mcmc)), car::recode ( names(formals(tam.pv.mcmc)), "'tamobj'='runModelObj'; 'Y'='runModelObj[[\"Y\"]]'; 'nplausible'='attr(runModelObj, \"n.plausible\")'"), sep =" = ", collapse = ", "), ")",sep="")
         }
         pv   <- eval(parse(text=do))
         pvL  <- reshape2::melt(pv$pv, id.vars = "pid", na.rm=TRUE)
         pvL[,"PV.Nr"] <- as.numeric(eatTools::removePattern(string = unlist(lapply(strsplit(as.character(pvL[,"variable"]),"\\."), FUN = function (l) {l[1]})), pattern = "PV"))
         pvL[,"group"] <- colnames(qMatrix)[as.numeric(eatTools::removePattern(string = unlist(lapply(strsplit(as.character(pvL[,"variable"]),"\\."), FUN = function (l) {l[2]})), pattern = "Dim"))+1]
         res <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = pvL[,"pid"], var2 = NA , type = "indicator", indicator.group = "persons", group = pvL[,"group"], par = "pv",  derived.par = paste("pv", pvL[,"PV.Nr"],sep=""), value = pvL[,"value"] , stringsAsFactors = FALSE)
         return(res)}

getTamEAPs <- function ( runModelObj, qMatrix, leseAlles = leseAlles) {
         if(leseAlles == FALSE ) {return(NULL)}
         eaps <- runModelObj[["person"]]                                        ### Achtung: im eindimensionalen Fall enthalten die Spaltennamen keine Benennung der Dimension
         eind1<- ncol(eaps) == 7                                                ### (uneinheitlich zu pvs, wo es immer eine Benennung gibt.)
         if(eind1 == TRUE) {                                                    ### Im eindimensionalen Fall muss Benennung ergaenzt werden
            cols <- grep("EAP$", colnames(eaps))                                ### zur Sicherheit werden hier zwei Indikatoren fuer Eindimensionalitaet genutzt. Fehlermeldung bei Widerspruch
            stopifnot(length(cols) == 2)                                        ### ggf. muss diese Passage nach Release neuerer TAM-Versionen korrigiert werden
            colnames(eaps)[cols] <- paste(colnames(eaps)[cols], ".Dim1", sep="")
         }
         eaps <- reshape2::melt(eaps, id.vars = "pid", measure.vars = grep("EAP", colnames(eaps)), na.rm=TRUE)
         eaps[,"group"] <- colnames(qMatrix)[as.numeric(eatTools::removePattern ( string = eatTools::halveString(string = as.character(eaps[,"variable"]), pattern = "\\.", first = FALSE)[,"X2"], pattern = "Dim"))+1]
         eaps[,"par"]   <- "est"
         eaps[grep("^SD.",as.character(eaps[,"variable"])),"par"]   <- "se"
         res  <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = eaps[,"pid"], var2 = NA , type = "indicator", indicator.group = "persons", group = eaps[,"group"], par = "eap", derived.par = eaps[,"par"], value = eaps[,"value"] , stringsAsFactors = FALSE)
         if (ncol(qMatrix)>2) { grp  <- names(runModelObj[["EAP.rel"]]) } else { stopifnot(length(unique(eaps[,"group"])) == 1); grp <- unique(eaps[,"group"]) }
         res  <- rbind(res, data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = NA, var2 = NA , type = "tech", indicator.group = "persons", group = grp, par = "eap", derived.par = "rel", value = runModelObj[["EAP.rel"]] , stringsAsFactors = FALSE) )
         return(res)}


getTamQ3 <- function(runModelObj, leseAlles, shw1, Q3, q3MinObs, q3MinType){
         if(leseAlles == FALSE || Q3 == FALSE) {return(NULL)}
         nObs <- NULL
         if ( !is.null(q3MinObs) ) {                                            ### untere Zeile: paarweise Anzahl Beobachtungen je Itempaar
              if ( q3MinObs > 1 ) {
                   nObs <- nObsItemPairs ( responseMatrix = runModelObj[["resp"]], q3MinType = q3MinType )
              }
              mat  <- tam.modelfit ( tamobj = runModelObj, progress = FALSE )
              matL <- reshapeQ3 (mat = mat$Q3.matr, q3MinObs = q3MinObs, nObs = nObs)
              if( nrow(matL)>0) {
                  res  <- data.frame ( model = attr(runModelObj, "analysis.name"), source = "tam", var1 = matL[,"Var1"], var2 = matL[,"Var2"] , type = "fixed",indicator.group = "items", group = paste(names(table(shw1[,"group"])), collapse="_"), par = "q3", derived.par = NA, value = matL[,"value"] , stringsAsFactors = FALSE)
              }  else  {
                  res  <- NULL
              }
         }
         return(res)}


getTamResults     <- function(runModelObj, omitFit, omitRegr, omitWle, omitPV, nplausible , ntheta , normal.approx, samp.regr, theta.model, np.adj, Q3=Q3, q3MinObs =  q3MinObs, q3MinType = q3MinType,
                     pvMethod , group, beta_groups , level , n.iter , n.burnin, adj_MH , adj_change_MH , refresh_MH, accrate_bound_MH,	sample_integers, theta_init, print_iter , verbose, calc_ic) {
         qMatrix<- attr(runModelObj, "qMatrix")
         qL     <- reshape2::melt(qMatrix, id.vars = colnames(qMatrix)[1], variable.name = "dimensionName", na.rm=TRUE)
         qL     <- qL[which(qL[,"value"] != 0 ) , ]
         varName<- colnames(qMatrix)[1]                                         ### untere Zeile: Standardfehler auslesen, falls vorhanden
         if( omitRegr == FALSE && !inherits(runModelObj, "tamBayes")) {
             txt <- capture.output ( regr <- tam.se(runModelObj))               ### Namen der Regressoren stehen nicht im tam-Output 'reg' drin, nur Ziffern
             stopifnot ( nrow(regr$beta) == ncol(attr(runModelObj, "Y") )+1)    ### die Namen muessen daher jetzt wieder aus den Spaltennamen der Y-Matrix rekonstruiert werden
             rownames(regr$beta) <- c("(Intercept)", colnames(attr(runModelObj, "Y")))
          } else {
             regr <- NULL
         }
    ### wenn PVs bayesianisch gezogen werden sollen ohne dass 'tam.mml' aufgerufen wurde, muessen alle Schritte bis zur PV-Ziehung nun uebersprungen werden
         if ( !inherits(runModelObj, "tamBayes") ) {leseAlles <- TRUE} else {leseAlles <- FALSE}
         ret    <- NULL                                                       ### Rueckgabeobjekt initialisieren, und untere Zeile: Itemparameter auslesen
         resItem<- getTamItempars(runModelObj=runModelObj, qL=qL, qMatrix=qMatrix, leseAlles = leseAlles)
         ret    <- rbind(ret, resItem[["shw1"]], resItem[["shw2"]])
    ### deskriptive Werte auslesen
         ret    <- rbind(ret, getTamDescriptives(runModelObj=runModelObj, qL=qL, qMatrix=qMatrix, leseAlles = leseAlles))
    ### Diskriminationswerte auslesen
         ret    <- rbind(ret, getTamDiscrim(runModelObj=runModelObj, qL=qL, qMatrix = qMatrix, leseAlles = leseAlles))
    ### 2pl-Diskriminationsparameter auslesen
         ret    <- rbind(ret, getTam2plDiscrim(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, regr = regr, omitRegr=omitRegr))
    ### Infit auslesen
         ret    <- rbind(ret, getTamInfit(runModelObj=runModelObj, qL=qL, qMatrix = qMatrix, leseAlles = leseAlles, omitFit = omitFit))
    ### Populationsparameter auslesen
         ret    <- rbind(ret, getTamPopPar(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles))
    ### Regressionsparameter auslesen
         ret    <- rbind(ret, getTamRegPar(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, omitRegr = omitRegr, regr=regr))
    ### Modellindizes auslesen
         ret    <- rbind(ret, getTamModInd(runModelObj=runModelObj, leseAlles = leseAlles))
    ### Personenparameter auslesen (WLEs)
         ret    <- rbind(ret, getTamWles(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, omitWle = omitWle))
    ### PVs auslesen
         tamArg <- as.list(match.call(definition = getTamResults))
         weg    <- which(names(tamArg) %in% c("tamobj", "Y", "runModelObj", "qMatrix", "leseAlles", "omitPV", "pvMethod", "omitFit", "omitRegr", "omitWle"))
         if ( length(weg)>0) {tamArg <- tamArg[-weg]}
         # tamArg <- lapply ( tamArg[2:length(tamArg)], eval )                  ### das klappt irgendwie nicht, also in den folgenden
         tamarg <- list()                                                       ### vier Zeilen auf die komplizierte Variante
         for ( i in 2:length(tamArg)) {
              tamarg[[names(tamArg)[i]]] <- eval(tamArg[[i]])
         }
         retPVs <- getTamPVs ( runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles, omitPV = omitPV, pvMethod = pvMethod, tam.pv.arguments = tamarg)
         ret    <- rbind(ret, retPVs)
    ### EAPs auslesen
         ret    <- rbind(ret, getTamEAPs(runModelObj=runModelObj, qMatrix=qMatrix, leseAlles = leseAlles))
    ### Q3 auslesen
         ret    <- rbind(ret, getTamQ3(runModelObj=runModelObj, leseAlles = leseAlles, shw1 = resItem[["shw1"]], Q3=Q3, q3MinObs=q3MinObs, q3MinType=q3MinType))
         return(ret)}


### Hilfsfunktion zur Vereinheitlichung der Q3-Matrix
reshapeQ3 <- function ( mat, q3MinObs, nObs ) {
             for (ii in 1:(nrow(mat)-1)) { mat[ii,ii:ncol(mat)] <- NA}          ### entferne alles oberhalb der Hauptdiagonale
             matL <- reshape2::melt ( mat , na.rm = TRUE)                                 ### das entfernt alle doppelten Eintraege
             if ( !is.null(nObs)) {                                             ### dass hier soll nur passieren, wenn Eintraege aus der Q3 Matrix ggf. entfernt werden
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

### Extraktorfunktionen
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
              re <- by(re, INDICES = re[,"model"], FUN = function (m) {
                    mw <- reshape2::dcast(m, var1~group+derived.par, value.var="value")
                    colnames(mw) <- gsub("_NA$", "_est", colnames(mw))          ### 'do' = domain names
                    do <- unique(eatTools::halveString(colnames(mw)[-1], "_", first=FALSE)[,1])
                    for ( u in do) {                                            ### 'xp' statt 'p', damit die Spaltensortierung einfacher klappt (wird am Ende wieder richtig zurueckbenannt)
                         mw[,paste0(u, "_xp")]   <- 2*(1-pnorm(abs(mw[,paste0(u,"_est")] / mw[,paste0(u,"_se")])))
                         mw[,paste0(u, "_ysig")] <- eatTools::num.to.cat(mw[,paste0(u, "_xp")], cut.points = c(0.001, 0.01, 0.05), cat.values = c("***", "**", "**", ""))
                    }
                    mw <- mw[,c(colnames(mw)[1], sort(colnames(mw)[-1]))]
                    colnames(mw) <- gsub("_ysig$", "_sig", gsub("_xp$", "_p", colnames(mw)))
                    if(length(do) ==1) {                                        ### Namen des Kompetenzbereichs aus Spaltennamen entfernen, wenn es nur einen gibt
                       colnames(mw) <- gsub(paste0("^", do, "_"), "",colnames(mw))
                    }
                    colnames(mw)[1] <- "parameter"
                    if(!is.null(digits)) {mw <- eatTools::roundDF(mw, digits =digits)}
                    return(mw)})
              return(re)
          }}

pvFromRes  <- function ( resultsObj, toWideFormat = TRUE, idVarName = NULL, verbose=TRUE) {
          pvRow<- intersect( which(resultsObj[,"par"] == "pv"),which(resultsObj[,"indicator.group"] == "persons"))
          if ( length ( pvRow ) == 0 ) {
               warning("'resultsObj' does not contain any pv values.")
               return ( NULL )
          }  else  {
             sel  <- resultsObj[pvRow, ]                                        ### Hotfix: ID namen identifizieren
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

### ID identifizieren (zur Kompatibilitaet mit aelteren Paketversionen
getIdVarName <- function ( id, idVarName, verbose=TRUE) {
          if (length( id ) == 0 ) {
              if ( is.null(idVarName)) { new <- "idstud"} else { new <- idVarName}
              if(verbose){warning(paste0("Cannot identify student identifier variable (possibly because 'resultsObj' was created by an older version of 'eatModel'). student id variable will be defaulted to '",new,"'."))}
              id <- new
          }
          return(id)}

itemFromRes<- function ( resultsObj ) {                                         ### Funktion wird so oft ausgefuehrt, wie es Modelle gibt
     ### hier muss "rbind.fill" genommen werden, denn 1pl und 2pl Modelle unterscheiden sich in den Spalten (bei 2pl gibt es zusaetzliche Diskriminationsspalten)
          res <- do.call(plyr::rbind.fill, by ( data = resultsObj, INDICES = resultsObj[,"model"], FUN = function ( mod ) {
                 sel  <- mod[intersect( which(mod[,"par"] %in% c("est", "estSlope", "Nvalid", "itemP", "ptBis", "itemDiscrim", "offset")),which(mod[,"indicator.group"] == "items")),]
                 if (nrow(sel)==0) {
                     return(NULL)
                 }  else  {
     ### gibt es DIF? wenn ja, wird das separat ausgelesen
                     isDif<- intersect(which(mod[,"type"] == "tech"), which(mod[,"par"] == "DIF.var"))
                     if ( length( isDif ) > 0 ) {                               ### untere Zeile: alle Variablen auslesen (quasi ein Hilfsobjekt)
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
                           sel      <- sel[which ( sel[,"par"] != "ptBis" ) , ] ### Hotfix: wenn DIF ausgegeben, wird keine ptBis berechnet
                           selDIF   <- do.call("rbind", by(selForDif, INDICES = selForDif[,"group"], FUN = function ( gr ) {
                                       res  <- reshape2::dcast ( gr , model+var1~par+derived.par, value.var = "value")
                                       mat  <- lapply( vars, FUN = function ( v ) { grep(paste0("_",v,"_"), res[,"var1"])})
                                       stopifnot (  all ( sapply(mat, length) == 1) )
                                       res[unlist(mat),"item"]  <- vars
                                       colnames(res) <- car::recode ( colnames(res) , "'est_infit'='infitDif'; 'est_se'='seDif'; 'est_NA'='estDif'")
                                       res[,"absDif"]<- abs ( res[,"estDif"]  * 2 )
                                       pval <- intersect(intersect(which(mod[,"type"] == "tech"), which(mod[,"par"] == "dif")), which(mod[,"derived.par"] == "p.value"))
                                       stopifnot (length(pval) == 1)
                                       pval <- mod[pval, "value"]               ### untere Zeile: adb = 'abs.dif.bound'; sdb = 'sig.dif.bound'
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
     ### Implementiere Formel nach Lord (1980) und ETS-Klassifikation von DIF; siehe Funktion equating.rasch aus 'eatRest'
                                                              ets   <- "A"
                                                              ets1  <- d[["absDif"]] > 0.43 & d[["absDif"]] < 0.64
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
     ### erstmal ohne schulformspezifische p-Werte (die kommen spaeter dazu)
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

q3FromRes<- function ( resultsObj, out = c("wide", "long" )) {
       out   <- match.arg(arg = out, choices = c("wide", "long" ))
       selM  <- by(data = resultsObj, INDICES = resultsObj[,"model"], FUN = function ( mr ) {
                sel  <- mr[which(mr[,"par"] == "q3"),]
                if ( nrow(sel)>0) {
                     if ( out == "wise") {
                          sel  <- reshape2::dcast( sel, var1~var2, value.var = "value")
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

### 19. Oktober 2011: Funktion wird "generisch". Wenn 'datei' ein string ist, werden PVs aus Conquest eingelesen,
### wenn Datei vom Typ 'mer' ist, werden PVs aus lmer generiert.
### 22. August 2014: Funktion hoert auf, generisch zu sein
get.plausible <- function(file, quiet = FALSE, forConquestResults = FALSE)  {   ### hier beginnt Einlesen fuer Plausible Values aus Conquest
                 input           <- scan(file,what="character",sep="\n",quiet=TRUE)
                 input           <- strsplit(eatTools::crop(gsub("-"," -",input) ) ," +")
    ### gibt die maximale Spaltenanzahl
                 n.spalten       <- max ( sapply(input,FUN=function(ii){ length(ii) }) )
                 input           <- data.frame( matrix( t( sapply(input,FUN=function(ii){ ii[1:n.spalten] }) ),length(input), byrow = FALSE), stringsAsFactors = FALSE)
                 pv.pro.person   <- sum (input[-1,1]==1:(nrow(input)-1) )       ### Problem: wieviele PVs gibt es pro Person? Kann nicht suchen, ob erste Ziffer ganzzahlig, denn das kommt manchmal auch bei Zeile 7/8 vor, wenn entsprechende Werte = 0.0000
                 n.person        <- nrow(input)/(pv.pro.person+3)               ### Anzahl an PVs pro Person wird bestimmt anhand der uebereinstimmung der ersten Spalte mit aufsteigenden 1,2,3,4...
                 weg             <- c(1, as.numeric( sapply(1:n.person,FUN=function(ii){((pv.pro.person+3)*ii-1):((pv.pro.person+3)*ii+1)}) ) )
                 cases           <- input[(1:n.person)*(pv.pro.person+3)-(pv.pro.person+2),1:2]
                 input.sel       <- input[-weg,]
                 n.dim <- dim(input.sel)[2]-1                                   ### Anzahl der Dimensionen
                 if(quiet == FALSE) {cat(paste(n.person,"persons and",n.dim,"dimensions(s) found.\n"))
                               cat(paste(pv.pro.person,"plausible values were drawn for each person on each dimension.\n"))}
                 ID              <- input[  (pv.pro.person + 3) *  (1:n.person) - (pv.pro.person + 2) ,2]
                 colnames(input.sel) <- c("PV.Nr", paste("dim.",1:(ncol(input.sel)-1),sep=""))
                 input.sel[,1]   <- gsub( " ", "0", formatC(input.sel[,1],width = max(nchar(input.sel[,1]))))
                 input.sel$ID    <- rep(ID, each = pv.pro.person)
                 is.na.ID        <- FALSE
                 if(is.na(input.sel$ID[1])) {                                   ### wenn keine ID im PV-File, wird hier eine erzeugt (Fall-Nr), da sonst reshapen misslingt
                    is.na.ID        <- TRUE                                     ### Die ID wird spaeter wieder geloescht. Um das machen zu koennen, wird Indikatorvariable erzeugt, die sagt, ob ID fehlend war.
                    input.sel$ID    <- rep( 1: n.person, each = pv.pro.person)
                 }
                 input.melt      <- reshape2::melt(input.sel, id.vars = c("ID", "PV.Nr") , stringsAsFactors = FALSE)
                 input.melt[,"value"] <- as.numeric(input.melt[,"value"])
                 input.wide      <- data.frame( case = gsub(" ", "0",formatC(as.character(1:n.person),width = nchar(n.person))) , reshape2::dcast(input.melt, ... ~ variable + PV.Nr) , stringsAsFactors = FALSE)
                 colnames(input.wide)[-c(1:2)] <- paste("pv.", paste( rep(1:pv.pro.person,n.dim), rep(1:n.dim, each = pv.pro.person), sep = "."), sep = "")
                 weg.eap         <- (1:n.person)*(pv.pro.person+3) - (pv.pro.person+2)
                 input.eap    <- input[setdiff(weg,weg.eap),]                   ### nimm EAPs und deren Standardfehler und haenge sie an Datensatz - all rows that have not been used before
                 input.eap    <- na.omit(input.eap[,-ncol(input.eap),drop=FALSE])## find EAPs and posterior standard deviations
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

### umgeschrieben fuer 2pl (mit fixiertem slope)
get.wle <- function(file)      {
            input <- eatTools::crop(scan(file, what = "character", sep = "\n", quiet = TRUE))
            input <- strsplit(input," +")
            n.spalten <- max ( sapply(input,FUN=function(ii){ length(ii) }) )   ### Untere Zeile gibt die maximale Spaltenanzahl:
            n.wle <- floor((n.spalten-1) / 4)                                   ### Dies minus eins und dann geteilt durch 4 ergibt Anzahl an WLEs (mehr oder weniger)
            input <- suppressWarnings(eatTools::asNumericIfPossible(data.frame( matrix( t( sapply(input,FUN=function(ii){ ii[1:n.spalten] }) ),length(input),byrow = FALSE), stringsAsFactors = FALSE), force.string = FALSE))
            valid <- na.omit(input)
            cat(paste("Found valid WLEs of ", nrow(valid)," person(s) for ", n.wle, " dimension(s).\n",sep=""))
            if (nrow(valid) != nrow(input)) { cat(paste("    ",nrow(input)-nrow(valid)," persons with missings on at least one latent dimension.\n",sep="")) }
            namen1<- c(rep ( x = c("n.solved", "n.total"), times = n.wle), rep ( x = c("wle", "std.wle"), times = n.wle))
            namen2<- rep(rep ( paste(".", 1:n.wle, sep=""), each = 2),2)        ### untere Zeile: wenn es keine Spalte 'case' gibt, wird die erste Spalte mit 'ID' benannt
            colnames(valid)[(ncol(valid)-length(namen2)):1] <- c("ID","case")[1:(ncol(valid)-length(namen2))]
            colnames(valid)[(ncol(valid)-length(namen2)+1):ncol(valid)] <- paste(namen1,namen2,sep="")
            return(valid)}

### liest Conquest-Outputfiles (*.shw) als R-Objekte ein (siehe auch P:\Aufgabenentwicklung\Grundschule\Daten\Misc\R-Routinen\R2Conquest.R)
### "dif.term" definiert dabei die Variable, nach der ggf. DIF-Analysen ausgelesen werden sollen, z.B. "item*sex". Wird "dif.term" nicht
### spezifiziert, werden keine DIF-Daten ausgegeben.
get.shw <- function(file, dif.term, split.dif = TRUE, abs.dif.bound = 0.6, sig.dif.bound = 0.3, p.value = 0.9) {
            all.output <- list();   all.terms <- NULL                           ### "dif.term" muss nur angegeben werden, wenn DIF-Analysen geschehen sollen.
            input.all <- scan(file,what="character",sep="\n",quiet=TRUE)        ### ginge auch mit:   input <- readLines(file)
            rowToFind <- c("Final Deviance","Total number of estimated parameters")
            rowToFind <- sapply(rowToFind, FUN = function(ii) {                 ### Find the rows indicated in "rowToFind"
                         row.ii <- grep(ii,input.all)                           ### get the parameter of desired rows
                         stopifnot(length(row.ii) == 1)
                         row.ii <- as.numeric(unlist(lapply (strsplit(input.all[row.ii], " +"), FUN=function(ll) {ll[length(ll)]}) ))
                         return(row.ii)})
            ind <- grep("TERM",input.all)                                       ### Wieviele Tabellen gibt es einzulesen?
            grenzen <- grep("An asterisk",input.all)
            if(length(ind)==0) {stop(paste("No TERM-statement found in file ",file,".\n",sep=""))}
            for (i in 1:length(ind)) {
                 term <- input.all[ind[i]];  steps <- NULL
                 doppelpunkt <- which( sapply(1:nchar(term),FUN=function(ii){u <- substr(term,ii,ii); b <- u==":"  }) )
                 term <- substr(term,doppelpunkt+2,nchar(term))
                 cat(paste("Found TERM ",i,": '",term,"' \n",sep=""))
                 all.terms <- c(all.terms,term)                                 ### Dies dient nur dazu, hinterher die Liste mit ausgelesenen Tabellen beschriften zu koennen.
                 bereich <- (ind[i]+6) : (grenzen[i] -2)                        ### Dies der Bereich, der ausgewaehlt werden muss
                 namen   <- c("No.", strsplit(input.all[bereich[1]-2]," +")[[1]][-1])
                 namen   <- gsub("\\^","",namen)
                 index   <- grep("CI",namen)                                    ### Wenn ein "CI" als Spaltenname erscheint, muessen daraus im R-Dataframe zwei Spalten werden!
                 if(length(index) > 0)  {
                    for (ii in 1:length(index)) {
                         namen  <- c(namen[1:index[ii]], "CI",namen[(index[ii]+1):length(namen)] )}}
                 input.sel  <- eatTools::crop( input.all[bereich] )             ### Textfile wird reduziert, und voranstehende und abschliessende Leerzeichen werden entfernt
                 input.sel  <- gsub("\\(|)|,"," ",input.sel)                    ### entferne Klammern und Kommas (wenn's welche gibt)
                 input.sel  <- gsub("\\*    ", "  NA", input.sel)               ### hier: gefaehrlich: wenn mittendrin Werte fehlen, wuerde stringsplit eine unterschiedliche Anzahl Elemente je Zeile finden
                 foo        <- strsplit(input.sel," +")                         ### und die fehlenden Elemente stets ans Ende setzen. Fatal!
                 maxColumns <- max(sapply(foo, FUN=function(ii){ length(ii)}))  ### Gefahr 2: '*' bezeichnet fixierte Parameter, die keinen Standardfehloeer haben. Manchmal steht aber trotzdem einer da (z.B. in DIF). Ersetzung soll nur stattfinden, wenn mehr als vier Leerzeichen hinterher
                 nDifferentColumns <- length( table(sapply(foo, FUN=function(ii){ length(ii)  })))
                 maxColumns <- which( sapply(foo, FUN=function(ii){ length(ii) == maxColumns  }) ) [1]
    ### untere Zeile: WICHTIG! wo stehen in der Zeile mit den meisten nicht fehlenden Werten Leerzeichen?
                 foo.2      <- which( sapply(1:nchar(input.sel[maxColumns]),FUN=function(ii){u <- substr(input.sel[maxColumns],ii,ii); b <- u==" "  }) )
                 foo.3      <- diff(foo.2)                                      ### zeige die Position des letzten Leerzeichens vor einem Nicht-Leerzeichen
                 foo.3      <- foo.2[foo.3 !=1]                                 ### suche nun in jeder Zeile von input.sel: ist das Zeichen zwei Stellen nach foo.3 ein Leerzeichen? Wenn ja: NA!
                 ESTIMATE   <- which( sapply(1:nchar(input.all[ind[i] + 4] ),FUN=function(ii){u <- substr(input.all[ind[i] + 4],ii,ii+7); b <- u=="ESTIMATE"  }) )
                 foo.3      <- foo.3[foo.3>(ESTIMATE-3)]                        ### Achtung: das alles soll aber nur fuer Spalten beginnen, die hinter ESTIMATE stehen! (missraet sonst fuer Produktterme, z.B. "item*sex")
                 if(nDifferentColumns>1) {
                    if(length(foo.3)>0) {                                       ### Und nochmal: das soll NUR geschehen, wenn es in mindestens einer Zeile nicht die vollstaendige (=maximale) Anzahl von Elementen gibt!
                       for (ii in 1:length(input.sel)) {                        ### also wenn nDifferentColumns groesser als EINS ist (kleiner darf es nicht sein)
                            for (iii in 1:length(foo.3)) {
                                 if(substr( input.sel[ii], foo.3[iii] + 2 , foo.3[iii] + 2 ) == " ") {input.sel[ii] <- paste(substr(input.sel[ii],1,foo.3[iii]), "NA", substring(input.sel[ii],foo.3[iii]+3) , sep="")}}}}
                    if(length(foo.3)==0) {cat(paste("There seem to be no values in any columns behind 'ESTIMATE'. Check outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep=""))}}
                 input.sel <- strsplit(input.sel," +")
                 if(length(input.sel[[1]]) == 0 ) {cat(paste("There seem to be no valid values associated with term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep=""))
                                                   all.terms <- all.terms[-i]}
                 if(length(input.sel[[1]]) > 0 ) {
                    referenzlaenge <- max (sapply( input.sel, FUN=function(ii ){  length(ii)    }) )
                    if(referenzlaenge < length(namen) ) {
                       cat(paste("Several columns seem to be empty for term '",all.terms[length(all.terms)],"' in file: '",file,"'.\n",sep=""))
    ### bloeder spezialfall: wenn dif-Analyse mit 'compute.fit=FALSE' gemacht wurde, fehlen die infit-Spalten ... die eigentlich notwendige zusaetzliche spalte 'add.column' wird dann nicht eingefuegt. Finde raus, ob das der Fall ist
                       head <- eatTools::crop(input.all[bereich[1]-2])          ### Ueberschrift
                       leerz<- gregexpr(" ", head)[[1]]                         ### suche: wo beginnt der zweite Block aufeinanderfolgender leerzeichen?
                       leerd<- which ( diff ( leerz) > 1 )[2]
                       vgl  <- length(strsplit ( eatTools::crop(substr(input.all[bereich[1]], 1, leerd)), split = " +")[[1]])
                       if ( vgl == 4 ) {
                            namen <- c(namen[1:2], "add.column1", namen[3:(referenzlaenge-1)])
                       }  else  {
                            referenzlaenge <- length(namen)
                       }
                    }
    ### Ende spezialfall
                    if(referenzlaenge > length(namen) ) {
                       if(referenzlaenge == length(namen) + 1) {
                          cat(paste("There seem to be one more column than columns names. Expect missing column name before 'ESTIMATE'. \nCheck outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep=""))
                          ind.name <- which(namen == "ESTIMATE")
                          namen    <- c(namen[1:ind.name-1], "add.column",namen[ind.name:length(namen)])}
                       if(referenzlaenge >  length(namen) + 1) {
                          cat(paste("There seem to be more columns than names for it. Check outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep=""))
                          namen<- c(namen, rep("add.column",referenzlaenge-length(namen) )) }}
                    input.sel  <- t(sapply(input.sel, FUN=function(ii){ c(ii, rep(NA,referenzlaenge-length(ii))) }))
                    colnames(input.sel) <- namen                                ### untere Zeile: entferne eventuelle Sternchen und wandle in Dataframe um!
                    input.sel  <- suppressWarnings(eatTools::asNumericIfPossible(data.frame( gsub("\\*","",input.sel), stringsAsFactors = FALSE), force.string = FALSE))
                    results.sel<- data.frame(input.sel,filename=file,stringsAsFactors = FALSE)
                    if(is.na(as.numeric(results.sel$ESTIMATE[1]))) {cat(paste("'ESTIMATE' column in Outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"' does not seem to be a numeric value. Please check!\n",sep=""))}
                    if(!missing(dif.term)) {                                    ### Der absolute DIF-Wert ist 2 * "Betrag des Gruppenunterschieds". Fuer DIF muessen ZWEI Kriterien erfuellt sein:
                       if(all.terms[length(all.terms)] == dif.term) {           ### Der absolute DIF-Wert muss groesser als 'abs.dif.bound' (z.B. 0.6) und zugleich signifikant groesser als 'sig.dif.bound' (z.B. 0.3) sein
                          cat(paste("Treat '",all.terms[length(all.terms)],"' as DIF TERM.\n",sep=""))
                          results.sel <- data.frame(results.sel,abs.dif = 2*results.sel$ESTIMATE,stringsAsFactors=FALSE)
                          konfNiveau  <- round(100*p.value)                     ### Das bedeutet, fuer Werte groesser 0.6 darf 0.3 NICHT im 90 bzw. 95%-Konfidenzintervall liegen. Nur dann haben wir DIF!
                          results.sel[,paste("KI.",konfNiveau,".u",sep="")] <- results.sel$abs.dif-2*abs(qnorm(0.5*(1-p.value)))*results.sel$ERROR
                          results.sel[,paste("KI.",konfNiveau,".o",sep="")] <- results.sel$abs.dif+2*abs(qnorm(0.5*(1-p.value)))*results.sel$ERROR
                          results.sel[,paste("sig.",konfNiveau,sep="")] <- ifelse(abs(results.sel[,"abs.dif"])>abs.dif.bound & abs(results.sel[,paste("KI.",konfNiveau,".u",sep="")])>sig.dif.bound & abs(results.sel[,paste("KI.",konfNiveau,".o",sep="")])>sig.dif.bound,1,0)
                          results.sel$filename <- file
                          if(split.dif==TRUE) {results.sel <- results.sel[1:(dim(results.sel)[1]/2),]
                                               if(dim(results.sel)[1]!=dim(results.sel)[1]) {warning("missing variables in DIF table.")}}}}
                 all.output[[i]] <- results.sel}}
              if(!missing(dif.term)) {if(sum(all.terms==dif.term)==0) {cat(paste("Term declarated as DIF: '",dif.term,"' was not found in file: '",file,"'. \n",sep=""))  }}
              names(all.output) <- all.terms
    ### ggf. Regressionsparameter einlesen!
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
    ### Kovarianz-/ Korrelationsmatrix einlesen: schwierig, also Trennen nach ein- vs. mehrdimensional. Eindimensional: zweimal "-----" zwischen Beginn und Ende des COVARIANCE-Statements
              korStart <- grep("COVARIANCE/CORRELATION MATRIX", input.all)
              korEnd   <- grep("An asterisk next", input.all)
              korEnd   <- min(korEnd[korEnd > korStart])
              korStriche <- grep("-----",input.all)
              korStriche <- korStriche[korStriche > korStart & korStriche < korEnd]
              if(length(korStriche) == 2) {                                     ### eindimensional!
                 varRow    <- grep("Variance", input.all)
                 variance  <- as.numeric( unlist( lapply(strsplit(input.all[varRow]," +"), FUN=function(ll) {ll[length(ll)]}) ) )
                 names(variance) <- "variance"
                 all.output$cov.structure <- variance
              }
              if(length(korStriche) > 2) {                                      ### mehrdimensional!
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
    ### Reliabilitaetsindices einlesen
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
            input <- strsplit( gsub("\\\t"," ",eatTools::crop(input)), "/\\*")  ### Hier ist es wichtig, gsub() anstelle von sub() zu verwenden! sub() loescht nur das erste Tabulatorzeichen
            ret   <- data.frame ( do.call("rbind", strsplit( eatTools::crop(unlist(lapply(input, FUN = function ( l ) {l[1]}))), " +")), stringsAsFactors = FALSE)
            nameI <- eatTools::crop(eatTools::removePattern ( eatTools::crop( eatTools::crop(unlist(lapply(input, FUN = function ( l ) {l[length(l)]}))), char = "item"), pattern = "\\*/"))
            ret   <- data.frame ( Case= as.numeric(ret[,1]), item = nameI, parameter= as.numeric(ret[,2]) ,stringsAsFactors = FALSE)
            return(ret)}

get.itn <- function(file)  {
            input <- scan(file, what = "character", sep="\n", quiet = TRUE)
            ind.1 <- grep("==========",input)
            items <- grep( "item:", input )
            diff.last <- ind.1[length(ind.1)-1] - items[length(items)] + 4
            items <- cbind(1:length(items),items,c(diff(items),diff.last))      ### dort wo diff(items) != 13 , ist das entsprechende Item partial credit. (Fuer das letzte Item ist das komplizierter, da length(diff(items))<length(items).    )
            ind.2 <- gregexpr(":", input[items[,2]])                            ### Folgende Zeilen dienen dazu zu pruefen, ob DIFs in der Tabelle vorkommen oder nicht (falls ja, dann gibt es zwei Doppelpunkte pro input[items[,2]]
            ind.3 <- unlist(ind.2)                                              ### Dann ist ind.3 auch doppelt so lang wie ind.2, weil jedes Element aus ind.2 ein Vektor mit zwei Elementen ist
            ind.3 <- matrix(ind.3,length(ind.2),byrow=T)
            item.namen <- substr(input[items[,2]], ind.3[,dim(ind.3)[2]]+1+nchar(as.character(items[,1])),100)
            item.namen <- gsub(" ","",item.namen)                               ### Leider funktioniert gsub() nicht fuer Klammern, da diese fuer regular expression reserviert sind, aber...
            item.namen <- gsub("\\)","",item.namen); item.namen <- gsub("\\(","",item.namen)
            if(dim(ind.3)[2]>1)                                                 ### kommen DIFs din vor? Ja, falls Bedingung TRUE
              {stopifnot(length(table(ind.3[,1]))==1)                           ### sollte 1 sein; da es immer dieselbe DIF-Variable mit ergo derselben Zeichenlaenge ist.
               dif.name <- rep(substr(input[items[,2]], 1, ind.3[,1]-1),(items[,3]-11))                          ### Auslesen der Variablennamen fuer DIF
               dif.value <- rep(as.numeric(substr(input[items[,2]], ind.3[,1]+1, ind.3[,1]+1)),(items[,3]-11))}  ### Auslesen des Wertes der DIF-Variablen
            zeilen <- list(); reihe <- NULL                                     ### Was geschieht oben? Die DIF-Variable wird fuer Item repetiert, und zwar zweimal, wenn es ein normales, dreimal, wenn es ein partial credit-Item ist. Die entsprechende Information steht in items[,3]; vgl.: rep(1:4,1:4)
            for (i in 1:dim(items)[1])                                          ### finde die Zeilen fuer jedes Item
                {zeilen[[i]] <- (items[i,2]+7) : (items[i,2]+ (items[i,3]-5) )  ### kein partial credit: beginne sieben Zeilen unter "item:" und ende bei acht Zeilen (= 13-5) unter "item:". Fuer partial credit, ende items[i,3]-5 Zeilen unter "items:"
                 cases       <- gsub("NA ","NA",input[zeilen[[i]]])             ### Untere Zeile: Korrektur, wenn die zwei Datenzeilen leere felder enthalten (NA wird nachtraeglich eingetragen)
                 cases <- gsub("_BIG_ ","NA",cases)
                 cases <- gsub("_BIG_","NA",cases)
                 if(length(table(sapply(1:length(cases),FUN=function(ii){length(unlist(strsplit(cases[ii]," +"))) }) ) )>1 )
                   {cases <- gsub("          ","    NA    ",cases)}             ### Perfekt! ueberall dort, wo zehn Leerzeichen infolge stehen, muss eine Auslassung sein! Hier wird ein Ersetzung gemacht!
                 cases       <- data.frame( matrix ( unlist( strsplit(eatTools::crop(gsub(" +"," ", cases))," ") ), nrow=length(zeilen[[i]]),byrow=T ) , stringsAsFactors=F)
                 ind         <- grep("\\)",cases[1,]); cases[,ind] <- gsub("\\)","",cases[,ind] )
                 cases       <- data.frame(cases[,1:(ind-1)],matrix(unlist(strsplit(cases[,6],"\\(")),nrow=length(zeilen[[i]]),byrow=T),cases[,-c(1:ind)],stringsAsFactors=F)
                 for(jj in 1:ncol(cases)) {cases[,jj] <- as.numeric(cases[,jj])}
                 colnames(cases) <- c("Label","Score","Abs.Freq","Rel.Freq","pt.bis","t.value","p.value",paste(rep(c("PV1.Avg.","PV1.SD."),((ncol(cases)-7)/2) ),rep(1:((ncol(cases)-7)/2),each=2),sep=""))
                 threshold.zeile   <- input[items[i,2]+2]; threshold <- NULL; delta <- NULL
                 bereich <- ifelse( (items[i,3]-12)<1,1,(items[i,3]-12))        ### Sicherheitsbedingung, falls Variable nur eine Kategorie hat
                 if((items[i,3]-12)<1) {cat(paste("Item",i,"hat nur eine Antwortkategorie.\n"))}
                 for (j in 1: bereich )
                     {threshold  <- c(threshold ,as.numeric(substr(threshold.zeile,  6*j+16,6*j+21)))
                      delta      <- c(delta,     as.numeric(substr(input[items[i,2]+3],6*j+13,6*j+18)))}
                 while(length(threshold) < nrow(cases)) {threshold <- c(threshold,NA)}
                 while(length(delta) < nrow(cases)) {delta <- c(delta,NA)}
                 item.p <- NA                                                   ### Manchmal kann kein p-wert bestimmt werden. Wenn doch, wird das NA ueberschrieben
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

### Funktion generiert und schreibt Syntax. Voraussetzung: Datensatz ist da.
### gen.syntax wird von defineModel aufgerufen,  keine user-level Funktion!
### Funktion ist von Restrukturierung ("defineModel" statt "prep.conquest") erstmal nicht betroffen, alles bleibt wie's ist
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
                   if(is.null(Title))   {                                       ### wenn kein Titel gesetzt, erstelle ihn aus Sys.getenv()
                      all.inf  <- Sys.getenv()
                      Title    <- paste("Analysis name: ",Name, ", User: ",all.inf["USERNAME"],", Computername: ",all.inf["COMPUTERNAME"],", ", R.version$version.string , ", Time: ",date(),sep="")}
                   converge <- paste("0",substring(as.character(converge+1),2),sep="")
                   deviancechange <- paste("0",substring(as.character(deviancechange+1),2),sep="")
                   syntax    <- gsub("####hier.title.einfuegen####",Title,mustersyntax)
                   if(is.null(n.plausible))   {n.plausible <- 0}  ; if(is.na(n.plausible))     {n.plausible <- 0}
                   if(n.plausible == 0 )     {                                  ### wenn Anzahl PVs = 0 oder NULL, loesche Statement; andernfalls: setze Anzahl zu ziehender PVs ein!
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
                   beginn    <- NULL                                            ### setze "beginn" auf NULL. Wenn DIF-Variablen spezifiziert sind, wird "beginn" bereits
                   if(length(namen.all.hg)>0)    {                              ### untere Zeile: wieviele "character" haben Hintergrundvariablen?
                     all.hg.char.kontroll <- all.hg.char
                     all.hg.char <- sapply(namen.all.hg, FUN=function(ii) {max(nchar(as.character(na.omit(daten[,ii]))))})
                     stopifnot(all(all.hg.char == all.hg.char.kontroll))        ### Trage nun die Spalten in das Format-Statement ein: Fuer ALLE expliziten Variablen
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
                         ind.model <- grep("model item", syntax)                ### Aendere model statement
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
                   if(length(all.Names[["weight.var"]])>0)  {                   ### Method wird erst hier gesetzt, weil sie davon abhaengt, ob es ein HG-Modell gibt
                      ind.4   <- grep("caseweight",syntax)
                      syntax[ind.4] <- paste( syntax[ind.4], " ", tolower(all.Names[["weight.var"]]),";",sep="") }
                   if(length(all.Names[["group.var"]])>0) {
                       ind.3   <- grep("^group$",syntax)
                       stopifnot(length(ind.3) == 1)
                       syntax[ind.3] <- paste(eatTools::crop(paste( c(syntax[ind.3], tolower(all.Names[["group.var"]])), collapse=" ")),";",sep="")
                       ### gebe gruppenspezifische Descriptives
                       add.syntax.pv  <- as.vector(sapply(all.Names[["group.var"]], FUN=function(ii) {paste("descriptives !estimates=pv, group=",tolower(ii)," >> ", Name,"_",tolower(ii),"_pvl.dsc;",sep="")} ))
                       add.syntax.wle <- as.vector(sapply(all.Names[["group.var"]], FUN=function(ii) {paste("descriptives !estimates=wle, group=",tolower(ii)," >> ", Name,"_",tolower(ii),"_wle.dsc;",sep="")} ))
                       ind.3    <- grep("quit",syntax)
                       stopifnot(length(ind.3)==1)
                       syntax   <- c(syntax[1:(ind.3-1)],add.syntax.pv, add.syntax.wle, syntax[ind.3:length(syntax)]) }
                   if(is.null(beginn)) {beginn <- ID.char+1}
                   syntax[ind] <- paste(syntax[ind], "responses ",beginn,"-",beginn-1+var.char*ncol(data.frame(daten[,all.Names[["variablen"]]],stringsAsFactors = FALSE)),";",sep="")
                   if(var.char>1)  {                                            ### Items haben mehr als eine Spalte Stelligkeit (Conquest-Handbuch, S.177)
                      syntax[ind] <- paste(gsub(";","",syntax[ind]), " (a",var.char,");",sep="")}
                   score.statement <- .writeScoreStatementMultidim (data=daten, itemCols=all.Names[["variablen"]], qmatrix=model, columnItemNames = 1 ,use.letters=use.letters, allowAllScoresEverywhere = allowAllScoresEverywhere )
                   expected.nodes  <- nodes^(ncol(model)-1)
                   if(expected.nodes>3500 & method != "montecarlo") {cat(paste("Specified model probably will use ",expected.nodes," nodes. Choosen method ",method," may not appropriate. Recommend to use 'montecarlo' instead.\n",sep=""))}
                   ind <- grep("labels ",syntax)
                   stopifnot(length(ind)==1)
                   syntax <- c(syntax[1:ind],score.statement,syntax[(ind+1):length(syntax)])
                   if(length(all.Names[["HG.var"]])==0) {                       ### wenn kein HG-model, loesche entsprechende Syntaxzeilen
                      ind.2 <- grep("^regression$",syntax)
                      stopifnot(length(ind.2)==1)
                      syntax <- syntax[-ind.2]
                      ind.3 <- grep("export reg_coefficients",syntax)
                      stopifnot(length(ind.3)==1)
                      syntax <- syntax[-ind.3] }
                   if(length(all.Names[["group.var"]]) ==0) {                   ### wenn keine Gruppen definiert, loesche Statement
                      ind.3 <- grep("^group$",syntax)
                      stopifnot(length(ind.3)==1)
                      syntax <- syntax[-ind.3]}
                   if(length(all.Names[["weight.var"]]) ==0) {                  ### wenn keine Gewichte definiert, loesche Statement
                      ind.4 <- grep("^caseweight$",syntax)
                      stopifnot(length(ind.4)==1)
                      syntax <- syntax[-ind.4]}
                   if(match.arg(equivalence.table) == "NULL") {                 ### wenn keine Equivalence-Statement definiert, loesche Zeile
                      ind.5   <- grep("^equivalence",syntax)
                      stopifnot(length(ind.5)==1)
                      syntax <- syntax[-ind.5]}
                   if(is.null(seed)) {                                          ### wenn keine seed-Statement definiert, loesche Zeile
                      ind.7   <- grep("^set seed",syntax)
                      stopifnot(length(ind.7)==1)
                      syntax <- syntax[-ind.7]}
                   if(n.plausible == 0)     {                                   ### wenn Anzahl PVs = 0 oder NULL, loesche Statement
                      ind.6   <- grep("^show cases! estimates=latent", syntax)
                      stopifnot(length(ind.6) == 1)
                      syntax  <- syntax[-ind.6]}
                   if(anchored == FALSE) {ind.2 <- grep("anchor_parameter",syntax)# wenn keine ANKER gesetzt, loesche entsprechende Syntaxzeile
                                        syntax <- syntax[-ind.2]}
                   if(anchored == TRUE)  {ind.2 <- grep("^set constraints",syntax)# wenn ANKER gesetzt, setze constraints auf "none"
                                        if(match.arg(constraints) != "none") { cat("Anchorparameter were defined. Set constraints to 'none'.\n")}
                                        syntax[ind.2]  <- "set constraints=none;"}
                   if(!all(sapply(export, inherits, what="logical"))) {stop("All list elements of argument 'export' have to be of class 'logical'.")}
                   export <- as.list(userSpecifiedList ( l = export, l.default = export.default ))
                   weg <- names(export[which(export == FALSE)])
                   if(length(weg)>0)    {                                       ### hier wird, was nicht exportiert werden soll, aus Syntax geloescht.
                      for (ii in seq(along=weg) ) {
                           ind.x <- grep(paste("export ", weg[ii], sep=""), syntax)
                           stopifnot(length(ind.x) == 1)
                           syntax <- syntax[-ind.x]}}
                   if(export["history"] == TRUE)  {
                      if(!is.null(conquest.folder))  {
                         checkWhetherConquestExeExists(pkgname="eatModel")
                         cq.version <- getConquestVersion( path.conquest = conquest.folder, path.temp = pfad)
                         if(cq.version < date::as.date("1Jan2007") )   {
   									      ind.3 <- grep("^export history",syntax)               ### wenn Conquest aelter als 2007, soll history geloescht werden,
                           stopifnot(length(ind.3) == 1 )                       ### auch dann, wenn der Benutzer History ausgeben will
                           syntax <- syntax[-ind.3]
                         }
                      }
                      if(is.null(conquest.folder)) {warning("Conquest folder was not specified. Unable to detect Conquest version. When you propose to use 2005 version,\nhistory statement will invoke to crash Conquest analysis. Please remove history statement manually if you work with 2005 version.")} }
                   write(syntax,file.path(pfad,paste(Name,".cqc",sep="")),sep="\n")}

### das oeffnet ein Menue, um die conquest.exe zu finden, falls es sie nicht gibt
checkWhetherConquestExeExists <- function (pkgname) {
     root <- system.file(package = pkgname)
     if ( !file.exists(file.path(root, "exec"))) {dir.create(file.path(root, "exec"))}
     if ( !file.exists( system.file("exec", "console_Feb2007.exe", package = pkgname) )) {
           if ( !file.exists("i:/Methoden/00_conquest_console/console_Feb2007.exe") ) {
               packageStartupMessage("Cannot find conquest executable file. Please choose manually.")
               fname <- file.choose()
           }  else  {
               fname <- "i:/Methoden/00_conquest_console/console_Feb2007.exe"
           }
           if ( nchar(fname)>0) { foo <- file.copy(from = fname, to = file.path(root, "exec", "console_Feb2007.exe") ) }
     }
}

anker <- function(lab, prm, qMatrix, domainCol, itemCol, valueCol, multicore )  {
                  stopifnot(ncol(lab)==2)
                  if ( !ncol(prm) == 2 )   {                                    ### wenn itemliste nicht unique ... 'domain'-Spalte kann ausgelassen werden
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
                  if(!is.null(multicore) && multicore == TRUE) {                ### nur fuer multicore umweg ueber capture.output gehen, damit die messages auch im multicore-betrieb kommen, sonst landen sie nicht auf der konsole
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

### columnItemNames         ... in welcher Spalte der q-Matrix stehen Itemnamen?
### columnsDimension        ... in welchen Spalten der Q-Matrix stehen die Dimensionen?
###                             Default: in erster Spalte stehen Itemnamen, in allen uebrigen Spalten stehen Indikatoren fuer Dimensionen
.writeScoreStatementMultidim <- function(data, itemCols, qmatrix, columnItemNames = 1 ,columnsDimensions = -1, use.letters=use.letters , allowAllScoresEverywhere) {
            n.dim      <- (1:ncol(qmatrix) )[-columnItemNames]                  ### diese Spalten bezeichnen Dimensionen. untere Zeile: Items, die auf keiner Dimension laden, werden bereits in prep.conquest entfernt. hier nur check
            stopifnot(length( which( rowSums(qmatrix[,n.dim,drop = FALSE]) == 0))==0)
      	    if(length(setdiff(names(eatTools::tableUnlist(qmatrix[,-1, drop = FALSE])), c("0","1"))) > 0 )  {
               cat("Found unequal factor loadings for at least one dimension. This will result in a 2PL model.\n")
               for (u in 2:ncol(qmatrix)) {qmatrix[,u] <- as.character(round(qmatrix[,u], digits = 3))}
            }                                                                   ### obere Zeile: Identifiziere Items mit Trennschaerfe ungleich 1.
            stopifnot(all(qmatrix[,1] == itemCols))                             ### untere Zeile: Items im Datensatz, aber nicht in Q-Matrix? wird bereits in prep.conquest behandelt
            cat(paste("Q matrix specifies ",length(n.dim)," dimension(s).\n",sep=""))
            stopifnot(length(setdiff(colnames(data[,itemCols]),  qmatrix[,columnItemNames]) )==0)
            unique.patter <- qmatrix[which(!duplicated(do.call("paste", qmatrix[,-1, drop = FALSE] ))), -1, drop = FALSE]
            colnames(unique.patter) <- paste("Var",1:ncol(unique.patter), sep="")## obere Zeile: Finde alle uniquen Pattern in qmatrix! Jedes unique Pattern muss in Conquest einzeln adressiert werden!
            score.matrix  <- data.frame(score=1, unique.patter, matrix(NA, nrow= nrow(unique.patter), ncol=length(itemCols), dimnames=list(NULL, paste("X",1:length(itemCols),sep=""))),stringsAsFactors = FALSE)
            scoreColumns  <- grep("^Var",colnames(score.matrix))
            for (i in 1:length(itemCols))  {                                    ### gebe alle Items auf den jeweiligen Dimensionen
               qmatrix.i    <- qmatrix[qmatrix[,columnItemNames] == itemCols[i],]## auf welcher Dimension laedt Variable i? Untere Zeile: in diese Zeile von score.matrix muss ich variable i eintragen
               matchRow     <- which(sapply ( 1:nrow(score.matrix) , function(ii) {all ( as.numeric(qmatrix.i[,n.dim]) == as.numeric(score.matrix[ii,scoreColumns])) }))
               stopifnot(length(matchRow) == 1)
               matchColumn  <- min(which(is.na(score.matrix[matchRow,])))       ### in welche spalte von Score.matrix muss ich variable i eintragen?
               stopifnot(length(matchColumn) == 1)
               score.matrix[matchRow,matchColumn] <- i
		        }
            rowsToDelete <- which(is.na(score.matrix[, max(scoreColumns) + 1])) ### welche Zeilen in Score.matrix koennen geloescht werden?
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
            ### Pruefung, ob "tranformation" des score-statements ok ist
            control <- lapply(kollapse.string,FUN=function(ii) {eval(parse(text=paste("c(",gsub("-",":",ii),")",sep="")))})
            if (!all(unlist(lapply(1:length(control), FUN=function(ii) {all(kollapse[[ii]] == control[[ii]])})))) {
                cat("Error in creating score statement.\n")
			      }
            score.matrix <- data.frame(prefix="score",score.matrix[,c(1,scoreColumns)],items="! items(",kollapse.string=unlist(kollapse.string),suffix=");",stringsAsFactors=F)
            score.statement <- sapply(1:nrow(score.matrix), FUN=function(ii) { paste(score.matrix[ii,],collapse=" ")})
            return(score.statement) }


### Hilfsfunktion fuer .writeScoreStatementMultidim()
fromMinToMax <- function(dat, score.matrix, qmatrix, allowAllScoresEverywhere, use.letters)    {
                all.values <- plyr::alply(as.matrix(score.matrix), .margins = 1, .fun = function(ii) {sort(names(eatTools::tableUnlist(dat[,na.omit(as.numeric(ii[grep("^X", names(ii))])), drop = FALSE])) ) })
                if ( length(all.values) > 1) {                                  ### obere Zeile: "alply" ersetzt "apply"! http://stackoverflow.com/questions/6241236/force-apply-to-return-a-list
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

### Hilfsfunktionen fuer gen.syntax
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

### Funktion komplett neu geschrieben, 1. Dezember 2011; nutzt Funktion "table.muster"
desk.irt <- function(daten, itemspalten, na=NA,percent=FALSE,reduce=TRUE,codebook=list(datei=NULL,item=NULL,value=NULL,lab=NULL, komp=NULL), quiet = FALSE ) {
             daten <- eatTools::makeDataFrame(daten)
             if(!missing(itemspalten)) {daten <- daten[,itemspalten,drop=FALSE]}
             if (is.na(na[1])==FALSE) {                                         ### wenn spezifiziert, werden hier missings recodiert
                 recode.statement <- paste(na,"= NA",collapse="; ")
                 daten            <- data.frame(sapply(daten,FUN=function(ii) {car::recode(ii,recode.statement)}),stringsAsFactors=FALSE)
             }
             specific.codes <- lapply(daten,function(ii){NULL})                 ### definiert ggf. Spezifische Codes, nach denen je Variable gesucht werden soll
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

### angelehnt an das Skript von Alexander Robitzsch, "R_Skalierung.odt", Seite 11
item.diskrim <- function(daten, itemspalten, streng = TRUE) {
                 if(!missing(itemspalten))  {daten <- daten[,itemspalten]}      ### Trennschaerfe ist eigentlich Korrelation des Items mit dem Summenscore ohne dieses Item.
                 trenn <- suppressWarnings(eatTools::pwc(daten))                ### Dieses macht die Option "streng = T"; die andere berechnet Korrelation mit Summenscore einschliesslich dieses Items
                 if(streng) {return(data.frame(item.name=trenn[,"item"],item.diskrim = trenn[,"partWholeCorr"],stringsAsFactors = FALSE))} else {return(data.frame(item.name=trenn[,"item"],item.diskrim = trenn[,"corr"],stringsAsFactors = FALSE))}}

### foo <- prepRep( T.t1t2, dfrT1P, dfrT2P)
prepRep <- function ( calibT2, bistaTransfT1, bistaTransfT2, makeIdsUnique = TRUE) {
           if ( !inherits(calibT2, "transfBista" )) { stop("'calibT2' object must be of class 'transfBista'.\n")}
           if ( !inherits(bistaTransfT1, "transfBista" )) { stop("'bistaTransfT2' object must be of class 'transfBista'.\n")}
           if ( !inherits(bistaTransfT2, "transfBista") ) { stop("'bistaTransfT2' object must be of class 'transfBista'.\n")}
           if (!nrow(calibT2[["itempars"]]) < nrow(bistaTransfT1[["itempars"]])) { stop("Mismatch between 'calibT2' and 'bistaTransfT1'. \n")}
           if (!nrow(calibT2[["itempars"]]) < nrow(bistaTransfT2[["itempars"]])) { stop("Mismatch between 'calibT2' and 'bistaTransfT2'. \n")}
     ### check: heissen die ID-Variablen etc. in beiden Datensaetzen gleich? ... falls nicht, misslingt unten das 'rbind' ... ggf. neue ID (falls nicht identisch in beiden Datensaetzen)
           idT1<- unique(bistaTransfT1[["all.Names"]][which(bistaTransfT1[["all.Names"]][,"par"] == "ID"),"derived.par"])
           idT2<- unique(bistaTransfT2[["all.Names"]][which(bistaTransfT2[["all.Names"]][,"par"] == "ID"),"derived.par"])
           stopifnot(length(idT1)==1, length(idT2)==1)
           if ( idT1 != idT2 ) {
                warning(paste0("ID variables do not match between t1 and t2. ID for t1: '",idT1,"'. ID for t2: '",idT2,"'. \n    IDs will be unified with '",idT1,"'."))
                recStat <- paste ( "'", idT2 , "' = '", idT1, "'", sep="")
                colnames ( bistaTransfT2[["personpars"]] ) <- car::recode ( colnames ( bistaTransfT2[["personpars"]] ), recStat)
           }
     ### finde Spalten mit Linkingfehlern
           lc  <- colnames( calibT2[["personpars"]] ) [grep("^linking", colnames(calibT2[["personpars"]]) )]
           if(length(lc)==0) { stop("No columns with linking error information found in 'calibT2'.\n")}
     ### benenne spalten in 'trend...' um
           lcn <- paste("trend", eatTools::removePattern(string = lc, pattern = "linking"), sep="")
           colnames( calibT2[["personpars"]] ) [grep("^linking", colnames(calibT2[["personpars"]]) )] <- lcn
     ### suche Spalten zum Mergen
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
     ### Reduziere Kalibrierungs-'datensatz' auf das Noetigste
           red <- calibT2[["personpars"]][,c(toM,  lcn)]
           red <- red[!duplicated(red),]
           dat1<- data.frame ( trend = "T1" , merge ( bistaTransfT1[["personpars"]], red, by = toM, all = TRUE))
     ### checks (sollten eigentlich ueberfluessig sein)
           stopifnot ( nrow(dat1) == nrow(bistaTransfT1[["personpars"]]))
           dat2<- data.frame ( trend = "T2" , merge ( bistaTransfT2[["personpars"]], red, by = toM, all = TRUE))
           stopifnot ( nrow(dat2) == nrow(bistaTransfT2[["personpars"]]))
     ### IDs unique machen (wenn gewuenscht)
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
           it[,"est"] <- rowSums(it[,c("est", "estOffset")], na.rm = TRUE)      ### untere Zeilen: wenn 1pl und 2pl gemeinsam im resultsobjekt auftauchen, gibt es fuer 1pl keinen
           if ( !"estSlope" %in% colnames(it) ) { it[,"estSlope"] <- 1 }        ### slope parameter; die werte sind NA. Zum Plotten muessen sie daher fuer das Raschmodell auf 1 gesetzt werden
           if ( length(which(is.na(it[,"estSlope"]))) > 0) { it[which(is.na(it[,"estSlope"])), "estSlope"] <- 1 }
           eapA<- eapFromRes (resultsObj)                                       ### eap fuer alle; muss wideformat haben!!!
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
     ### Plotten findet fuer jedes Item separat statt
           pl  <- by ( data = it, INDICES = it[,c("model", "item")], FUN = function ( i ) {
                  xlm <- c(i[["est"]]+2, i[["est"]]-2)
                  anf <- -6                                                     # anf <- if ( min(xlm) < -4 ) { anf <- floor(min(xlm)) } else { anf  <- -4}
                  ende<- 6                                                      # ende<- if ( max(xlm) >  4 ) { ende<- ceiling(max(xlm)) } else { ende <- 4}
                  x   <- seq ( anf, ende, l = 400)
                  y   <- exp( i[["estSlope"]]*x - i[["est"]] ) / (1+exp( i[["estSlope"]]*x - i[["est"]] ))
                  plot (x, y, type = "l", main = paste("Item '",as.character(i[["item"]]),"'\n\n",sep=""), xlim = c(-6,6), ylim = c(0,1), xlab = "theta", ylab = "P(X=1)", col = "darkred", cex = 8, lwd = 2)
                  graphics::mtext( paste("Model = ",i[["model"]],"  |  Dimension = ",i[["dimension"]], "  |  difficulty = ",round(i[["est"]], digits = 3),"  |  Infit = ",round(i[["infit"]], digits = 3),"\n",sep=""))
                  eap <- eapA[intersect ( which (eapA[,"dimension"] == i[["dimension"]]) , which (eapA[,"model"] == i[["model"]])),]
                  if ( inherits(defineModelObj, "defineMultiple")) {            ### Problem: je nachdem ob modelle gesplittet wurden oder nicht, muss der Itemdatensatz woanders gesucht werden ... Hotfix
                       woIst<- which ( lapply ( defineModelObj, FUN = function ( g ) {   g[["analysis.name"]] == i[["model"]] }) == TRUE)
                       stopifnot(length(woIst) == 1)
                       dat  <-defineModelObj[[woIst]][["daten"]]
                  }  else  {
                       dat  <- defineModelObj[["daten"]]
                  }
     ### Hotfix: ID namen identifizieren
                  id  <- unique(resultsObj[intersect(which(resultsObj[,"type"] == "tech"), which(resultsObj[,"par"] == "ID")),"derived.par"])
                  stopifnot(length(id)==1)
                  prbs<- na.omit ( merge ( dat[,c( "ID", as.character(i[["item"]]))], eap[,c( id, "EAP")], by.x = "ID", by.y = id))
                  anz <- round ( nrow(prbs) / personsPerGroup ) + 1             ### mindestens 'personsPerGroup' Personen pro Gruppe
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
           mat  <- data.frame ( iter = 1:length(ind), as.numeric(eatTools::crop(substring(input[ind], 13))))
           if(omitUntil>0)  {
              dc<- mat[-c(1:omitUntil),2]                                       ### 'dc' = 'deviance chance'
           } else {
              dc<- mat[,2]
           }
           if ( change ){
              dc<- diff(dc)
              yl<- "Deviance Change"                                            ### labels der y-Achse definieren
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


