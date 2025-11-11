### Hilfsfunktion zur bestimmung von thurstone thresholds fuer GPCM in mirt 
###  Wahrscheinlichkeit, mindestens Kategorie k zu erreichen
### a ist diskriminationsparameter, ds sind die schwellen, m sind die anzahl der schwellen, k ist die konkrete schwelle 
P_get_k <- function(theta, k, m, ds, a) {
    num   <- exp(sapply(k:m, FUN = function(t) {a * sum(theta - ds[1:t])} ) )
    denom <- sum(exp(sapply(0:m, FUN = function(t) {if (t == 0) {return(0)} else {return(a * sum(theta - ds[1:t]))} })))
    return(sum(num) / denom)  }

plotDevianceMirt <- function ( mirt.obj, omitUntil = 1) {
           sysInfo  <- Sys.info()  
           if(sysInfo[["sysname"]] == "Linux") {
              x11(width = 800/72, height = 600/72)
           } else  {
              windows(width = 800/72, height = 600/72)
           }
           dev <- (-1) * diff( (-2) * mirt.obj@Internals$collectLL ) 
           mat <- data.frame ( iter = 1:length(dev), dev = dev, stringsAsFactors = FALSE)
           if(omitUntil>0)  {
              dc<- mat[-c(1:omitUntil),2]                                       ### 'dc' = 'deviance chance'
           } else {
              dc<- mat[,2]
           }
           dc  <- data.frame ( nr=omitUntil + 1:length(dc), dc)
           xm  <- ceiling( max(dc[,1])/10 )*10
           xt   <- NULL
           for ( i in c( 1:30 ) ){xt <- c ( xt, (xm/10) %% i==0 )  }
           xt   <- max ( which ( xt ) )
           cex  <- 0.85 - ( length(dc[,1]) / 1000 )
           if ( cex < 0.40 ) {
                cex <- 0.40
           }
           titel<- "Deviance Change Plot\n"
           par(mar = c(5.1, 4.1, 8.1, 2.1))                                     ### mehr Abstand zwischen titel und plot, damit mehrzeilige ueberschriften reinpassen
           plot ( dc[,1], dc[,2], type="o",
                main=titel,  xlab="Iteration",
                xlim=c(min(dc[,1]),max(dc[,1])),  xaxp=c(0,xm,xt),
                ylab="deviance change", pch=20, cex=cex, lwd=0.75, mar = c(5, 4, 10, 2) + 0.1)
           si   <- devtools::session_info(pkgs = "eatModel")
           si   <- si[["packages"]][which(si[["packages"]][,"package"] == "eatModel"),]
           inf  <- Sys.getenv()
           sysi <- Sys.info()
           sys  <- sessionInfo()
           if(inherits(try(cpu  <- benchmarkme::get_cpu(), silent=TRUE ),"try-error"))  {cpu <- list()}
           if(inherits(try(ram  <- benchmarkme::get_ram(), silent=TRUE ),"try-error"))  {ram <- list()}
           stri <- paste0("'eatModel', version ", si[["loadedversion"]], ", build ",si[["date"]], ", user: ",sysi[["user"]], ", computername: ", ifelse(sysi[["sysname"]] == "Linux", sysi[["nodename"]], inf["COMPUTERNAME"]), "\nsystem: ", sys[["running"]], ", cpu: ", cpu[["model_name"]], ", cores: ",cpu[["no_of_cores"]], ", RAM: ",capture.output(ram))
           stri <- paste0("method = '",mirt.obj@Options[["method"]], "', quadrature points = ",mirt.obj@Options[["quadpts"]], ", final log Likelihood = ",round(mirt.obj@Fit[["logLik"]],digits = 1), ", dims = '",paste(mirt.obj@Model[["factorNames"]], collapse="', '"), "'. Elapsed time: ",timeFormat(mirt.obj@time[[1]],digits = 1, format="s"),"\n", stri, "\n")
           graphics::mtext(stri)
           graphics::abline( a=0, b=0 )
           dcr  <- dc[dc[,2]<0,]
           graphics::points( dcr[,1], dcr[,2], pch=20, cex=cex, col="red")  }

getMirtRegPar <- function(runModelObj, qMatrix) {
   coefs <- coef(runModelObj,  IRTpars = TRUE, printSE = TRUE)[["lr.betas"]]
   if(!is.null(coefs)) {
      for(i in 1:length(coefs)) {colnames(coefs[[i]]) <- colnames(qMatrix)[-1]} ### das jetzt in die Ergebnisstruktur pressen
      mrg <- merge(data.frame(var1 = rownames(coefs[[1]]), coefs[[1]], stringsAsFactors = FALSE), data.frame(var1 = rownames(coefs[[2]]), coefs[[2]], stringsAsFactors = FALSE), by="var1", all=TRUE, suffixes = c("_est", "_se"))
      mrgL<- reshape2::melt(mrg, id.vars = "var1", na.rm=TRUE)
      mrgL<- suppressWarnings(tidyr::separate(mrgL, col = "variable", into = c("group", "derived.par")) |> dplyr::mutate_at(.vars = "derived.par", .funs = car::recode, recodes = "'est'=NA"))
      ret <- data.frame ( par="est", model = attr(runModelObj, "defineModelObj")[["analysis.name"]],source = "mirt", var2 = NA, type = "regcoef", indicator.group = NA,mrgL, stringsAsFactors = FALSE)
      return(ret)
   }}

getMirtPopPar <- function(runModelObj=runModelObj, qMatrix=qMatrix) {
   vals <- mod2values(runModelObj)
   vals <- vals[grep("^COV", vals[,"name"]),]
   dims <- data.frame(do.call("rbind", strsplit(eatTools::halveString(vals[,"name"], "_", first=FALSE)[,2], split = "")), stringsAsFactors = FALSE)
   for(i in 1:ncol(dims)) {dims[,i] <- colnames(qMatrix)[as.numeric(dims[,i]) + 1]}
   colnames(dims) <- paste0("var", 1:ncol(dims))
   dims[which(dims[,2] == dims[,1]),2] <- NA
   dims[,"par"]   <- car::recode(dims[,2], "NA='var'; else = 'correlation'")
   ret  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "mirt", dims, type = "random", indicator.group = NA, group = "persons", derived.par = NA, value = vals[,"value"], stringsAsFactors = FALSE)
   return(ret)}

adaptSkelForAnchor<- function (allNam, skel, anch, qmat, slope, irtmodel, est.slopegroups){
     ### pro dimension checken: wenn es mindestens ein item mit explizit fixierter Trennschaerfe oder eines des Typs 'Rasch'
     ### gibt, muessen die latenten Varianzen frei geschaetzt werden
       dims   <- unlist(lapply(qmat[,-1, drop=FALSE], FUN = function (d) {      ### objekt slope ist fixSlopeMat       
          it  <- qmat[which(d != 0),1]
          it1 <- intersect(it, irtmodel[which(irtmodel[,2]== "Rasch"),1])       ### untere Zeile, hotfix fuer 'it2'! wenn es mehr als zwei spalten im slope data.frame gibt, heissen die 
          it2 <- intersect(it, slope[["ori"]][,allNam[["slopeMatItemCol"]]])    ### nicht zwangslaeufig "item" und "slope" ... bzw. die item-spalte ist nicht zwangslaeufig die erste
          if(length(it1) > 0 || length(it2) > 0 ) {ret <- TRUE} else {ret <- FALSE}
          return(ret)}))
       dimc   <- which(dims == TRUE)
       if(length(dimc)>0) {
          for(g in dimc) {
              ind1 <- grep(paste0("COV_", paste(rep(g,2), collapse="")), skel[,"name"])
              ind2 <- grep("GroupPars", skel[,"class"])
              ind3 <- intersect(ind1, ind2)
              stopifnot(length(ind3)==1)
              skel[ind3,"est"] <- TRUE
          }
       }
     ### standardmaessig alle Kovarianzen frei schaetzen lassen
       ind1   <- grep("^COV", skel[,"name"])
       ind2   <- data.frame(do.call("rbind", strsplit(substr(skel[,"name"], nchar(skel[,"name"])-1, nchar(skel[,"name"])), split="")))
       ind3   <- which(ind2[,1] != ind2[,2])
       ind4   <- intersect(ind1, ind3)
       skel[ind4,"est"] <- TRUE
     ### est.slopegroups: setze die werte in parnum-Spalte gleich, fuer die ein gemeinsamer diskriminationsparameter geschaetzt werden soll 
       if(!is.null(est.slopegroups)){                                           ### est.slopegroups hat immer nur zwei spalten, deshalb sollten hier die spalten numerisch adressiert werden, damit es egal ist, wie der user sie benannt hat
          grps <- by(data= est.slopegroups, INDICES = est.slopegroups[,2], FUN = function(g) {g[,"item"]})
          neu  <- as.numeric(as.factor(names(grps)))
          for(i in 1:length(grps)){
              ind1 <- eatTools::whereAre(grps[[i]], skel[,"item"], verbose=FALSE)
              ind2 <- grep("^a", skel[,"name"])
              ind3 <- which(skel[,"est"] == TRUE)
              ind4 <- intersect(intersect(ind1, ind2), ind3)
              if(length(ind4) >0) {
                 skel[ind4,"const"] <- neu[i]
              }
          }
       }
     ### welche Dimensionen enthalten mindestens 1 verankertes Item? deren latente Mittelwerte muessen frei geschaetzt werden     
       itemsA <- qmat[which(qmat[,"item"] %in% anch[["ank"]][,"item"]),]
       if(nrow(itemsA)>0) {
          toAnk  <- which(colSums(itemsA[,-1, drop=FALSE])>0)
          stopifnot(length(toAnk)>0)
          skel[which(skel[,"name"] %in% paste0("MEAN_", toAnk)),"est"] <- TRUE
       }
     ### wenn es ein hintergrundmodell gibt, muessen die intercepte der dimensionen mit mindestens einem verankerten Item frei geschaetzt werden
       lr     <- grep("^BETA", skel[,"item"]) 
       if(length(lr)>0 && nrow(itemsA)>0) {                                     ### das unter soll ja nur passieren, wenn es hgm UND verankerte items gibt 
          ind1<- grep("(Intercept)", skel[,"name"])
          ind2<- unlist(lapply(names(toAnk), FUN = function (nam) {grep(nam, skel[,"name"])}))
          ind3<- intersect(ind1, ind2)
          skel[ind3,"est"] <- TRUE
       }
     ### jetzt in einer schleife nach und nach alle Ankerparameter ins skeleton eintragen
       for(i in anch[["ank"]][,"item"]) {
          subI <- anch[["ank"]][which(anch[["ank"]][,"item"] == i),]
          subS <- slope[["ori"]][which(slope[["ori"]][,allNam[["slopeMatItemCol"]]] == i),]
          if(is.null(subS) || nrow(subS)==0) {                                  ### in diesem fall: es gibt keinen slope. Das darf fuer Verankerung nur geschehen, wenn Rasch oder PCM
             subIrt <- irtmodel[which(irtmodel[,1] == i),]
             if(subIrt[,2] != "Rasch") {
                stop(paste0("Item '",i,"': You cannot anchoring difficulty parameter without also specifying slope for non-Raschtype items. '",i,"' has type '",unique(subIrt[,2]),"'."))
             } else {
                slp <- 1
             }   
          } else {
             slp <- subS[,allNam[["slopeMatValueCol"]]]
          }  
     ### slope-Parameter in skeleton eintragen und Wert in der "est"-Spalte auf FALSE setzen 
          stopifnot(nrow(skel[intersect(intersect(which(skel[,"item"]==i), grep("^a", skel[,"name"])), which(skel[,"value"] != 0)),])==1)
          skel[intersect(intersect(which(skel[,"item"]==i), grep("^a", skel[,"name"])), which(skel[,"value"] != 0)),"est"] <- FALSE
          skel[intersect(intersect(which(skel[,"item"]==i), grep("^a", skel[,"name"])), which(skel[,"value"] != 0)),"value"] <- slp
          if(unique(skel[which(skel[,"item"] == i),"class"]) == "dich") {       ### dichotome Items: nur einen parameter uebertragen 
             stopifnot(nrow(subI) == 1)
             skel[intersect(which(skel[,"item"] == i),which(skel[,"name"] == "d")), "value"] <- (-1) * subI[,"parameter"] * slp
             skel[intersect(which(skel[,"item"] == i),which(skel[,"name"] == "d")), "est"]   <- FALSE
          } else {                                                              ### polytom
     ### das was in metarep auf Anregung von Janines Mail passiert, muss hier nicht geschehen, da hier nicht einzelnen Stufen
     ### parametrisiert werden und nicht ein "globaler" Itemparameter und die STufen dann als Differenz davon (glaube ich zumindest) 
             for(z in 1:nrow(subI)) {
                 skel[intersect(which(skel[,"item"] == i), which(skel[,"name"] == paste0("d", eatTools::removeNonNumeric(subI[z,"category"])))),"est"] <- FALSE
                 skel[intersect(which(skel[,"item"] == i), which(skel[,"name"] == paste0("d", eatTools::removeNonNumeric(subI[z,"category"])))),"value"] <- (-1) * subI[z,"parameter"] * slp
             }
          }
       }
       return(skel)}


isItPCM <- function(eql) {
       if(inherits(attr(eql[["results"]], "runModelAttributes")[["defineModelObj"]][["irtmodel"]], "data.frame")) {
          pcm <- any(grepl("pcm", attr(eql[["results"]], "runModelAttributes")[["defineModelObj"]][["irtmodel"]][,2], ignore.case=TRUE))
       } else {
          pcm <- FALSE
       }   
       if(isFALSE(pcm)) {
          isPCM <- FALSE                                                        ### das hier zur kompatibilitaet mit aelteren Versionen, wo pcm-Abfrage noch nicht stattfand
       } else {                                                                 ### und ansonsten NA zurueckgegeben wird, was dann zur Fehlermeldung fuehrt
          isPCM <-  pcm || unique(attr(eql[["results"]], "runModelAttributes")[["defineModelObj"]][["irtmodel"]] %in% c("PCM", "PCM2", "GPCM", "GPCM.groups"))
       }
       return(isPCM)}

### Hilfsfunktion fuer defineModel()
prepItemTypeMirt <- function(irtmodel, allNam, qMatrix) {
    if(!all(allNam[["variablen"]] %in% irtmodel[,1])) {                         ### alle items im datensatz auch in der Liste
      cli::cli_abort(c("All items in the data set must be included in the 'irtmodel' data.frame. Missing items:",  "x"=  paste0("'",paste(setdiff(allNam[["variablen"]],irtmodel[,1]), collapse="', '"), "'")))
    }
    if(length(irtmodel[,1]) != length(unique(irtmodel[,1]))) {                  ### liste unique?
      stop("Item identifiers in the first column of 'irtmodel' are not unique.")
    }
    allow<- c("Rasch", "gpcm", "1PL", "2PL", "3PL", "gpcmIRT")                  ### unerlaubte eintraege in der model spalte
    diff <- setdiff(irtmodel[,2] , allow)
    if(length(diff)>0) {
      cli::cli_warn(c("Unexpected entries in the 'type of response mode' column in the 'irtmodel' data.frame. Unexpected entries:",  "i"=paste0("'", paste(diff, collapse="', '"), "'")))
    }
    if(length(allNam[["variablen"]]) != nrow(irtmodel) ||  !all(allNam[["variablen"]] == irtmodel[,1])) {                           
      irtmodel <- irtmodel[which(irtmodel[,1] %in% allNam[["variablen"]]),]     ### nur die items in irtmodel lassen, die es im Datensatz gibt 
      irtmodel <- irtmodel[match(allNam[["variablen"]], irtmodel[,1]),]         ### reihenfolge anpassen
    }
    if("gpcmIRT" %in% irtmodel[,2] && ncol(qMatrix) != 2) {                     ### gpcmIRT nur fuer eindimensionale Modelle 
      stop("'gpcmIRT' is only allowed for uni-dimensional models.")
    }
    return(irtmodel)}   


getMirtInfit <- function(runModelObj, qL, qMatrix) {
   infit<- itemfit(runModelObj, fit_stats = "infit")
   infit<- merge(infit, qL[,-match("value", colnames(qL))],  by.x = "item", by.y = colnames(qMatrix)[1], all = TRUE)
   ret  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "mirt", var1 = infit[,"item"], var2 = NA , type = "fixed", indicator.group = "items",   group = infit[,"dimensionName"], par = "est",  derived.par = "infit", value = infit[,"infit"], stringsAsFactors = FALSE)
   return(ret) }

getMirtResults <- function(runModelObj, omitFit, omitRegr, omitWle, omitPV) {
         qMatrix<- attr(runModelObj, "defineModelObj")[["qMatrix"]]
         qL     <- eatTools::facToChar(reshape2::melt(qMatrix, id.vars = colnames(qMatrix)[1], variable.name = "dimensionName", na.rm=TRUE))
         qL     <- qL[which(qL[,"value"] != 0 ) , ]
         varName<- colnames(qMatrix)[1]
    ### deskriptive Werte auslesen
         ret    <- getTamDescriptives(runModelObj=runModelObj, qL=qL, qMatrix=qMatrix, leseAlles = TRUE, software = "mirt")
    ### (saemtliche) itemkennwerte auslesen 
         ret    <- rbind(ret, getMirtItempars(runModelObj=runModelObj, qL=qL, qMatrix=qMatrix, leseAlles = TRUE, software = "mirt"))
    ### Diskriminationswerte auslesen: hier kann dieselbe Funktion verwendet werden wie fuer TAM
         ret    <- rbind(ret, getTamDiscrim(runModelObj=runModelObj, qL=qL, qMatrix = qMatrix, leseAlles = TRUE, software="mirt"))
    ### Infit auslesen
         ret    <- rbind(ret, getMirtInfit(runModelObj=runModelObj, qL=qL, qMatrix = qMatrix))
    ### Populationsparameter auslesen (tbd)
         ret    <- rbind(ret, getMirtPopPar(runModelObj=runModelObj, qMatrix=qMatrix))
    ### Regressionsparameter auslesen
         ret    <- rbind(ret, getMirtRegPar(runModelObj=runModelObj, qMatrix=qMatrix))
    ### Modellindizes auslesen
         # ret    <- rbind(ret, getTamModInd(runModelObj=runModelObj, leseAlles = leseAlles))
    ### Personenparameter auslesen (WLEs)
         beg    <- Sys.time()
         ret    <- rbind(ret, getMirtWles(runModelObj=runModelObj, qMatrix=qMatrix,  omitWle = omitWle))
         diffe <- Sys.time() - beg
         if(as.numeric(diffe) > 0.2) {message(paste0("Getting WLEs calling fscores(method=\"WLE\") from getMirtWles: ", timeFormat(diffe)))}
    ### PVs auslesen
         beg    <- Sys.time()
         ret    <- rbind(ret, getMirtPVs ( runModelObj=runModelObj, qMatrix=qMatrix, omitPV = omitPV))
         diffe <- Sys.time() - beg
         if(as.numeric(diffe) > 0.2) {message(paste0("Getting PVs calling fscores from getMirtPVs: ", timeFormat(diffe)))}
    ### EAPs auslesen
         ret    <- rbind(ret, getMirtEAPs(runModelObj=runModelObj, qMatrix=qMatrix))
    ### Q3 auslesen
# weiter ab hier
         #beg    <- Sys.time()
         #ret    <- rbind(ret, getTamQ3(runModelObj=runModelObj, leseAlles = leseAlles, shw1 = resItem[["shw1"]], Q3=Q3, q3MinObs=q3MinObs, q3MinType=q3MinType))
         #diffe  <- Sys.time() - beg
         #if(as.numeric(diffe) > 0.2) {message(paste0("Getting Q3 statistic calling tam.modelfit from getTamQ3: ", timeFormat(diffe)))}
         return(ret)}         

getMirtItempars <- function(runModelObj, qL, qMatrix, leseAlles, software) {
      coefs <- coef(runModelObj, IRTpars = TRUE, printSE = TRUE) 
      coef2 <- coef(runModelObj, IRTpars = TRUE, simplify = TRUE)$items         ### notwendig fuer thurstonian thresholds 
      items <- attr(runModelObj, "defineModelObj")[["allNam"]][["variablen"]]
      ret   <- do.call("rbind", lapply(items, FUN = function (i) {              ### untere Zeile: bei mehreren dimensionen gibt es soviele slope-parameter pro item wie es dimensionen gibt
               rowEst <- which(rownames(coefs[[i]])=="par")                     ### bei between item dimensionality modellen will ich aber nur den auswaehlen, der zu der dimension gehoert,
               rowSE  <- which(rownames(coefs[[i]])=="SE")                      ### auf der das Item laedt
               colSlo <- grep("^a", colnames(coefs[[i]]), value=TRUE)[which(subset(qMatrix, item == i)==1)-1]
               if(length(colSlo) != 1) {browser()}                              ### das sollte nie passieren
               if("b" %in% colnames(coefs[[i]])) {                              ### dichotomer Fall (Rasch und 2pl)
                  est    <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "mirt",  var1 = i, var2 = NA , type = "fixed", indicator.group = "items", group = subset(qL, item==i)[,"dimensionName"],  par = car::recode(coefs[[i]][rowSE,"b"],"NA='offset';else='est'"),  derived.par = NA, value = coefs[[i]][rowEst,"b"], stringsAsFactors = FALSE)
                  se     <- est |> dplyr::mutate(derived.par = "se", value = coefs[[i]][rowSE,"b"])
               } else {                                                         ### polytomer Fall (pcm und gpcm)
                  se     <- NULL                                                ### initialisieren, damit es am ende ge-rbinded werden kann
                  cols   <- grep("^b", colnames(coefs[[i]]), value=TRUE)
                  est    <- do.call("rbind", lapply(cols, FUN = function(col1) {
                            est1 <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "mirt",  var1 = i,  var2 = paste0("Cat", eatTools::removeNonNumeric(col1))  , type = "fixed", indicator.group = "items", group = subset(qL, item==i)[,"dimensionName"],  par = car::recode(coefs[[i]][rowSE,col1],"NA='offset'; else='est'"),  derived.par = NA, value = coefs[[i]][rowEst,col1], stringsAsFactors = FALSE)
                            se1  <- est1 |> dplyr::mutate(derived.par = "se", value = coefs[[i]][rowSE,col1])
                            return(rbind(est1, se1))}))
               }
               slope  <- est[1,] |> dplyr::mutate(par = "estSlope", value = coefs[[i]][rowEst,colSlo])
               slopeSE<- est[1,] |> dplyr::mutate(par = "estSlope", derived.par = "se",value = coefs[[i]][rowSE,colSlo]) |> eatTools::na_omit_selection(varsToOmitIfNA = "value")
    ### Thurstonian thresholds ausgeben lassen
    ### im dichotomen modell werden thresholds analog zur itemparametertransformaton bestimmt: wegen der abweichenden mirt 
    ### parametrisierung kann man nicht dieselbe Funktion nehmen, obwohl es konzeptuell aequivalent ist 
               if("b" %in% colnames(coefs[[i]])) {
                  thurs <- coefs[[i]][rowEst,"b"] + log(0.625/(1-0.625))
               } else {
                  pars  <- coef2[which(rownames(coef2) == i), , drop=TRUE]
                  ds    <- pars[grep("^(d|b)\\d+", names(pars))] |> as.numeric()
                  ds    <- ds[!is.na(ds)]
                  bounds<- c(-6, 6)                                             
                  while(isTRUE(inherits(try(thurs  <- sapply(1:length(ds), FUN = function(k) { uniroot( function(th) P_get_k(theta = th, k=k, m = length(ds), ds=ds, a = coefs[[i]][rowEst,colSlo]) - 0.625, interval = bounds)$root  }) , silent=TRUE ),"try-error"))) {
                     bounds[1] <- bounds[1] - 1                                 ### theta bounds vergroessern wenn itemparameter ausserhalb des ranges und es daher fehlermeldung gibt 
                     bounds[2] <- bounds[2] + 1
                     cat(paste0("Extend theta bounds for item '",i,"' to ", bounds[1], ", ",  bounds[2], "\n"))
                  }
               }
               thurs  <- data.frame(model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "mirt",  var1 = i,var2 = paste0("Cat", 1:length(thurs)), type="fixed", indicator.group = "items", group = subset(qL, item==i)[,"dimensionName"],  par = "est",  derived.par = "thurstone", value = thurs, stringsAsFactors = FALSE)
               ret    <- rbind(est, se, slope, slopeSE, thurs)
               return(ret)}))
      ismis <- which(is.na(ret[,"var2"]))
      if(length(ismis)>0) {ret[ismis,"var2"] <- "Cat1"}
      return(ret)}         
      
      
getMirtWles <- function(runModelObj, qMatrix, omitWle) {
         if(omitWle == TRUE ) {return(NULL)}
         if(ncol(qMatrix) ==2) {                                                ### fuer eindimensionale MOdelle funktion direkt aufrufen
            wle  <- fscores(runModelObj,method="WLE", verbose=FALSE, full.scores.SE = TRUE)
         } else {                                                               ### fuer mehrdimensionale Modelle ueber batch-aufruf, um die konsole nicht mit verwirrenden mirt fehlermeldungen zu bombardieren
            wle  <- mirtWlesMultidim(runModelObj)
         }
         ids  <- attr(runModelObj, "personID")
         wleL <- reshape2::melt(data.frame(ID = ids, wle,stringsAsFactors = FALSE), id.vars = "ID",  na.rm=TRUE)
         dims <- colnames(qMatrix)[-1]
         wleL[,"group"] <- dims[as.numeric(eatTools::removeNonNumeric(stringr::str_remove(as.character(wleL[,"variable"]), pattern = "^SE_")))]
         wleL[,"derived.par"] <- car::recode(eatTools::halveString(as.character(wleL[,"variable"]), "_")[,2], "NA='est'; else = 'se'")
         res  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "mirt", var1 = wleL[,"ID"], var2 = NA,  type = "indicator", indicator.group = "persons",  group = wleL[,"group"], par = "wle", derived.par = wleL[,"derived.par"], value = wleL[,"value"], stringsAsFactors = FALSE)
         return(res)}           

mirtWlesMultidim <- function(m) {
    sysInfo  <- Sys.info()
    nams     <- sample(1:99999, size = 3, replace=FALSE)
    filenam1 <- paste0("F", nams[1], ".rds")
    filenam2 <- paste0("W", nams[2], ".rds")
    filenam3 <- paste0("S", nams[3], ".R")
    saveRDS(m, file = file.path(tempdir(), filenam1), compress="bzip2")
    syn      <- c("library(mirt)", paste0("m <- readRDS(\"",normalizePath(file.path(tempdir(), filenam1), winslash = "/"),"\")"),  "wle  <- fscores(m,method=\"WLE\", verbose=FALSE, full.scores.SE = TRUE)",   paste0("saveRDS(wle, file = \"",normalizePath(file.path(tempdir(), filenam2), winslash = "/"),"\")")  )
    write(syn, file = file.path(tempdir(), filenam3), sep="\n")
    setwd(tempdir())
    if(sysInfo[["sysname"]] == "Linux") {
       cat("#!/bin/bash\n",paste0( "R CMD BATCH --vanilla ",filenam3," syntax.Rout\n"), file = "run_analysis.sh")
       Sys.chmod("run_analysis.sh", mode = "0755")                              ### macht die Datei ausfuehrbar
       system("./run_analysis.sh")
    } else {
       rcmd<- file.path(Sys.getenv()[["R_HOME"]], "bin/x64/Rcmd.exe")           ### Pfad der 'rcmd' angeben
       pfad<- tempdir()
       bat <- c(substr(pfad, 1, 2), paste0("cd ", normalizePath (pfad), "\\"), paste0("CALL ", "\"",normalizePath(rcmd),"\" BATCH --vanilla ",filenam3," syntax.rout"), "exit")
       write(bat,file.path(pfad, "start.bat"), sep="\n")
       system(file.path(pfad, "start.bat"), intern=FALSE, show.output.on.console = FALSE, wait=TRUE, invisible = FALSE)
    }
    wle <- readRDS(file.path(tempdir(), filenam2))
    return(wle)} 

getMirtPVs <- function ( runModelObj, qMatrix, omitPV) {
         if(omitPV == TRUE ) {return(NULL)}
         pv   <- fscores(runModelObj, plausible.draws = attr(runModelObj, "defineModelObj")[["n.plausible"]])
         id   <- attr(runModelObj, "personID")
         dims <- colnames(qMatrix)[-1]
         pvL  <- do.call("rbind", lapply(1:length(pv), FUN = function(imp) {
                 long <- reshape2::melt(data.frame(ID = id, pv[[imp]], stringsAsFactors = FALSE), id.vars = "ID") |> dplyr::mutate(derived.par = paste0("pv", imp))
                 long[,"group"] <- dims[as.numeric(long[,"variable"])]
                 return(long)}))
         res  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "mirt", var1 = pvL[,"ID"], var2 = NA , type = "indicator", indicator.group = "persons", par = "pv", pvL[,c("group", "derived.par", "value")], stringsAsFactors = FALSE)
         return(res)}

getMirtEAPs <- function ( runModelObj, qMatrix) {
         eap  <- fscores(runModelObj,method = "EAP", full.scores.SE=TRUE)
         id   <- attr(runModelObj, "personID")
         dims <- colnames(qMatrix)[-1]
         eapL <- reshape2::melt(data.frame(ID = id, eap,stringsAsFactors = FALSE), id.vars = "ID",  na.rm=TRUE)
         eapL[,"group"] <- dims[as.numeric(eatTools::removeNonNumeric(stringr::str_remove(as.character(eapL[,"variable"]), pattern = "^SE_")))]
         eapL[,"derived.par"] <- car::recode(eatTools::halveString(as.character(eapL[,"variable"]), "_")[,2], "NA='est'; else = 'se'")
         res  <- data.frame ( model = attr(runModelObj, "defineModelObj")[["analysis.name"]], source = "mirt", var1 = eapL[,"ID"], var2 = NA, type = "indicator", indicator.group = "persons", par = "eap", eapL[,c("group", "derived.par", "value")],stringsAsFactors = FALSE)
         return(res)}

