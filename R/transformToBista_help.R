### called by transformToBista()

adaptEatRepVersion <- function ( x ) {
  if ( inherits(x, "data.frame"))  {
    return ( x )
  }  else  {
    x <- x[[1]][[1]]
    stopifnot ( inherits(x, "data.frame") )
    return(x)
  } }

### ----------------------------------------------------------------------------

createLinkingErrorObject <- function (itempars, years) {
  if(is.null(itempars) || nrow(itempars) == 0 || !"dimension" %in% colnames(itempars)) {return(NULL)}
  depVars <- character()
  if("linkingError" %in% colnames(itempars)) {depVars <- c(depVars, "value")}
  if("linkingErrorTransfBista" %in% colnames(itempars)) {depVars <- c(depVars, "valueTransfBista")}
  if(all(c("traitLevel", "linkingErrorTraitLevel") %in% colnames(itempars))) {depVars <- c(depVars, "traitLevel")}
  if(length(depVars) == 0) {return(NULL)}                                       ### die Zeilen bis hierhin kommen von Karolines KI
  res <- do.call("rbind", by(data = itempars, INDICES = itempars[,"dimension"], FUN = function (d) {
    r1 <- do.call("rbind", lapply(depVars, FUN = function (av) {
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

### ----------------------------------------------------------------------------

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

### called by createItemVeraObj() ----------------------------------------------

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

checkUserCuts <- function(itFrame, cuts){
       if(!itFrame[1,"dimension"] %in% names(cuts) ) {stop(paste("Cannot found dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.",sep=""))}
       mat1<- match( unique(itFrame[,"dimension"]), names(cuts))
       if(!"values" %in% names(cuts[[mat1]]) ) {stop(paste("'cuts' must be a named list. Cannot found 'values' element for dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
       if(length(cuts[[mat1]])>1) {
          if ( !"labels" %in% names(cuts[[mat1]]) ) {stop(paste("'cuts' must be a named list. Cannot found 'labels' element for dimension '",itFrame[1,"dimension"],"' in the 'cuts' list.\n",sep=""))}
       }
       return(mat1)}

transform625Dichotom <- function(itFrame, equatingList, resMD) {
       slp1<- grep("slope", colnames(itFrame), value=TRUE, ignore.case=TRUE)
       slp1<- slp1[!grepl("(^se|_se$|\\.se$|se$|stderr|std\\.err|error)", slp1, ignore.case=TRUE)]
       slp2<- grep("est", colnames(itFrame), value=TRUE, ignore.case=TRUE)
       if(length(slp1) < 2 || length(slp2) ==0) {                               ### wenn es nur eine spalte mit "slope" im Namen gibt, ist das die estimator-spalte,
          slp <- slp1                                                           ### selbst wenn sie nicht mit "est" benannt ist. Gibt es zusaetzlich noch einen
       } else {                                                                 ### Standardfehler des slopes, gibt es mehrere Spalten mit "slope" im Namen; in diesem
          slp <- intersect(slp1, slp2)                                          ### Falle soll die ausgewaehlt werden, die zusaetzlich "est" im Namen traegt
       }
       stopifnot(length(slp) %in% 0:1)                                          ### es darf nur eine oder gar keine slope Spalte geben
       if(length(slp)==0) {
          itFrame[,"estTransf625"]<- itFrame[,"estTransf"] + log(0.625/(1-0.625))## alte Variante (1pl)
       } else {
          if(any(itFrame[,slp] < 0)) {
             warning("The slope parameters for some items are less than 0. These item parameters cannot be meaningfully transformed into the educational standards metric.")
          }                                                                     ### untere Zeile: Variante fuer 2pl entsprechend Karolines AI-Agent (Mail Karoline, 11.06.2026, 9.24 Uhr, bzw. "https://github.com/weirichs/eatModel/issues/34#issuecomment-4829725821")
          software <- attr(equatingList[["results"]], "runModelAttributes")[["software"]]
          if(is.null(software) || length(software) == 0) {
             if(!is.null(resMD) && "source" %in% colnames(resMD)) {
                software <- unique(na.omit(resMD[,"source"]))
             } else {
                if("source" %in% colnames(equatingList[["results"]])) {
                   software <- unique(na.omit(equatingList[["results"]][,"source"]))
                }
             }
          }
          if(length(software) != 1) {stop("Cannot uniquely identify estimation software for 2PL item parameter transformation.")}
          if(!software %in% c("tam", "mirt")) {stop("'software' must be either TAM or mirt.")}
          if(software == "tam") {
             itFrame[,"estTransf625"]<- (itFrame[,"estTransf"] + log(0.625/(1-0.625))) / itFrame[,slp]
          } else {                                                              ### obere Zeile: tam; untere Zeile: mirt
             itFrame[,"estTransf625"]<- itFrame[,"estTransf"] + log(0.625/(1-0.625)) / itFrame[,slp]
          }
       }
       return(itFrame)}

computeMeanSD_focPop <- function(weights, pv, id, dimname) {
       if(is.null(weights) ) {
          msdF<- eatRep::repMean ( datL = pv, ID = id, imp = "imp", dependent = "valueTransfBista", na.rm = TRUE, verbose = FALSE, progress = FALSE)
          msdF<- adaptEatRepVersion(msdF)
       } else {
          pvF <- eatTools::mergeAttr ( pv, weights , by.x = id, by.y = colnames(weights)[1], all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "cases", xName = paste0("plausible values for dimension ",dimname), yName = "weights", verbose = c("match", "dataframe"))
          mis <- which(is.na(pvF[,colnames(weights)[2]]))
          if(length(mis) > 0 ) {                                                ### missings in the weights frame are not allowed
             cat(paste ( "Found ",length(mis)," missing values in the 'weights' frame.\n    Cases with missing values on weighting variable will be ignored for transformation.\n",sep=""))
             pvF <- pvF[-mis,]
          }
          msdF<- eatRep::repMean ( datL = pvF, ID = id, imp = "imp", wgt = colnames(weights)[2], dependent = "valueTransfBista", na.rm = TRUE, verbose = FALSE, progress = FALSE)
          msdF<- adaptEatRepVersion(msdF)
       }
       msdFok <- c(msdF[intersect(which(msdF[,"parameter"] == "mean"), which(msdF[,"coefficient"] == "est")),"value"], msdF[intersect(which(msdF[,"parameter"] == "sd"), which(msdF[,"coefficient"] == "est")),"value"])
       return(msdFok)}

computeTraitLevelSEs <- function(itFrame, cutsMis, equatingList, isPCM, cuts, mat1, mat, msdFok, traitLevel) {
       le      <- NULL                                                          ### initialisieren
       if(cutsMis == FALSE & !is.null ( equatingList[["items"]] ) & !isPCM) {
          cts <- c( -10^6, cuts[[mat1]][["values"]], 10^6)
          le  <- do.call("rbind", lapply ( (length(cts)-1):1 , FUN = function ( l ) {
                 kmp<- c(cts[l], cts[l+1])       ### Linkingfehler fuer einzelnen Kompetenzintervalle; absteigend wie bei karoline
                 a1 <- sum ( dnorm ( ( kmp - mat[,5]) / mat[,6] ) * c(-1,1) / mat[,6] )
                 a2 <- sum ( dnorm ( ( kmp - msdFok[1]) / msdFok[2] ) * c(-1,1) / msdFok[2] )
    ### Achtung! der 'mutmassliche Fehler' kann auch auftreten, wenn das Equating zuvor durchgeschleift wurde und deshalb gar keine Linkingfehler berechnet werden koennen
                 if(a2 == 0 ) {cat("Error during the transformation of linking errors. If the equating was previously looped through - meaning no equating took place in the strict sense - this message is non-critical and can be ignored.\n")}
                 del<- ( (  a1^2 + a2^2 ) * (unique(itFrame[,"linkingErrorTransfBista"])^2) / 2  )^0.5
                 del<- data.frame ( traitLevel = attr(traitLevel, "cat.values")[l], linkingErrorTraitLevel = del )
                 return(del)}))
    ### ggf. weg! 'linkingErrorTraitLevel' ergibt ja fuer Items keinen Sinn, nur fuer Personenparameter
          ori <- colnames(itFrame)
          chk <- unique(le[,"traitLevel"]) %in% unique(itFrame[,"traitLevel"])
          if(length( which(chk == FALSE)) > 0) {
             warning(paste("Model '",unique(itFrame[,"model"]),"', dimension '",unique(itFrame[,"dimension"]),"': No items on trait level(s) '",paste( unique(le[,"traitLevel"])[which(chk == FALSE)], collapse = "', '"), "'.", sep=""))
          }
          itFrame <- eatTools::mergeAttr ( itFrame, le, by = "traitLevel", sort = FALSE, all.x = TRUE, all.y = FALSE,  setAttr = FALSE, unitName = "trait levels", xName = "item parameter list", yName = "linking error list", verbose = "match")
          itFrame <- itFrame[,c(ori, "linkingErrorTraitLevel")]
       }
       itFrame <- itFrame |> dplyr::mutate(refMean= mat[,3], refSD = mat[,4], refTransfMean=mat[,5], refTransfSD= mat[,6])
       return(list(itFrame=itFrame, le=le))}

transformPersonParameter <- function(cutsMis, pv, cuts, mat1, isPCM, equatingList, le, itFrame, mod, dimname, weights, id) {
       if(isFALSE(cutsMis) ) {pv[,"traitLevel"]   <- eatTools::num.to.cat(x = pv[,"valueTransfBista"], cut.points = cuts[[mat1]][["values"]], cat.values = cuts[[mat1]][["labels"]])}
       pv[,"dimension"]  <- pv[,"group"]
       if(isFALSE(cutsMis) && !isPCM && !is.null ( equatingList[["items"]]) ) {
          if(!is.null(le)) {
             warning("Skip check whether all competence levels are occupied (due to bayesian plausible values imputation).")
          } else {
             chk <- unique(le[,"traitLevel"]) %in% unique(pv[,"traitLevel"])
             if(length( which(chk == FALSE)) > 0) {warning(paste("Model '",unique(itFrame[,"model"]),"', dimension '",unique(itFrame[,"dimension"]),"': No plausible values on trait level(s) '",paste( unique(le[,"traitLevel"])[which(chk == FALSE)], collapse = "', '"), "'.", sep=""))}
          }
       }
       if(!isPCM && !is.null ( equatingList[["items"]]) && !is.null(itFrame) && "linkingErrorTransfBista" %in% colnames(itFrame)) {
          leTransf <- unique(na.omit(itFrame[,"linkingErrorTransfBista"]))
          stopifnot ( length(leTransf) %in% 0:1)
          pv[,"linkingError"] <- equatingList[["items"]][[mod]][[dimname]][["eq"]][["descriptives"]][["linkerror"]]
          pv[,"linkingErrorTransfBista"] <- if(length(leTransf)==0) {NA_real_} else {leTransf}
       }
       ori <- colnames(pv)                             ### nur wenn untere Bedingung == TRUE, gibt es das Objekt 'le', das gemergt werden soll
       if(cutsMis == FALSE && !is.null ( equatingList[["items"]]) && exists("le", inherits = FALSE) ) {
          pv  <- eatTools::mergeAttr ( pv, le, by = "traitLevel", sort = FALSE, all.x = TRUE, all.y = FALSE, setAttr = FALSE, unitName = "trait levels", xName = "plausible values", yName = "linking error list", verbose = "match")
          pv  <- pv[,c(ori, "linkingErrorTraitLevel")]
       }
    ### ggf. Gewichte an Personenframe mit dranhaengen
       if(!is.null(weights)) {
          pv  <- merge ( pv, weights , by.x = id, by.y = colnames(weights)[1], all.x = TRUE, all.y = FALSE)
       }
       return(pv)}
