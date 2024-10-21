### called by different functions like getResults(), plotICC(),
### transformToBista() or equat1pl()

### called by getRestuls() and plotICC() ---------------------------------------

eapFromRes <- function (resultsObj, idVarName = NULL, verbose = TRUE){
  checkmate::assert_data_frame(resultsObj)
  checkmate::assert_character(idVarName, null.ok = TRUE)
  checkmate::assert_logical(verbose, len = 1)
  #
  eapRo <- setdiff(intersect(which(resultsObj[,"par"] == "eap"),
                             which(resultsObj[,"indicator.group"] == "persons")),
                   which(resultsObj[,"derived.par"] == "rel"))
  id <- unique(resultsObj[intersect(which(resultsObj[,"type"] == "tech"),
                                    which(resultsObj[,"par"] == "ID")),"derived.par"])
  id <- getIdVarName(id, idVarName, verbose=verbose)
  if(length(eapRo) == 0){
    warning("'resultsObj' does not contain any eap values.")
    return(NULL)
  }  else  {
    sel  <- resultsObj[eapRo,]
    sel  <- do.call("rbind", by(sel, INDICES = sel[,c("model", "group")], FUN = function ( gr ) {
      res  <- reshape2::dcast ( gr , model+group+var1~derived.par, value.var = "value")
      colnames(res)[-c(1:2)] <- c(id, "EAP", "SE.EAP")
      weg  <- match(c("model", id), colnames(res))
      res  <- data.frame(res[,c("model", id)], dimension = as.character(gr[1,"group"]),
                         res[,-weg,drop=FALSE], stringsAsFactors = FALSE)
      return(res)}))
    return(sel)
  }
}

### called by getRestuls() and plotICC() ---------------------------------------

wleFromRes <- function(resultsObj, idVarName = NULL, verbose=TRUE) {
  checkmate::assert_data_frame(resultsObj)
  checkmate::assert_character(idVarName, null.ok = TRUE)
  checkmate::assert_logical(verbose, len = 1)
  #
  wleRo <- setdiff(intersect(which(resultsObj[,"par"] %in% c("wle","NitemsSolved", "NitemsTotal")),
                             which(resultsObj[,"indicator.group"] == "persons")),
                   which(resultsObj[,"derived.par"] == "rel"))
  if(length(wleRo) == 0){
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

### called by getRestuls(), transformToBista() and plotICC() -------------------

pvFromRes <- function(resultsObj, toWideFormat = TRUE, idVarName = NULL, verbose=TRUE) {
  checkmate::assert_data_frame(resultsObj)
  checkmate::assert_character(idVarName, null.ok = TRUE)
  lapply(c(toWideFormat, verbose), checkmate::assert_logical, len = 1)
  #
  pvRow <- intersect(which(resultsObj[,"par"] == "pv"),
                     which(resultsObj[,"indicator.group"] == "persons"))
  if(length(pvRow) == 0){
    warning("'resultsObj' does not contain any pv values.")
    return(NULL)
  }  else  {
    sel  <- resultsObj[pvRow, ]
    id   <- unique(resultsObj[intersect(which(resultsObj[,"type"] == "tech"),
                                        which(resultsObj[,"par"] == "ID")),"derived.par"])
    id   <- getIdVarName(id, idVarName, verbose=verbose)
    if (toWideFormat == TRUE ) {
      sel  <- do.call("rbind", by(sel, INDICES = sel[,c("model","group")], FUN = function(gr){
        res  <- reshape2::dcast(gr, model+var1~derived.par, value.var = "value")
        colnames(res)[2] <- id
        weg  <- match(c("model", id), colnames(res))
        res  <- data.frame(res[,c("model", id)], dimension = as.character(gr[1,"group"]),
                           res[,-weg,drop=FALSE], stringsAsFactors = FALSE)
        return(res)}))
    }  else  {
      sel  <- sel[,c("model", "var1", "group", "derived.par", "value")]
      recSt<- paste ( "'var1'='",id,"'; 'derived.par'='imp'",sep="")
      colnames(sel) <- car::recode ( colnames(sel), recSt)
    }
    return(sel)
  }  }

### called by getRestuls(), equat1pl(), transformToBista() and prepareAndCheckEatModelObject()

itemFromRes<- function(resultsObj){
  checkmate::assert_data_frame(resultsObj)
  #
  res <- do.call(plyr::rbind.fill, by(data = resultsObj, INDICES = resultsObj[,"model"],
                                      FUN = function(mod){
    sel <- mod[intersect(which(mod[,"par"] %in% c("est", "estSlope", "Nvalid", "itemP", "ptBis", "itemDiscrim", "offset")),
                         which(mod[,"indicator.group"] == "items")),]
    if (nrow(sel)==0) {
      return(NULL)
    }  else  {
      isDif <- intersect(which(mod[,"type"] == "tech"), which(mod[,"par"] == "DIF.var"))
      if(length(isDif) > 0) {
        vars     <- mod[intersect(which(mod[,"type"] == "tech"),
                                  which(mod[,"par"] == "variablen")),"derived.par"]
        itemList <- do.call("rbind", lapply ( vars, FUN = function ( v ) {
          ind <- grep( paste0("_",v,"_"), sel[,"var1"])
          it  <- sort ( unique ( sel[ind,"var1"]))
          if(length(it)>2) {
            warning(paste0("DIF variable '",mod[isDif,"derived.par"],
                           "' seems to have more than two categories. To date, this is not supported by 'eatModel'."))
          }
          return ( data.frame ( item = v, dif = it[1], weg = it[length(it)] , stringsAsFactors = FALSE) ) }))
        weg      <- eatTools::whereAre ( itemList[,"weg"], sel[,"var1"], verbose=FALSE)
        forDif   <- eatTools::whereAre ( itemList[,"dif"], sel[,"var1"], verbose=FALSE)
        stopifnot(length( intersect(weg, forDif)) == 0 )
        selForDif<- sel[forDif, ]
        sel      <- sel[-c(weg, forDif) , ]
        sel      <- sel[which ( sel[,"par"] != "ptBis" ) , ]
        selDIF   <- do.call("rbind", by(selForDif, INDICES = selForDif[,"group"], FUN = function(gr){
          res  <- reshape2::dcast ( gr , model+var1~par+derived.par, value.var = "value")
          mat  <- lapply( vars, FUN = function ( v ) { grep(paste0("_",v,"_"), res[,"var1"])})
          stopifnot (  all ( sapply(mat, length) == 1) )
          res[unlist(mat),"item"]  <- vars
          colnames(res) <- car::recode ( colnames(res) , "'est_infit'='infitDif'; 'est_se'='seDif'; 'est_NA'='estDif'")
          res[,"absDif"]<- abs ( res[,"estDif"]  * 2 )
          pval <- intersect(intersect(which(mod[,"type"] == "tech"), which(mod[,"par"] == "dif")),
                            which(mod[,"derived.par"] == "p.value"))
          stopifnot (length(pval) == 1)
          pval <- mod[pval, "value"]
          adb  <- mod[intersect(intersect(which(mod[,"type"] == "tech"),
                                          which(mod[,"par"] == "dif")),
                                which(mod[,"derived.par"] == "abs.dif.bound")),"value"]
          sdb  <- mod[intersect(intersect(which(mod[,"type"] == "tech"),
                                          which(mod[,"par"] == "dif")),
                                which(mod[,"derived.par"] == "sig.dif.bound")),"value"]
          res[,paste("CI__", pval ,"__lb",sep="")] <- res[,"absDif"] - 2*abs(qnorm(0.5*(1-pval))) * res[,"seDif"]
          res[,paste("CI__", pval ,"__ub",sep="")] <- res[,"absDif"] + 2*abs(qnorm(0.5*(1-pval))) * res[,"seDif"]
          res  <- data.frame ( res, do.call("rbind", apply(res[,c("absDif", "seDif",
                                                                  paste("CI__",pval,"__lb",sep=""),
                                                                  paste("CI__",pval,"__ub",sep=""))],
                                                           MARGIN = 1, FUN = function(d){
            check <- all ( !is.na(d) )
            if(check == TRUE) {
              crit1 <- d[["absDif"]] > adb
              crit2 <- !all(sort(c(d[[paste("CI__",pval,"__lb",sep="")]], sdb,
                                   d[[paste("CI__",pval,"__ub",sep="")]]), index.return = TRUE)$ix == 1:3)
              if ( crit1 == TRUE & crit2 == TRUE) { res <- 1 }  else { res <- 0}
              ets   <- "A"
              ets1  <- d[["absDif"]] > 0.43
              ets2  <- !all(sort(c(d[[paste("CI__",pval,"__lb",sep="")]], 0 ,
                                   d[[paste("CI__",pval,"__ub",sep="")]]), index.return = TRUE)$ix == 1:3 )
              if ( ets1 == TRUE & ets2 == TRUE) { ets <- "B" }
              etsC1 <- d[["absDif"]] > 0.64
              etsC2 <- !all(sort(c(d[[paste("CI__",pval,"__lb",sep="")]], 0.43,
                                   d[[paste("CI__",pval,"__ub",sep="")]]), index.return = TRUE)$ix == 1:3 )
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
        colnames(res) <- car::recode(colnames(res),"'var1'='item'; 'est_infit'='infit'; 'est_outfit'='outfit'; 'est_se'='se'; 'est_NA'='est'; 'estSlope_se'='seSlope'; 'estSlope_NA'='estSlope'; 'offset_NA'='estOffset'; 'Nvalid_NA'='Nvalid'; 'ptBis_NA'='ptBis'; 'itemP_NA'='itemP'; 'itemDiscrim_NA'='itemDiscrim'")
        cols <- c("Nvalid", "itemP", "itemDiscrim", "est", "estOffset", "se",
                  "estSlope", "seSlope", "infit","outfit", "ptBis")
        drin1<- which(cols %in% colnames(res))
        drin2<- grep("ptBis_", colnames(res))
        drin3<- grep("itemP", colnames(res))
        res  <- data.frame(res[,c("model", "item")], dimension = as.character(gr[1,"group"]),
                           res[,c(cols[drin1], colnames(res)[drin2],
                                  setdiff(colnames(res)[drin3], cols[drin1])),drop=FALSE],
                           stringsAsFactors = FALSE)
        return(res)}))
      if ( length( isDif ) > 0 ) {
        ciCo<- colnames(selDIF)[grep("^CI__", colnames(selDIF))]
        sel <- merge(sel, selDIF[,c("item", "model", "estDif", "seDif", "infitDif",
                                    "absDif", ciCo, "difIndex", "ETS")],
                     by=c("item","model"), all=TRUE)
      }
      return(sel)
    }
  }))
  return(res)}

### called by getRestuls() and addQ3() -----------------------------------------

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

### called by transformToBista() and other fromRes-functions  ------------------

getIdVarName <- function ( id, idVarName, verbose=TRUE) {
  if (length( id ) == 0 ) {
    if ( is.null(idVarName)) { new <- "idstud"} else { new <- idVarName}
    if(verbose){warning(paste0("Cannot identify student identifier variable (possibly because 'resultsObj' was created by an older version of 'eatModel'). student id variable will be defaulted to '",new,"'."))}
    id <- new
  }
  return(id)}

### not called  ----------------------------------------------------------------

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

### not called  ----------------------------------------------------------------

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

### not called  ----------------------------------------------------------------

wleRelFromRes <- function(resultsObj) {
  ret <- resultsObj[intersect(which(resultsObj[,"derived.par"] == "rel"), which(resultsObj[,"par"] == "wle")),c("model", "group", "value")]
  colnames(ret) <- car::recode(colnames(ret), "'value'='rel'; 'group'='domain'")
  return(ret)}

### not called  ----------------------------------------------------------------

eapRelFromRes <- function(resultsObj) {
  ret <- resultsObj[intersect(which(resultsObj[,"derived.par"] == "rel"), which(resultsObj[,"par"] == "eap")),c("model", "group", "value")]
  colnames(ret) <- car::recode(colnames(ret), "'value'='rel'; 'group'='domain'")
  return(ret)}
