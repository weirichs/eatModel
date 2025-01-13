equat1pl<- function(results, prmNorm, item = NULL, domain = NULL, testlet = NULL,
                    value = NULL, excludeLinkingDif = TRUE, difBound = 1, iterativ = FALSE,
                    method = c("Mean.Mean", "Haebara", "Stocking.Lord", "robust", "Haberman"),
                    itemF = NULL, domainF = NULL, testletF = NULL, valueF = NULL,
                    estimation=c("OLS", "BSQ", "HUB", "MED", "LTS", "L1", "L0"),
                    b_trim=Inf, lts_prop=.5) {

### checking/assering the arguments --------------------------------------------
  estimation <- match.arg(estimation)
  method     <- match.arg(method)

  # data frame (from getResults() or df with parameters of focus group)
  results<- eatTools::makeDataFrame(results)
  lapply(c(itemF, domainF, testletF, valueF), checkmate::assert_vector, null.ok = TRUE)

  # columns for data frame 'prmNorm' with normed anchor item parameters
  lapply(c(item, domain, testlet, value), checkmate::assert_vector, null.ok = TRUE)

  checkmate::assert_numeric(difBound, len = 1)

  # logicals
  lapply(c(excludeLinkingDif, iterativ), checkmate::assert_logical, len = 1)

  # if method = "Haberman"
  lapply(c(b_trim, lts_prop), checkmate::assert_numeric, len = 1)

### function -------------------------------------------------------------------

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
          prmM<- allN[["prmNorm"]] [ which(allN[["prmNorm"]][,allN[["domain"]]] %in% unique(it[,"dimension"])) ,]
        }  else  {
          prmM<- allN[["prmNorm"]]
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
