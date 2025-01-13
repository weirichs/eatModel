### called by equat1pl()

### also called by transformItemParListIntoResults () --------------------------

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
  if(isTRUE(allF[["testlet"]] == "dimension")) {                       ### hier muss isTRUE stehen, weil es sonst fehlschlaegt, wenn xx == yy logical(0) ergibt
    message(paste0("'dimension' is not allowed for testlet column name in 'prmNorm'. Rename column to 'dimensionName'."))
    allF[["testlet"]] <- "dimensionName"
    colnames(prmNorm) <- car::recode(colnames(prmNorm), "'dimension'='dimensionName'")
  }
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
  allF[["prmNorm"]] <- prmNorm
  return(allF)}

### ----------------------------------------------------------------------------

### hilfsfunktion fuer equat1pl: transformiert itemparameterliste in results-objekt, wenn nicht das Rueckgabeobjekt von getResults()
### die Fokusparameter enthaelt, sondern ein einfacher data.frame
transformItemParListIntoResults <- function(results, itemF, domainF, testletF, valueF){
  allF <- checkItemParLists(prmNorm =results, item = itemF, domain = domainF, testlet = testletF, value = valueF)
  if (!is.null(allF[["domain"]])) {
    dims <- names(table(allF[["prmNorm"]][,allF[["domain"]]]))
  }  else  {
    allF[["domain"]] <- "domaene"
    allF[["prmNorm"]][, allF[["domain"]]] <- dims <- "global"
  }
  allF[["prmNorm"]][,"model"] <- allF[["prmNorm"]][, allF[["domain"]]]
  weg  <- intersect ( colnames (allF[["prmNorm"]] ) , setdiff  ( c("item", "dimension", "est"), unlist(allF) ))
  if ( length ( weg ) > 0 )  {                                         ### damit keine Spalten durch 'recode' doppelt benannt werden, muessen spalten, die sich durch die Recodierung aendern
    allF[["prmNorm"]] <- allF[["prmNorm"]][, -match(weg, colnames(allF[["prmNorm"]]))]
  }                                                                    ### und zugleich schon im datensatz 'results' vergeben sind, raus
  allF <- allF[which(!sapply(allF, is.null))]
  toRec<- lapply(setdiff(names(allF),"prmNorm"), FUN = function ( ff ) { paste ( "'",allF[[ff]],"'='",car::recode(ff, "'item'='item'; 'domain'='dimension'; 'value'='est'"),"'",sep="")})
  toRec<- paste(toRec, collapse = "; ")
  colnames(allF[["prmNorm"]]) <- car::recode (colnames(allF[["prmNorm"]]), toRec)
  return(list(results=allF[["prmNorm"]], dims=dims))}

### ----------------------------------------------------------------------------

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

### ----------------------------------------------------------------------------

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

### ----------------------------------------------------------------------------

handleLinkingDif <- function(prmDim,prbl, eq, difBound, dif, method, excludeLinkingDif, iterativ,prmM, allN) {
  cat(paste ( "\nDimension '", prmDim[1,"dimension"], "': ", length( prbl), " of ", nrow( eq[["anchor"]]), " items with linking |DIF| > ",difBound," identified.\n",sep=""))
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

### ----------------------------------------------------------------------------

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

### ----------------------------------------------------------------------------

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

### also called by handleLinkingDif() ------------------------------------------

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

