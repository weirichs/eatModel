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
