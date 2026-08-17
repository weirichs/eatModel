### called by getConquestQ3() or getTamQ3(),
### which are called by getResults() e.g. getConquestResults() or getTamResults()

.q3TabLevels <- function(x) {
  if(is.factor(x)) {return(levels(x))}
  return(sort(unique(x[!is.na(x)])))}

.q3TabCodes <- function(x, levels) {
  if(is.factor(x)) {x <- as.character(x)}
  return(match(x, levels))}

.q3MinTabValue <- function(x, y, q3MinType) {
  xLev <- .q3TabLevels(x)
  yLev <- .q3TabLevels(y)
  nX   <- length(xLev)
  nY   <- length(yLev)
  if(nX < 2 || nY < 2) {return(0)}
  xCod <- .q3TabCodes(x, xLev)
  yCod <- .q3TabCodes(y, yLev)
  ok   <- !is.na(xCod) & !is.na(yCod)
  if(!any(ok)) {return(0)}
  tab <- matrix(tabulate(xCod[ok] + (yCod[ok] - 1L) * nX, nbins = nX * nY),
                nrow = nX, ncol = nY)
  if ( q3MinType == "singleObs" ) {
    minVal <- min(tab)
  }  else  {
    minVal <- min(c(colSums(tab), rowSums(tab)))
  }
  return(minVal)}

tabdt <- function(x, q3MinType){
  minVal <- .q3MinTabValue(x = x[[1]], y = x[[2]], q3MinType = q3MinType)
  ret <- data.frame ( Var1 = names(x)[1], Var2 = names(x)[2], minValue = minVal, stringsAsFactors = FALSE)
  return(ret)}


### Hilfsfunktion zur Bestimmung der Anzahl der Beobachtungen je Itempaar
nObsItemPairs <- function ( responseMatrix, q3MinType) {
        q3MinType <- match.arg(q3MinType, choices = c("singleObs", "marginalSum"))
        dat <- data.frame(responseMatrix, check.names = FALSE)
        nItems <- ncol(dat)
        if(nItems < 2) {
          return(data.frame(Var1 = character(0), Var2 = character(0), minValue = numeric(0), stringsAsFactors = FALSE))
        }
        itemNames <- colnames(dat)
        itemLevels <- lapply(dat, FUN = .q3TabLevels)
        nLevels <- lengths(itemLevels)
        itemCodes <- Map(f = .q3TabCodes, x = dat, levels = itemLevels)
        itemPairs <- utils::combn(seq_len(nItems), 2)
        minValue <- numeric(ncol(itemPairs))
        for(pp in seq_len(ncol(itemPairs))) {
          ii <- itemPairs[1, pp]
          jj <- itemPairs[2, pp]
          if(nLevels[ii] < 2 || nLevels[jj] < 2) {
            minValue[pp] <- 0
          } else {
            iiCodes <- itemCodes[[ii]]
            jjCodes <- itemCodes[[jj]]
            complete <- !is.na(iiCodes) & !is.na(jjCodes)
            if(!any(complete)) {
              minValue[pp] <- 0
            } else {
              tab <- matrix(tabulate(iiCodes[complete] + (jjCodes[complete] - 1L) * nLevels[ii],
                                     nbins = nLevels[ii] * nLevels[jj]),
                            nrow = nLevels[ii], ncol = nLevels[jj])
              if(q3MinType == "singleObs") {
                minValue[pp] <- min(tab)
              } else {
                minValue[pp] <- min(c(rowSums(tab), colSums(tab)))
              }
            }
          }
        }
        spl <- data.frame(Var1 = itemNames[itemPairs[1,]], Var2 = itemNames[itemPairs[2,]],
                          minValue = minValue, stringsAsFactors = FALSE)
        return(spl)}

### ----------------------------------------------------------------------------

reshapeQ3 <- function ( mat, q3MinObs, nObs ) {
  for (ii in 1:(nrow(mat)-1)) { mat[ii,ii:ncol(mat)] <- NA}          ### entferne alles oberhalb der Hauptdiagonale
  matL <- reshape2::melt ( mat , na.rm = TRUE)                       ### das entfernt alle doppelten Eintraege
  if ( !is.null(nObs)) {                                             ### dass hier soll nur passieren, wenn Eintraege aus der Q3 Matrix ggf. entfernt werden
    chk1 <- lapply(list(matL, nObs), FUN = function (l) {
      vals <- l[,3]
      chk  <- do.call("rbind", apply(l[,-ncol(l)], MARGIN = 1, FUN = function ( y ) { ret <- sort ( y); ret <- data.frame ( Var1 = ret[1], Var2 = ret[2], stringsAsFactors = FALSE); return(ret)}))
      l2   <- data.frame ( chk, X = l[,3], stringsAsFactors = FALSE)
      colnames(l2)[3] <- colnames(l)[3]
      return(l2)})
    matL <- eatTools::na_omit_selection(merge ( chk1[[1]], chk1[[2]], by = c("Var1", "Var2"), all = TRUE), "value")
    weg  <- which(matL[,"minValue"] < q3MinObs)
    if (length(weg)>0) { matL <- matL[-weg,]}
  }
  if ( nrow(matL) == 0 ) {
    cat("No observations left in Q3 matrix.\n")
  }
  return(matL)}

