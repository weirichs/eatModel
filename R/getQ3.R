### called by getConquestQ3() or getTamQ3(),
### which are called by getResults() e.g. getConquestResults() or getTamResults()

tabdt <- function(x, q3MinType){
  d   <- na.omit(data.frame (x[[1]], x[[2]]))
  if(nrow(d) < 2) {return(NULL)}
  tab <- Rfast::Table(x = d[,1], y = d[,2], names=FALSE)
  if ( q3MinType == "singleObs" ) {
    minVal <- min(tab)
  }  else  {
    minVal <- min(c(colSums(tab), rowSums(tab)))
  }
  ret <- data.frame ( Var1 = names(x)[1], Var2 = names(x)[2], minValue = minVal)
  return(ret)}

### Hilfsfunktion zur Bestimmung der Anzahl der Beobachtungen je Itempaar
nObsItemPairs <- function ( responseMatrix, q3MinType) {
        a   <- as.list(data.frame(responseMatrix))
        spl <- do.call("rbind", combinat::combn(x=a, m=2, fun = tabdt, simplify=FALSE, q3MinType=q3MinType))
        return(spl)}

### ----------------------------------------------------------------------------

reshapeQ3 <- function ( mat, q3MinObs, nObs ) {
  for (ii in 1:(nrow(mat)-1)) { mat[ii,ii:ncol(mat)] <- NA}
  matL <- reshape2::melt ( mat , na.rm = TRUE)
  if ( !is.null(nObs)) {
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
