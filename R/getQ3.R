### called by getConquestQ3() or getTamQ3(),
### which are called by getResults() e.g. getConquestResults() or getTamResults()

tabdt <- function(x, q3MinType){
  d   <- na.omit(data.frame (x[[1]], x[[2]]))
  if(nrow(d) < 2) {return(NULL)}
  tab <- Rfast::Table(x = d[,1], y = d[,2], names=FALSE)
  if(ncol(tab)==1) {
      minVal <- 0
  }  else  {
      if ( q3MinType == "singleObs" ) {
        minVal <- min(tab)
      }  else  {
        minVal <- min(c(colSums(tab), rowSums(tab)))
      }
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

