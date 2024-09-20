### called by getConquestQ3() or getTamQ3(),
### which are called by getResults() e.g. getConquestResults() or getTamResults()

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
