### called by getConquestQ3() or getTamQ3(),
### which are called by getResults() e.g. getConquestResults() or getTamResults()

#tabdt <- function(x, q3MinType){
#  d   <- na.omit(data.frame (x[[1]], x[[2]]))
#  if(nrow(d) < 2) {return(NULL)}
#  tab <- Rfast::Table(x = d[,1], y = d[,2], names=FALSE)
#  if(ncol(tab)==1) {
#      minVal <- 0
#  }  else  {
#      if ( q3MinType == "singleObs" ) {
#        minVal <- min(tab)
#      }  else  {
#        minVal <- min(c(colSums(tab), rowSums(tab)))
#      }
#  }
#  ret <- data.frame ( Var1 = names(x)[1], Var2 = names(x)[2], minValue = minVal)
#  return(ret)}


### Hilfsfunktion zur Bestimmung der Anzahl der Beobachtungen je Itempaar
#nObsItemPairs <- function ( responseMatrix, q3MinType) {
#        a   <- as.list(data.frame(responseMatrix))
#        spl <- do.call("rbind", combinat::combn(x=a, m=2, fun = tabdt, simplify=FALSE, q3MinType=q3MinType))
#        return(spl)}

### ----------------------------------------------------------------------------
fast_pairwise_tables <- function(data, useNA = c("no", "ifany")) {
  useNA <- match.arg(useNA)
  data <- as.matrix(data)
  n <- ncol(data)
  cn <- colnames(data)
  if (is.null(cn)) cn <- paste0("V", seq_len(n))
  has_values <- apply(data, 2, function(x) any(!is.na(x)))
  mins <- rep(NA_real_, n)
  maxs <- rep(NA_real_, n)
  mins[has_values] <- apply(data[, has_values, drop = FALSE], 2, min, na.rm = TRUE)
  maxs[has_values] <- apply(data[, has_values, drop = FALSE], 2, max, na.rm = TRUE)
  nlev <- maxs - mins + 1L
  codes <- data
  codes[, has_values] <- sweep(data[, has_values, drop = FALSE], 2, mins[has_values] - 1L, "-")
  storage.mode(codes) <- "integer"
  pairs <- combn(n, 2)
  npairs <- ncol(pairs)
  Var1   <- character(npairs)
  Var2   <- character(npairs)
  minVal <- numeric(npairs)
  minSum <- numeric(npairs)
  for (k in seq_len(npairs)) {
    i <- pairs[1, k]; j <- pairs[2, k]
    Var1[k] <- cn[i]
    Var2[k] <- cn[j]
    if (!has_values[i] || !has_values[j]) {
      minVal[k] <- NA_real_
      minSum[k] <- NA_real_
      next
    }
    ni <- nlev[i]; nj <- nlev[j]
    xi <- codes[, i]; xj <- codes[, j]
    ok <- !is.na(xi) & !is.na(xj)
    idx <- rep(NA_integer_, length(xi))
    idx[ok] <- (xi[ok] - 1L) * nj + xj[ok]
    tab <- tabulate(idx, nbins = ni * nj)
    mat <- matrix(tab, nrow = ni, ncol = nj, byrow = TRUE, dimnames = list(mins[i]:maxs[i], mins[j]:maxs[j]))
    if (useNA == "ifany" && (any(is.na(xi)) || any(is.na(xj)))) {
      na_i <- is.na(xi) & !is.na(xj)
      na_j <- !is.na(xi) & is.na(xj)
      na_both <- is.na(xi) & is.na(xj)
      row_na <- tabulate(xj[na_i], nbins = nj)
      col_na <- tabulate(xi[na_j], nbins = ni)
      corner_na <- sum(na_both)
      mat <- rbind(mat, "NA" = row_na)
      mat <- cbind(mat, "NA" = c(col_na, corner_na))
    }
    minVal[k] <- min(mat)
    minSum[k] <- min(c(rowSums(mat), colSums(mat)))
  }
  return(data.frame(Var1 = Var1, Var2 = Var2, minVal = minVal, minSum = minSum, stringsAsFactors = FALSE))
}

reshapeQ3 <- function ( mat, q3MinObs, q3MinType, nObs) {
             for (ii in 1:(nrow(mat)-1)) { mat[ii,ii:ncol(mat)] <- NA}          ### entferne alles oberhalb der Hauptdiagonale
             matL <- reshape2::melt ( mat , na.rm = TRUE)                       ### das entfernt alle doppelten Eintraege
             if(!is.null(nObs)) {                                               ### dass hier soll nur passieren, wenn Eintraege aus der Q3 Matrix ggf. entfernt werden
                chk1 <- lapply(list(matL, nObs), FUN = function (l) {           ### problem: beide data.frames sind nur mergebar, wenn die reihenfolge der var1-var2 in beiden datensaetzen gleich ist
                        chk  <- do.call("rbind", by(data=l, INDICES = l[,c("Var1", "Var2")], FUN = function(y) {
                                y[,"mrge"] <- paste(sort(as.character(unlist(y[,1:2]))), collapse="_")
                                return(y)}))
                        return(chk)})
                matL <- chk1[[2]] |> dplyr::select(-tidyselect::any_of(c("Var1", "Var2"))) |> merge(y = chk1[[1]], by = "mrge", all = FALSE) |> dplyr::select(-tidyselect::any_of("mrge"))
                if(q3MinType == "singleObs") {col <- "minVal"} else {col <- "minSum"}
                weg  <- which(matL[,col] < q3MinObs)
                if(length(weg)>0) {
                   cat(paste0("   Remove ",length(weg), " of ", nrow(matL), " Q3 values due to 'q3MinObs = ", q3MinObs, "' and 'q3MinType = \"",q3MinType, "\"'.\n"))
                   matL <- matL[-weg,]
                }
                if ( nrow(matL) == 0 ) {cat("   No observations left in Q3 matrix.\n") }
             }
             return(matL)}

