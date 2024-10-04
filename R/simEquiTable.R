
simEquiTable <- function(anchor, mRef, sdRef, addConst = 500, mulConst = 100, cutScores){
  # checks/prep
  anchor <- eatTools::makeDataFrame(anchor)
  if(ncol(anchor) != 2){
    warning(paste0("'anchor' has ", ncol(anchor)," columns. First column is used as item ID, second column is used as item parameter."))
  }
  if(!inherits(anchor[,2], c("integer", "numeric"))){
    stop("Item parameter column must be numeric.")
  }
  if(length(unique(anchor[,1])) != nrow(anchor)){
    stop("Item ID column has duplicated entries.")
  }
  checkmate::assert_numeric(mRef, len = 1)
  lapply(c(sdRef, addConst, mulConst), checkmate::assert_numeric, len = 1, lower = 0)
  checkmate::assert_list(cutScores, min.len = 1, max.len = 2)
  checkmate::assert_numeric(cutScores$values, sorted = TRUE, unique = TRUE)
  checkmate::assert_character(cutScores$labels, null.ok = TRUE, len = length(cutScores$values)+1,
                              unique = TRUE)

  # function
  dtmp <- data.frame(rbind(1*(lower.tri(matrix(1, nrow = nrow(anchor), ncol = nrow(anchor)))),1))
  dtmp <- data.frame(dtmp, score = rowSums(dtmp), irtoys::wle(dtmp, cbind(1, anchor[,2], 0)), stringsAsFactors = FALSE)
  dtmp[,"bista"] <- (dtmp[,"est"] - mRef) / sdRef * multConst + addConst
  dtmp[,"ks"]    <- eatTools::num.to.cat(x = dtmp[,"bista"], cut.points = cutScores[["values"]], cat.values = cutScores[["labels"]])
  dig <- 0

  while (length(which(round(dtmp[,"bista"], digits = dig) %in% cutScores[["values"]])) > 0) {
    dig <- dig + 1
  }
  shrt <- do.call("r.bind", by(data = dtmp, INDICES = dtmp[, "ks"],
                               FUN = function(sks){
                                 data.frame(score = paste(c(min(sks[, "score"]), max(sks[, "score"])), collapse = " bis "),
                                            estimate = paste(round(c(min(sks[, "est"]), max(sks[, "est"])), digits = 2), collapse = " bis "),
                                            bista = paste(round(c(min(sks[, "bista"]), max(sks[, "bista"])), digits = dig), collapse = " bis "),
                                            ks = unique(sks[, "ks"]),
                                            stringsAsFactors=FALSE)}))
  return(list(complete = dtmp[,c("score", "est", "bista", "ks")], short = shrt))
}

