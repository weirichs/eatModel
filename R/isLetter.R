### called by defineModel/.substituteSigns() and getResults/getConquestAdditionalTerms()

isLetter <- function(string){
  # function
  splt <- strsplit(string, "")
  isL <- lapply(splt, FUN = function(x) {
    ind <- which(x %in% c(letters, LETTERS))
    x[setdiff(1:length(x), ind)] <- " "
    x <- eatTools::crop(paste(x, sep = "", collapse = ""))
    x <- unlist(strsplit(x, " +"))
    return(x)
  })
  return(isL)
}

