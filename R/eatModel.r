### diverse functions that can be used on their own or are called by
### multiple other functions which could't be usefully grouped.

### called by .substituteSigns() and getConquestAdditionalTerms () -------------

isLetter <- function(string){
  # checks
  checkmate::assert_character(string)

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

### not called ----------------------------------------------------------------------------

simEquiTable <- function ( anchor, mRef, sdRef, addConst = 500, multConst = 100, cutScores) {
                anchor<- eatTools::makeDataFrame(anchor)
                if ( ncol(anchor) != 2) {
                     warning(paste0("'anchor' has ",ncol(anchor)," columns. First column is used as item ID, second column is used as item parameter."))
                }
                if(!inherits(anchor[,2], c("integer", "numeric"))) {stop("Item parameter column must be numeric.")}
                if(length(unique(anchor[,1])) != nrow(anchor)) {stop("Item ID column has duplicated entries.")}
                dtmp  <- data.frame(rbind(1*(lower.tri(matrix(1, nrow = nrow(anchor), ncol = nrow(anchor)))),1))
                dtmp  <- data.frame(dtmp, score = rowSums(dtmp) , irtoys::wle(dtmp, cbind(1, anchor[,2], 0)), stringsAsFactors = FALSE)
                dtmp[,"bista"] <- (dtmp[,"est"] - mRef) / sdRef * multConst + addConst
                dtmp[,"ks"]    <- eatTools::num.to.cat ( x = dtmp[,"bista"], cut.points = cutScores[["values"]], cat.values = cutScores[["labels"]])
                dig   <- 0
                while ( length(which(round(dtmp[,"bista"], digits = dig) %in% cutScores[["values"]])) > 0) {dig <- dig + 1}
                shrt  <- do.call("rbind", by ( data = dtmp, INDICES = dtmp[,"ks"], FUN = function ( sks ) { data.frame ( score = paste(c(min(sks[,"score"]), max(sks[,"score"])), collapse=" bis "), estimate = paste(round(c(min(sks[,"est"]), max(sks[,"est"])),digits=2), collapse=" bis "), bista = paste(round(c(min(sks[,"bista"]), max(sks[,"bista"])),digits=dig), collapse=" bis "), ks=unique(sks[,"ks"]), stringsAsFactors=FALSE)}))
                return(list ( complete = dtmp[,c("score", "est", "bista", "ks")], short = shrt))}

### not called -----------------------------------------------------------------

getConquestVersion <- function ( path.conquest , path.temp , asDate = TRUE ) {
    wd <- path.temp
		f <- file.path ( wd , "delete.cqc" )
		write ( "quit;" , f )
		f <- normalizePath ( f )
		path.conquest <- normalizePath ( path.conquest )
		cmd <- paste ( "\"", path.conquest, "\" \"", f , "\"" , sep ="")
		r <- NULL
		suppressWarnings(try ( r <- system ( command = cmd , intern = TRUE ) , silent = TRUE ))
		file.remove ( f )
		if ( !is.null ( r ) ) {
				r <- r[1]
				r <- sub ( "ConQuest build: " , "" , r )
				r <- gsub ( "\\s+" , "-" , r )
				if ( asDate ) r <- date::as.date(r)
		}
		return (r)}

### not called -----------------------------------------------------------------

prepRep <- function ( calibT2, bistaTransfT1, bistaTransfT2, makeIdsUnique = TRUE) {
           if ( !inherits(calibT2, "transfBista" )) { stop("'calibT2' object must be of class 'transfBista'.\n")}
           if ( !inherits(bistaTransfT1, "transfBista" )) { stop("'bistaTransfT2' object must be of class 'transfBista'.\n")}
           if ( !inherits(bistaTransfT2, "transfBista") ) { stop("'bistaTransfT2' object must be of class 'transfBista'.\n")}
           if (!nrow(calibT2[["itempars"]]) < nrow(bistaTransfT1[["itempars"]])) { stop("Mismatch between 'calibT2' and 'bistaTransfT1'. \n")}
           if (!nrow(calibT2[["itempars"]]) < nrow(bistaTransfT2[["itempars"]])) { stop("Mismatch between 'calibT2' and 'bistaTransfT2'. \n")}
           idT1<- unique(bistaTransfT1[["all.Names"]][which(bistaTransfT1[["all.Names"]][,"par"] == "ID"),"derived.par"])
           idT2<- unique(bistaTransfT2[["all.Names"]][which(bistaTransfT2[["all.Names"]][,"par"] == "ID"),"derived.par"])
           stopifnot(length(idT1)==1, length(idT2)==1)
           if ( idT1 != idT2 ) {
                warning(paste0("ID variables do not match between t1 and t2. ID for t1: '",idT1,"'. ID for t2: '",idT2,"'. \n    IDs will be unified with '",idT1,"'."))
                recStat <- paste ( "'", idT2 , "' = '", idT1, "'", sep="")
                colnames ( bistaTransfT2[["personpars"]] ) <- car::recode ( colnames ( bistaTransfT2[["personpars"]] ), recStat)
           }
           lc  <- colnames( calibT2[["personpars"]] ) [grep("^linking", colnames(calibT2[["personpars"]]) )]
           if(length(lc)==0) { stop("No columns with linking error information found in 'calibT2'.\n")}
           lcn <- paste("trend", eatTools::removePattern(string = lc, pattern = "linking"), sep="")
           colnames( calibT2[["personpars"]] ) [grep("^linking", colnames(calibT2[["personpars"]]) )] <- lcn
           merg<- c("group", "imp", "traitLevel", "dimension")
           frms<- list ( calibT2=calibT2, bistaTransfT1=bistaTransfT1, bistaTransfT2=bistaTransfT2 )
           toM <- unique(unlist(lapply ( names(frms), FUN = function ( l.Name ) {
                  l    <- frms[[l.Name]]
                  drin <- merg %in% colnames(l[["personpars"]])
                  fehlt<- merg[which(drin==FALSE)]
                  if (!all(drin == TRUE)) { warning(paste0("Column(s) '",paste(fehlt, collapse = "', '"), "' are unexpectedly missing in '",l.Name,"'."))}
                  keep <- merg[which(drin==TRUE)]
                  return(keep)})))
           if(length(toM)==0) { stop("Merging impossible.\n")}
           red <- calibT2[["personpars"]][,c(toM,  lcn)]
           red <- red[!duplicated(red),]
           dat1<- data.frame ( trend = "T1" , merge ( bistaTransfT1[["personpars"]], red, by = toM, all = TRUE))
           stopifnot ( nrow(dat1) == nrow(bistaTransfT1[["personpars"]]))
           dat2<- data.frame ( trend = "T2" , merge ( bistaTransfT2[["personpars"]], red, by = toM, all = TRUE))
           stopifnot ( nrow(dat2) == nrow(bistaTransfT2[["personpars"]]))
           if ( makeIdsUnique == TRUE ) {
                dat1[, paste(idT1, "unique", sep="_")] <- paste(dat1[, "trend"], dat1[, idT1], sep="_")
                dat2[, paste(idT1, "unique", sep="_")] <- paste(dat2[, "trend"], dat2[, idT1], sep="_")
           }
           return(rbind ( dat1, dat2))}
