checkLinking <- function(design, blocks=NULL, bookletColumn=NULL, verbose=FALSE) {
      if(!all(sapply(design, inherits, what="character"))) {design <- data.frame(lapply(design, as.character), stringsAsFactors=FALSE)}
      if(!is.null(bookletColumn)) {
          if(length(bookletColumn) != 1) {stop("Argument 'bookletColumn' must be of length 1.")}
          book  <- eatTools::existsBackgroundVariables(dat = design, variable=bookletColumn)
      } else {
          book <- "bl"                                                          ### dieses komplizierte Zeug, weil der Name der Bookletspalte nicht bereits
          while (book %in% colnames(design)) {book <- paste0(book, sample(0:9,1))}# im Designobjekt vergeben sein darf
          design[,book] <- paste0("T", 1:nrow(design))
      }                                                                         ### untere zeile: langformat mit na.rm = TRUE erspart einem das 'removeNAd()'
      colnames(design)[-match(book, colnames(design))] <- paste0("Pos", 1:(ncol(design)-1))
      desL  <- eatTools::facToChar(reshape2::melt(design, id.vars = book, na.rm=TRUE, variable.name="blockPos", value.name="blockName"))
      if(verbose) {message("Special characters and spaces will be removed from names in 'blocks' and prefix 'B' will be added.")}
      desL[,"blockName"] <- paste0("B", eatTools::gsubAll(desL[,"blockName"], old=c(" ", "[[:punct:]]"), new=c("", "")))
      ### an welcher Position kommt welcher Block mit welcher Haeufigkeit vor?
      blPos <- reshape2::dcast(desL, blockName~blockPos, value.var="blockName", fun.aggregate=length)
      names(blPos)[2:ncol(blPos)] <- paste0("Pos", 1:(ncol(blPos)-1))
      if(!is.null(blocks)) {
          stopifnot(is.vector(blocks))
          if(!inherits(blocks, "character")) {stop("'block' must be a character vector.")}
          blocks <- paste0("B", eatTools::gsubAll(blocks, old=c(" ", "[[:punct:]]"), new=c("", "")))
          stopifnot(length(intersect(desL[,"blockName"], blocks)) > 0)
          notInDe<- setdiff(blocks,unique(desL[,"blockName"]))
          if(length(notInDe)>0) {warning(paste0(length(notInDe), " elements of 'blocks' vector missing in design."))}
          desL   <- desL[which(desL[,"blockName"] %in% blocks),]
          blPos <- blPos[blPos$blockName %in% blocks,]
          blPos <- blPos[order(blPos$blockName),]
      }
    ### welche Kombinationen von Bloecken kommen wie oft vor?
      kombs <- data.frame(do.call("rbind", combinat::combn(unique(desL[,"blockName"]),2, simplify=FALSE)), stringsAsFactors=FALSE)
      kombs[,"pairFreq"] <- apply(kombs, MARGIN = 1, FUN = function(z){ sum(as.vector(by(data = desL, INDICES = desL[,book], FUN = function (b) { all(z %in% b[,"blockName"])})))})
      names(kombs)[1:2] <- c("block.X", "block.Y")
      kombs <- kombs[order(kombs$pairFreq, decreasing = TRUE),]
      desW  <- reshape2::dcast(desL, as.formula(paste0(book, " ~ blockPos")), value.var="blockName")
      dat   <- do.call(plyr::rbind.fill, apply(desW, MARGIN = 1, FUN = simDat, booklet = book))
      link  <- checkLink(dataFrame = dat[,-1, drop = FALSE], remove.non.responser = TRUE, verbose = TRUE)
      res <- list(completelyLinked=link, occuringBlockCombinations=kombs, blockPositions =blPos)
      if(verbose) {
        print(kombs)
        print(res[["blockPositions"]])
      }
      return(res)
      }

### Hilfsfunktion fuer 'checkLinking'
simDat <- function ( z, booklet ) {                                             ### erzeugt Datensatz aus einer Zeile des Designs
          if ( !length(na.omit(z[-match(booklet, names(z))])) == length(unique(na.omit(z[-match(booklet, names(z))] ))) ) { stop("Blocks are not unique in each line.\n")}
          items<- as.vector(sapply(na.omit(z[-match(booklet, names(z))]), FUN= function ( i ) { paste(i, 1:3, sep="_")}))
          pers <- paste(z[[booklet]], 11:22, sep="_")                           ### Funktion muss also ueber "apply" aufgerufen werden!
          mat  <- data.frame ( id = pers, matrix ( sample ( 0:1, size = length(pers) * length(items), replace = TRUE), ncol = length(items), nrow = length(pers)))
          colnames(mat)[-1] <- items
          return(mat)}

checkLink <- function(dataFrame, remove.non.responser = FALSE, sysmis = NA, verbose = TRUE)   {
             if(!is.na(sysmis))  {
               na <- which(is.na(dataFrame))
               if(length(na)>0)  {
                  warning(paste0("'",sysmis,"' was specified to denote 'sysmis' in the data. ",length(na)," 'NA'-values were found in the dataset anyway. \n         Hence, ",sysmis," and 'NA' will be handled as 'sysmis'."))
               }
               dataFrame <- as.data.frame(lapply(dataFrame, FUN=function(ii) {car::recode(ii, paste(sysmis,"= NA", collapse = "; ") ) } ) )
             }
             if ( remove.non.responser == TRUE ) {
                na <- which( rowSums(is.na(dataFrame)) == ncol ( dataFrame ) )
                if(length(na)>0) {
                   dataFrame <- dataFrame[-na,]
                   if(verbose == TRUE ) {cat(paste("Remove ",length(na)," cases with missing on all items.\n", sep = ""))}
                }
             }
             non.missing.cases <- lapply(dataFrame, FUN=function(ii) {which(!is.na(ii))})
             all.cases <- non.missing.cases[[1]]
             i <- 2
             total.abbruch     <- FALSE
             while( (i < length(non.missing.cases) + 1 ) & !total.abbruch )  {
                  if(length( intersect(all.cases,non.missing.cases[[i]])) > 0 )  {
                     all.cases <- unique(c(all.cases, non.missing.cases[[i]] ) )
                  }  else   {
                     overlap        <- FALSE
                     remain.columns <- length(non.missing.cases) + 1 - i
                     ii             <- 1
                     while (overlap == FALSE & ii < remain.columns )  {
                           non.missing.cases <- non.missing.cases[c(setdiff(1:length(non.missing.cases),i),i)]
                          if(length( intersect(all.cases,non.missing.cases[[i]])) > 0 ) {overlap <- TRUE}
                           ii <- ii + 1
                     }
                     if (overlap == FALSE) {total.abbruch <- TRUE}
                     if (overlap == TRUE)  {all.cases <- unique(c(all.cases, non.missing.cases[[i]] ) ) }
                  }
                  i <- i + 1
             }
             if (length(all.cases) != nrow(dataFrame))   {
                if (verbose == TRUE) {cat("WARNING! Dataset is not completely linked.\n") }
                if ( remove.non.responser == TRUE ) {
                     missed <- setdiff ( 1:nrow(dataFrame), all.cases)
                     misFra <- reshape2::melt ( data.frame ( id = 1:length(missed), dataFrame[missed,]), id.vars = "id", na.rm=TRUE)
                     cat ( paste ( "W A R N I N G !   Dataset is NOT completely linked (even if cases with missings on all items are removed).\n                  ",length(missed)," cases unconnected. Following items are unconnected: \n",sep=""))
                     cat("                  "); cat ( paste ( unique(as.character(misFra[,"variable"])), collapse = ", ")); cat("\n")
                }
                return(FALSE)
             }
             if (length(all.cases) == nrow(dataFrame))   {
                if (verbose == TRUE) {cat("Dataset is completely linked.\n") }
                return(TRUE)
             }  }
