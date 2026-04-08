simEquiTable <- function ( anchor, item = NULL, cat = NULL, value = NULL, mRef, sdRef, addConst = 500, multConst = 100, cutScores) {
                lapply(list(mRef, sdRef, addConst, multConst), checkmate::assert_numeric, len=1,any.missing=FALSE)
                anchor<- eatTools::makeDataFrame(anchor)
                if(ncol(anchor) == 2) {                                         ### dichotome Aequivalenztabelle
                   if(!inherits(anchor[,2], c("integer", "numeric"))) {stop("Item parameter column must be numeric.")}
                   if(length(unique(anchor[,1])) != nrow(anchor)) {stop("Item ID column has duplicated entries.")}
                   if(length(which(is.na(anchor[,2]))) > 0) {
                       warning(paste0(length(which(is.na(anchor[,2]))), " missing values in item parameter values. Missings will be removed."))
                       anchor <- na.omit(anchor)
                   }
    ### temporaeren Datensatz mit allen moeglichen Summenscores erzeugen
                   dtmp  <- data.frame(rbind(1*(lower.tri(matrix(1, nrow = nrow(anchor), ncol = nrow(anchor)))),1))
                   dtmp  <- data.frame(dtmp, score = rowSums(dtmp) , irtoys::wle(dtmp, cbind(1, anchor[,2], 0)), stringsAsFactors = FALSE)
                } else {                                                        ### polytome Aequivalenztabelle
                   if(is.null(item) || is.null(value) || is.null(cat)) { stop("If 'anchor' has more than two columns (partial credit model), 'item', 'cat', and 'value' columns must be specified explicitly.") }
                   allF  <- list(item=item, cat=cat, value = value)             ### erzeuge dichotome pseudo-items, und zwar so viele wie es maximale summenwerte gibt
                   allF  <- lapply(allF, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = anchor, variable=ii)})
                   items <- by(anchor, INDICES = anchor[,allF[["item"]]], FUN = function(i) {0:nrow(i)})
                   items1<- sum(unlist(lapply(items, max)))
                   dtmp  <- data.frame(rbind(1*(lower.tri(matrix(1, nrow = items1, ncol = items1))),1))
                   pcitem<- items[which(sapply(items, length)>2)]
                   diItem<- setdiff(names(items), names(pcitem))
                   for(i in names(pcitem)) {
                       ncat <- length(pcitem[[i]])-1
                       dtmp <- data.frame(rowSums(dtmp[,(ncol(dtmp)-ncat+1):ncol(dtmp)]), dtmp[,1:(ncol(dtmp)-ncat)], stringsAsFactors = FALSE)
                       colnames(dtmp)[1] <- i}
                   colnames(dtmp)[(1+length(pcitem)):ncol(dtmp)] <- diItem
                   txt   <- capture.output(def3  <- dtmp |> dplyr::mutate(id = paste0("P", 1:nrow(dtmp))) |> defineModel(items = colnames(dtmp), id = "id",  irtmodel = "PCM",software="tam", anchor=it, itemCol = "item", valueCol = "est", catCol = "category", verbose=FALSE))
                   run3  <- runModel(def3)
                   res3  <- getResults(run3, omitPV=TRUE, verbose=FALSE)
                   dtmp  <- wleFromRes(res3) |> dplyr::rename(est=wle_est, score = NitemsSolved, n = NitemsTotal)
                }
                dtmp[,"bista"] <- (dtmp[,"est"] - mRef) / sdRef * multConst + addConst
                dtmp[,"ks"]    <- eatTools::num.to.cat ( x = dtmp[,"bista"], cut.points = cutScores[["values"]], cat.values = cutScores[["labels"]])
                dtmp  <- data.frame(dtmp[sort(dtmp[,"score"],decreasing=FALSE,index.return=TRUE)$ix,])
    ### jetzt noch die shortversion der Aequivalenztabelle erzeugen ... dazu rauskriegen, auf wieviele Stellen gerundet werden darf ... begonnen wird bei 0, dann wird geschaut, ob die gerundeten cuts identisch mit Kompetenzstufenschwelle sind. Falls ja, wird eine dezimalstelle dazugetan
                dig   <- 0
                while ( length(which(round(dtmp[,"bista"], digits = dig) %in% cutScores[["values"]])) > 0) {dig <- dig + 1}
                shrt  <- do.call("rbind", by ( data = dtmp, INDICES = dtmp[,"ks"], FUN = function ( sks ) { data.frame ( score = paste(c(min(sks[,"score"]), max(sks[,"score"])), collapse=" bis "), estimate = paste(round(c(min(sks[,"est"]), max(sks[,"est"])),digits=2), collapse=" bis "), bista = paste(round(c(min(sks[,"bista"]), max(sks[,"bista"])),digits=dig), collapse=" bis "), ks=unique(sks[,"ks"]), stringsAsFactors=FALSE)}))
                return(list ( complete = dtmp[,c("score", "est", "bista", "ks")], short = shrt))}
