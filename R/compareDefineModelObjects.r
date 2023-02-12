compareDefineModelObjects <- function (ref, tar, round = TRUE, digits = 3) {
       mult <- unlist(lapply(list(ref, tar), FUN = function (x) {length(grep("ultiple", class(x)))}))
       if(!sum(mult) %in% c(0,2) ) { message(paste0("Objects are incomparable. '", car::recode(which(mult == 1), "1='ref'; 2='tar'"), "' is multiple, '", car::recode(which(mult == 0), "1='ref'; 2='tar'"), "' is single.")); return()}
       if (sum(mult) == 2) {                                                    
           if (!length(ref) == length(tar)) {message(paste0("Objects are incomparable. 'ref' includes ", length(ref), " models, 'tar' includes ", length(tar), " models.")); return()}
           nams <- lapply(list(ref, tar), FUN = function (x) {sapply(x, FUN = function (y) {y[["analysis.name"]]})})
           if ( !all(sort(nams[[1]]) == sort (nams[[2]])) ) { message(paste0("Cannot match incomparable 'analysis.name': \n    'ref' = '",paste(sort(nams[[1]]), collapse="', '"), "'. \n    'tar' = '",paste(sort(nams[[2]]), collapse="', '"), "'.")); return()}
           ord  <- nams[[1]]
       }  else  {
           ord  <- 1
       }                                                                        
       co   <- lapply( ord, FUN = function(i) {                                 
               if(inherits(i, "character")) {                                   
                   nams <- sapply(list(ref, tar), FUN = function (x) {which(sapply(x, FUN = function (y) {y[["analysis.name"]]==i}))})
                   dats <- list(ref = ref[[nams[1]]], tar = tar[[nams[2]]])     
               }  else  {                                                       
                   dats <- list(ref=ref, tar=tar)
               }})
       vgl  <- lapply(co, FUN = function (v) {                                  
               stopifnot(length(v) == 2)                                        
               foo<- lapply(c("software", "constraint", "irtmodel", "n.plausible", "estVar", "pvMethod", "fitTamMmlForBayesian"), FUN = function (strng) {if ( !v[[1]][[strng]] == v[[2]][[strng]]) {message(paste0("'ref': ",strng," = '",v[[1]][[strng]], "', 'tar': ",strng," = '",v[[2]][[strng]],"'."))} })
               d  <- lapply(v, FUN =function(v1) { reshape2::melt(v1[["daten"]],id.vars ="ID", na.rm=TRUE, variable.name="item")})
               e  <- eatTools::mergeAttr(d[[1]], d[[2]], by=c("ID", "item"), all=TRUE, setAttr=FALSE, xName="reference", yName="target", verbose=TRUE)
               for ( j in unique(c(names(v[[1]][["all.Names"]]), names(v[[2]][["all.Names"]])))) {
                     eq <- all.equal(sort(v[[1]][["all.Names"]][[j]]), sort(v[[2]][["all.Names"]][[j]]))
                     if(!isTRUE(eq)) {message(paste0("Mismatch for entry '",j,"' of 'all.Names' object: '",paste( eq, collapse= "'. '"),"'"))}
               }
               qm <- compareQmatrices(qm1=v[[1]][["qMatrix"]], qm2=v[[2]][["qMatrix"]])
               if ( !is.null(v[[1]][["anchor"]]) ) {
                     stopifnot(!is.null(v[[2]][["anchor"]]))
                     ank<- eatTools::mergeAttr(data.frame ( item = v[[1]][["all.Names"]][["variablen"]], ank = v[[1]][["anchor"]][,2], stringsAsFactors = FALSE),  data.frame ( item = v[[2]][["all.Names"]][["variablen"]], ank = v[[2]][["anchor"]][,2], stringsAsFactors = FALSE), by="item", all=TRUE, setAttr=FALSE, xName="reference anchor item set", yName="target anchor item set", verbose=TRUE)
                     if(isFALSE(round)) {
                         apv<- all.equal(ank[,2], ank[,3])
                     }  else  {
                         apv<- all.equal(round(ank[,2], digits=digits), round(ank[,3], digits=digits))
                     }
                     if(!isTRUE(apv)) {message(paste0("Anchor parameter mismatch: '",paste( apv, collapse= "'. '"),"'"))}
               }
               if(!is.null(v[[1]][["est.slopegroups"]]) || !is.null(v[[2]][["est.slopegroups"]]) ) {
                     eq <- all.equal(v[[1]][["est.slopegroups"]], v[[2]][["est.slopegroups"]])
                     if(!isTRUE(eq)) {message(paste0("Slope group mismatch: '",paste( eq, collapse= "'. '"), "'"))}
               }
               if(!is.null(v[[1]][["guessMat"]])) {message("Guessing matrices comparison not yet implemented.")}
               nam<- unique(c(names(v[[1]][["control"]]), names(v[[2]][["control"]])))
               for ( j in nam) {
                     eq <- all.equal(v[[1]][["control"]][[j]], v[[2]][["control"]][[j]])
                     if(!isTRUE(eq)) {message(paste0("Control mismatch for object '",j,"': '",paste( eq, collapse= "'. '"), "'"))}
               }
               if(!is.null(v[[1]][["fixSlopeMat"]]) || !is.null(v[[2]][["fixSlopeMat"]]) ) {message("Fix slope matrices comparison not yet implemented.")}
               })}

compareQmatrices <- function (qm1, qm2){
       if ( ncol(qm1) != ncol(qm2) || !all(sort(colnames(qm1)) == sort(colnames(qm2))) ) {
            message(paste0("Colnames of Q matrices differ: \n   reference: '", paste(colnames(qm1), collapse="', '"), "\n      target: '", paste(colnames(qm2), collapse="', '")))
            if( ncol(qm1) == ncol(qm2) ) {                                      
                colnames(qm2) <- eatTools::recodeLookup(colnames(qm2), data.frame ( alt = colnames(qm2), neu = colnames(qm1), stringsAsFactors = FALSE))
            }
       }
       qL   <- lapply(list(qm1, qm2), FUN = function (qmlist) {reshape2::melt(qmlist, id.vars = colnames(qm2)[1], variable.name = "dim", na.rm=TRUE)})
       qLM  <- eatTools::mergeAttr(qL[[1]], qL[[2]], by=colnames(qL[[1]])[1:2], all=TRUE, setAttr=FALSE, xName="reference Q matrix", yName="target Q matrix", verbose=TRUE)
       if(!isTRUE(all.equal(qLM[,"value.x"], qLM[,"value.y"]))) {message(paste0("Q matrix mismatch: '",paste( all.equal(qLM[,"value.x"], qLM[,"value.y"]), collapse= "'. '"))) }
       }

