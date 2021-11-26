# equatingList: The object returned by \code{\link{equat1pl}}, for measurement occasion 1 vs. 3!
# multiEquatError_output: The object returned by \code{\link{multiEquatError}}
replaceLinkingError <- function(equatingList, multiEquatError_output, verbose = TRUE, digits = 4) {
    ### checks
       if ( !is.list(multiEquatError_output)&&is.null(dim(multiEquatError_output)) ) {stop("'multiEquatError_output' must be a list of data.frames.")}
    ### dimensionen in der Equatingliste. Diese Dimensionen muessen auch in dem 'multiEquatError_output'-Objekt vorhanden sein
       dims <- lapply(equatingList[["items"]], names)
       drin <- unlist(dims) %in% names(multiEquatError_output)
       if(!all(drin)) { stop(paste0("Dimension(s) '",paste(dims[which(drin==FALSE)],collapse="', '"),"' are missing in the 'multiEquatError_output'."))}
       for (d in 1:length(dims)) {
            oldLE <- equatingList[["items"]][[names(dims)[d]]][[dims[[d]]]][["eq"]][["descriptives"]][["linkerror"]]
            if ( verbose) {message(paste0("Dimension '",dims[[d]],"': Replace old linking error ",round(oldLE, digits = digits), " with ", round(multiEquatError_output[[dims[[d]]]][["le1223"]], digits = digits)))}
            equatingList[["items"]][[names(dims)[d]]][[dims[[d]]]][["eq"]][["descriptives"]][["linkerror"]] <- multiEquatError_output[[dims[[d]]]][["le1223"]]
       }
       return(equatingList)}


### Ergaenzen in hilfeseite von multiEquatError!
# im Anschluss an das allererste Beispiel dort!

# direct linking 1 to 3 (1 is reference)
# it1 <- itemFromRes(results[[1]])
# eq1.vs.3 <- equat1pl(results[[3]], prmNorm = it1[,c("item", "est")], difBound = 0.64, iterativ = TRUE)
# eq1.vs.3 <- replaceLinkingError (equatingList=eq1.vs.3, multiEquatError_output=lErrors)

# zweites Beispiel
# direct linking 1 to 3 (1 is reference)
# it1 <- itemFromRes(results2[[1]])
# eq1.vs.3 <- equat1pl(results2[[3]], prmNorm = it1[,c("item", "est")], difBound = 0.64, iterativ = TRUE)
# eq1.vs.3 <- replaceLinkingError (equatingList=eq1.vs.3, multiEquatError_output=lErrors2)
