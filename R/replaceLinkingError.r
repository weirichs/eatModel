# equatingList: The object returned by \code{\link{equat1pl}}, for measurement occasion 1 vs. 3!
# multiEquatError_output: The object returned by \code{\link{multiEquatError}}
replaceLinkingError <- function(equatingList, multiEquatError_output, verbose = TRUE, digits = 4) {
    ### checks
       if ( !is.list(multiEquatError_output)&&is.null(dim(multiEquatError_output)) ) {stop("'multiEquatError_output' must be a list of data.frames.")}
       checkmate::assert_logical(verbose, len = 1)
       checkmate::assert_numeric(digits, len = 1)
    ### dimensionen in der Equatingliste. Diese Dimensionen muessen auch in dem 'multiEquatError_output'-Objekt vorhanden sein
       dims <- lapply(equatingList[["items"]], names)
       drin <- unlist(dims) %in% names(multiEquatError_output)
       if(!all(drin)) { stop(paste0("Dimension(s) '",paste(dims[which(drin==FALSE)],collapse="', '"),"' are missing in the 'multiEquatError_output'."))}
       for (d in 1:length(dims)) {
            oldLE <- equatingList[["items"]][[names(dims)[d]]][[dims[[d]]]][["eq"]][["descriptives"]][["linkerror"]]
            if ( verbose) {message(paste0("Dimension '",dims[[d]],"': Replace old linking error ",round(oldLE, digits = digits), " with ", round(multiEquatError_output[[dims[[d]]]][["le3221"]], digits = digits)))}
            equatingList[["items"]][[names(dims)[d]]][[dims[[d]]]][["eq"]][["descriptives"]][["linkerror"]] <- multiEquatError_output[[dims[[d]]]][["le3221"]]
       }
       return(equatingList)}


