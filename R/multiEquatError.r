# wird auf namespace exportiert, soll sowohl data.frames als auch results-objekte verarbeiten koennen
# inputobjekte heissen unspezifisch x1, x2, x3, weil es ja sowohl data.frames als auch results-objekte sein koennen
multiEquatError <- function (eq.1_2, eq.2_3, eq.1_3, dependentDIF =FALSE, verbose=TRUE ){
       lapply(list(dependentDIF, verbose=TRUE), checkmate::assert_logical, any.missing = FALSE, len = 1)
       alls <- list(eq.1_2, eq.2_3, eq.1_3)
       dims <- lapply(alls, FUN = function (x) {names(x[["items"]])})
       stopifnot(all(dims[[1]] == dims[[2]]) && all(dims[[1]] == dims[[3]]))
       dims <- dims[[1]]
       LEs  <- lapply(dims, FUN = function (d) {                                ### fuer jede Dimension wird die Berechnung separat durchgefuehrt
               l21        <- eq.1_2[["items"]][[d]][[1]][["eq"]][["descriptives"]]
               l31        <- eq.1_3[["items"]][[d]][[1]][["eq"]][["descriptives"]]
               l32        <- eq.2_3[["items"]][[d]][[1]][["eq"]][["descriptives"]]
               method12   <- eq.1_2[["items"]][[d]][[1]][["method"]]             ### feststellen, welche Linkingmethode (mean-mean, Haebara, ... verwendet wurde)
               method13   <- eq.1_3[["items"]][[d]][[1]][["method"]]
               method23   <- eq.2_3[["items"]][[d]][[1]][["method"]]
               if((method12 != method13) || (method23 != method13)) {warning("Inconsistent linking methods between 1-2, 2-3, and 1-3.")}
               trend3221  <- eq.2_3[["items"]][[d]][[1]][["eq"]][["B.est"]][[method23]] + eq.1_2[["items"]][[d]][[1]][["eq"]][["B.est"]][[method12]]
               trend31    <- eq.1_3[["items"]][[d]][[1]][["eq"]][["B.est"]][[method13]]
               if(isFALSE(dependentDIF)) {                                      ### hier: independent DIF
                  le3221  <- sqrt(l21$linkerror^2 + l32$linkerror^2)
               } else {                                                         ### hier: dependent DIF, aber erstmal ohne testlets
                  is3221  <- intersect(eq.1_2[["items"]][[d]][[1]][["eq"]][["anchor"]][["item"]], eq.2_3[["items"]][[d]][[1]][["eq"]][["anchor"]][["item"]])
                  dif21   <- subset(eq.1_2[["items"]][[d]][[1]][["eq"]][["anchor"]], item %in% is3221) |> dplyr::mutate(diff = TransfItempar.Gr1 - Itempar.Gr2)
                  dif32   <- subset(eq.2_3[["items"]][[d]][[1]][["eq"]][["anchor"]], item %in% is3221) |> dplyr::mutate(diff = TransfItempar.Gr1 - Itempar.Gr2)
                  if(!is.null(eq.1_2[["items"]][[d]][[1]][["eq"]][["ntl"]])) {  ### es gibt Testlets: jackknife bestimmung der Kovarianz
                     testl   <- merge(dif21[, c("item", "diff")], dif32[, c("item", attr(eq.2_3, "allN")[["testlet"]], "diff")], by="item")
                     check1  <- merge(dif21[, c("item", attr(eq.1_2, "allN")[["testlet"]])], dif32[, c("item", attr(eq.2_3, "allN")[["testlet"]])], by="item")
                     if( !all(check1[,2] == check1[,3])) {warning(paste0("Dimension '",d,"': Inconsistent testlet definition between measurement 1 vs. 2 and 2 vs. 3."))}
                     cov1223a<- unlist(lapply(unique(testl[,attr(eq.2_3, "allN")[["testlet"]]]), FUN = function (tl) {
                                tl <- testl[testl[,attr(eq.2_3, "allN")[["testlet"]]] != tl,]
                                return(cov(tl[,"diff.x"], tl[,"diff.y"]))})) |> mean()
                  } else {                                                      ### es gibt keine Testlets
                     cov1223a<- cov(dif21[,"diff"], dif32[,"diff"])
                  }
                  if(cov1223a < 0) {
                     if(verbose) {message(paste0("Dimension '",d,"': Negative covariance of ",round(cov1223a, digits=3)," between DIF for measurement 1 vs. 2, and DIF for measurement 2 vs. 3. Covariance is set to 0 to avoid underestimation of the true linking error."))}
                     cov1223a <- 0
                  }
                  cov1223 <- sqrt(cov1223a)/sqrt(length(dif21))
                  le3221  <- sqrt(l21$linkerror^2 + l32$linkerror^2 - 2*cov1223^2)
               }
               le31 <- l31$linkerror
               if(verbose) {message(paste0("\nDimension '",d,"': Direct linking error of the three combinations of measurements: \n",eatTools::print_and_capture(eatTools::roundDF(plyr::rbind.fill(data.frame ( mzp1 = 1, mzp2 = 2, l21), data.frame ( mzp1 = 1, mzp2 = 3, l31, chained = le3221, mzp_1.vs.3= trend31, mzp_1.vs.2.vs.3 = trend3221), data.frame ( mzp1 = 2, mzp2 = 3, l32)), digits = 3), spaces = 5)))}
               res <- data.frame(trend31 = trend31, trend3221 = trend3221, le31 = le31, le3221 = le3221)
               return(res)})
       names(LEs)  <- dims
       return(LEs)}


replaceLinkingError <- function(equatingList, multiEquatError_output, verbose = TRUE, digits = 4) {
    ### checks
       if ( !is.list(multiEquatError_output) && is.null(dim(multiEquatError_output)) ) {stop("'multiEquatError_output' must be a list of data.frames.")}
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


