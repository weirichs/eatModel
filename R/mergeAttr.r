# Funktion temporaer nach eatModel geschrieben, da sie noch nicht in einem cran-paket enthalten ist

### mergen mit Attributen, das kann 'merge()' nicht:
### http://stackoverflow.com/questions/20306853/maintain-attributes-of-data-frame-columns-after-merge
mergeAttr <- function ( x, y, by = intersect(names(x), names(y)), by.x = by, by.y = by, all = FALSE, all.x = all, all.y = all, sort = TRUE, suffixes = c(".x",".y"), setAttr = TRUE, onlyVarValLabs = TRUE, homoClass = TRUE, unitName = "unit", xName = "x", yName = "y") {
     ### das muessen data.frames sein
             if(length(class(x))>1 || !"data.frame" %in% class(x)) {
                message(paste0("'",xName,"' must be a data.frame. Convert '",xName,"' to data.frame by applying `as.data.frame(x)`.") )
                x <- as.data.frame(x)
             }
             if(length(class(y))>1 || !"data.frame" %in% class(y)) {
                message("'",yName,"' must be a data.frame. Convert '",yName,"' to data.frame by applying `as.data.frame(y)`.")
                y <- as.data.frame(y)
             }
             byvars<- data.frame ( x=by.x, y=by.y, clx = sapply(x[,by.x,drop=FALSE], class), cly = sapply(y[,by.y,drop=FALSE], class), stringsAsFactors = FALSE)
     ### pruefen, ob die level der by-variablen in dem anderen datensatz enthalten sind
             levs  <- apply(X=byvars, MARGIN = 1, FUN = function (v) {
                      nix <- setdiff(unique(y[,v[["y"]]]), unique(x[,v[["x"]]]))
                      if(length(nix)>0) {cat(paste0(length(nix), " ",unitName,"(s) of merging variable '",v[["y"]],"' from data set '",yName,"' not included in data set '",xName,"'.\n"))}
                      niy <- setdiff(unique(x[,v[["x"]]]), unique(y[,v[["y"]]]))
                      if(length(niy)>0) {cat(paste0(length(niy), " ",unitName,"(s) of merging variable '",v[["x"]],"' from data set '",xName,"' not included in data set '",yName,"'.\n"))} })
     ### pruefen, ob die level der by-variablen unique sind
             if ( nrow(byvars)>1) {
                   xby <- unlist(plyr::alply(x, .margins = 1, .fun = function (z) {paste(as.vector(unlist(set.col.type(z[,byvars[,"x"]], col.type = list("character" = byvars[,"x"])))),collapse="_")}))
                   yby <- unlist(plyr::alply(y, .margins = 1, .fun = function (z) {paste(as.vector(unlist(set.col.type(z[,byvars[,"y"]], col.type = list("character" = byvars[,"y"])))),collapse="_")}))
             }  else  {
                   xby <- x[,byvars[1,"x"]]
                   yby <- y[,byvars[1,"y"]]
             }
             if ( length(xby) != length(unique(xby))) { cat("Merging levels are not unique in data set '",xName,"'.\n")}
             if ( length(yby) != length(unique(yby))) { cat("Merging levels are not unique in data set '",yName,"'.\n")}
     ### von allen by-variablen die Klassen homogenisieren, falls gewuenscht
             for ( i in 1:nrow(byvars) ) {
                   if ( length(unique(unlist(byvars[i,c("clx", "cly")]))) > 1 ) {
                        if ( isTRUE(homoClass)) {
                            cat(paste0("   Merging variable pair '", paste(unlist(byvars[i,c("x", "y")]), collapse = "'<==>'"), "' has different classes: '", paste(unlist(byvars[i,c("clx", "cly")]), collapse = "'<==>'"),"'. Classes will be homogenized to 'character'.\n   Use 'homoClass = FALSE' to suppress this behavior.\n"))
                            if ( byvars[i,"clx"] != "character" ) { x[, byvars[i,"x"]] <- as.character(x[, byvars[i,"x"]]) }
                            if ( byvars[i,"cly"] != "character" ) { y[, byvars[i,"y"]] <- as.character(y[, byvars[i,"y"]]) }
                        }  else  {
                            cat(paste0("   Merging variable pair '", paste(unlist(byvars[i,c("x", "y")]), collapse = "'<==>'"), "' has different classes: '", paste(unlist(byvars[i,c("clx", "cly")]), collapse = "'<==>'"),"'.\n   Use 'homoClass = TRUE' to homogenize classes.\n"))
                        }
                   }
             }
     ### jetzt mergen und DANACH die Attribute rekonstruieren
             datM  <- merge ( x=x, y=y, by.x=by.x, by.y=by.y, all=all, all.x=all.x, all.y=all.y, sort=sort, suffixes =suffixes)
             if ( isTRUE(setAttr) ) {
                   dats<- list(x=x, y=y)
                   for ( d in names(dats)) {
                         for ( v in colnames(dats[[d]])) {
                               vsuf <- paste0(v, suffixes[2])
                               if ( vsuf %in% colnames(datM) ) {
                                    if ( onlyVarValLabs == FALSE ) {
                                         if(!is.null(attributes(dats[[d]][,v]))) {attributes(datM[,vsuf]) <- attributes(dats[[d]][,v])}
                                    }  else  {
                                         if(!is.null(attr(dats[[d]][,v], "varLabel"))) {attr(datM[,vsuf], "varLabel") <- attr(dats[[d]][,v], "varLabel")}
                                         if(!is.null(attr(dats[[d]][,v], "valLabel"))) {attr(datM[,vsuf], "valLabel") <- attr(dats[[d]][,v], "valLabel")}
                                    }
                               }  else  {
                                    if ( v %in% colnames(datM) ) {
                                         if ( onlyVarValLabs == FALSE ) {
                                              if(!is.null(attributes(dats[[d]][,v]))) {attributes(datM[,v]) <- attributes(dats[[d]][,v])}
                                         }  else  {
                                              if(!is.null(attr(dats[[d]][,v], "varLabel"))) {attr(datM[,v], "varLabel") <- attr(dats[[d]][,v], "varLabel")}
                                              if(!is.null(attr(dats[[d]][,v], "valLabel"))) {attr(datM[,v], "valLabel") <- attr(dats[[d]][,v], "valLabel")}
                                         }
                                    }
                               }
                         }
                   }
             }
             return(datM)}
