multiEquatError <- function (x1, x2, x3, difBound = 1, dependentDIF =FALSE, testletStr = NULL, verbose=TRUE ){
       liste<- list(x1, x2, x3)
       chk1 <- checkInput(liste)
       if ( "derived.par" %in% colnames(x1) ) {                                 
            obj  <- prepareAndCheckEatModelObject(liste, difBound=difBound, verbose=verbose)
            link <- lapply(names(obj), FUN = function (dim ) { tripleEquatError(e1=obj[[dim]][[1]][,c("item", "est")], e2=obj[[dim]][[2]][,c("item", "est")], e3= obj[[dim]][[3]][,c("item", "est")], dependentDIF=dependentDIF, testletStr=testletStr, difBound=difBound, calledForEatModel = TRUE)})
            names(link) <- names(obj)
       } else {                                                                 
            link <- tripleEquatError(e1=x1, e2=x2, e3=x3, dependentDIF=dependentDIF, testletStr=testletStr, difBound=difBound, calledForEatModel = FALSE)
       }
       return(link) }


tripleEquatError <- function(e1, e2, e3, dependentDIF, testletStr, difBound, calledForEatModel) {
      chk1 <- checkInputConsistency(e1=e1,e2=e2,e3=e3,testletStr=testletStr)    
      if (isFALSE(calledForEatModel)) {                                         
           chk2 <- checkIPD(e1=e1,e2=e2,e3=e3,testletStr=chk1, difBound=difBound)
           if ( length(chk2[["e1_weg"]])>0) {e1 <- e1[-which(e1[,1] %in% chk2[["e1_weg"]]),]}
           if ( length(chk2[["e3_weg"]])>0) {e3 <- e3[-which(e3[,1] %in% chk2[["e3_weg"]]),]}
      }
      el21 <- equ.rasch(e2,e1)
      el31 <- equ.rasch(e3,e1)
      el32 <- equ.rasch(e3,e2)
      l21 <- el21$descriptives
      l31 <- el31$descriptives
      l32 <- el32$descriptives
      trend3221 <- el32$B.est$Mean.Mean + el21$B.est$Mean.Mean
      trend31 <- el31$B.est$Mean.Mean
    if(is.null(chk1)) {                                                         
     if(dependentDIF) {                                                         
       is3221 <- intersect(el21$anchor$item, el32$anchor$item)                  
       dif21 <- el21$anchor$TransfItempar.Gr1[match(is3221, el21$anchor$item)] - el21$anchor$Itempar.Gr2[match(is3221, el21$anchor$item)]
       dif32 <- el32$anchor$TransfItempar.Gr1[match(is3221, el32$anchor$item)] - el32$anchor$Itempar.Gr2[match(is3221, el32$anchor$item)]
        cov1223a <- cov(dif21, dif32)
        if(cov1223a < 0) {cov1223a <- 0}
        cov1223 <- sqrt(cov1223a)/sqrt(length(dif21))
        le3221 <- sqrt(l21$linkerror^2 + l32$linkerror^2 - 2*cov1223^2)
      } else {
        le3221 <- sqrt(l21$linkerror^2 + l32$linkerror^2)
      }
      le31 <- l31$linkerror
  } else {
      e12a <- merge(e1, e2, by = "item", all=TRUE)
      e12 <- merge(chk1, e12a, by="item", all.x=FALSE, all.y=TRUE)
      e12 <- e12[,c(2,3,4,1)]
      e13a <- merge(e1, e3, by = "item", all=TRUE)
      e13 <- merge(chk1, e13a, by="item", all.x=FALSE, all.y=TRUE)
      e13 <- e13[,c(2,3,4,1)]
      e23a <- merge(e2, e3, by = "item", all=TRUE)
      e23 <- merge(chk1, e23a, by="item", all.x=FALSE, all.y=TRUE)
      e23 <- e23[,c(2,3,4,1)]
      eli12 <- equa.rasch.jk(e12)
      eli13 <- equa.rasch.jk(e13)
      eli23 <- equa.rasch.jk(e23)
      l12 <- eli12$descriptives$linkerror.jackknife
      le31 <- eli13$descriptives$linkerror.jackknife
      l23 <- eli23$descriptives$linkerror.jackknife
      if(dependentDIF) {
        eli12$pars.data[,2] <- eli12$pars.data[,2]+eli12$descriptives$shift
        eli23$pars.data[,2] <- eli23$pars.data[,2]+eli23$descriptives$shift
        e1223 <- stats::na.omit(merge(eli12$pars.data, eli23$pars.data, by=c("unit","item"),all=TRUE))
        e1223$dif12 <- e1223[,3] - e1223[,4]
        e1223$dif23 <- e1223[,5] - e1223[,6]
        res1 <- NULL
        for(un in e1223$unit) {
          p1223 <- e1223[e1223[,1] != un,]
          res1 <- c(cov(p1223$dif12,p1223$dif23),res1)
        }
        cov1223a <- mean(res1, na.rm=TRUE)
        if(cov1223a < 0) {cov1223a <- 0}
        cov1223 <- sqrt(cov1223a)/sqrt(length(p1223$dif12))
        if((l12^2 + l23^2 - 2*cov1223^2) < 0) {
          le3221 <- sqrt(l12^2 + l23^2)
        } else {
          le3221 <- sqrt(l12^2 + l23^2 - 2*cov1223^2)
        }
      } else {
        le3221 <- sqrt(l12^2 + l23^2)
      }
  }
  res <- data.frame(trend31 = trend31, trend3221 = trend3221, le31 = le31, le3221 = le3221)
  return(res)
}

checkInput <- function(inputlist) {
     cls <- lapply(inputlist, class)
     if(!all.equal(cls[[1]], cls[[2]], cls[[3]])) {stop("'x1', 'x2', and 'x3' must have the same class.")}
     if(!is.data.frame(inputlist[[1]])) {stop("'x1', 'x2', and 'x3' must be of class 'data.frame'.")} }
     

prepareAndCheckEatModelObject <- function ( liste, difBound, verbose ) {
       uni  <- sapply(liste, FUN = function ( x ) { "correlation" %in% x[,"par"] })
       chk1 <- which(uni == TRUE)
       if ( length(chk1)>0) {stop("Following model(s) do not seem to be unidimensional: ",paste(chk1, collapse = ", "),".")}
       doms <- lapply(liste, FUN = function (o) {na.omit(setdiff(unique(o[,"group"]), "persons"))})
       comps<- data.frame(combinat::combn(1:length(doms),2))
       chk2 <- sapply(comps, FUN = function ( col ) { if(!all.equal(doms[[col[1]]], doms[[col[2]]])) {stop(paste0("Domains between ",col[1]," and ",col[2]," do not match."))}})
       its  <- lapply(liste, itemFromRes)
       dims <- unique(its[[1]][,"dimension"])
       its  <- lapply(dims, FUN = function (d) {  list(t1 = its[[1]][which(its[[1]][,"dimension"] == d),], t2 = its[[2]][which(its[[2]][,"dimension"] == d),], t3 = its[[3]][which(its[[3]][,"dimension"] == d),]) })
       names(its) <- dims
       chk3 <- lapply(dims, FUN = function ( d ) {
               bridges <- data.frame(combinat::combn(1:length(its[[d]]),2), stringsAsFactors = FALSE)
               bridges <- bridges[,which(sapply(bridges, diff)==1)]
               chk4    <- lapply(bridges, FUN = function (b ) {
                          com<- length(intersect( its[[d]][[b[1]]][,"item"], its[[d]][[b[2]]][,"item"]))
                          if ( com == 0 ) {stop(paste0("No common items for linking mzp ",b[1]," to mzp ", b[2]," for domain '",d,"'."))}
                          if ( com < 10 ) {cat(paste0("Warning: Only ",com," common items for linking mzp ",b[1]," to mzp ", b[2]," for domain '",d,"'.\n"))}
                          })})
       ref <- itemFromRes(liste[[1]])
       if (verbose) {
           eq  <- equat1pl(liste[[2]], prmNorm = ref[,c("item", "est")], difBound = difBound, iterativ = TRUE)
       }  else  {
           txt <- capture.output(eq  <- equat1pl(liste[[2]], prmNorm = ref[,c("item", "est")], difBound = difBound, iterativ = TRUE))
       }
       stru<- lapply(eq[["items"]], names)
       for ( d in 1:length(stru)) {
           weg <- setdiff(eq[["items"]][[names(stru)[d]]][[stru[[d]]]][["info"]][["itemExcluded"]], "")
           if ( length(weg)>0) {
                its[[stru[[d]]]][[1]] <- its[[stru[[d]]]][[1]][-match(weg, its[[stru[[d]]]][[1]][,"item"]),]
           }
       }
       ref2<- itemFromRes(liste[[2]])
       if (verbose) {
           eq  <- equat1pl(liste[[3]], prmNorm = ref2[,c("item", "est")], difBound = difBound, iterativ = TRUE)
       }  else  {
           txt  <- capture.output(eq  <- equat1pl(liste[[3]], prmNorm = ref2[,c("item", "est")], difBound = difBound, iterativ = TRUE))
       }
       stru<- lapply(eq[["items"]], names)
       for ( d in 1:length(stru)) {
           weg <- setdiff(eq[["items"]][[names(stru)[d]]][[stru[[d]]]][["info"]][["itemExcluded"]], "")
           if ( length(weg)>0) {
                its[[d]][[3]] <- its[[d]][[3]][-match(weg, its[[d]][[3]][,"item"]),]
           }
       }
       return(its)}
       
checkInputConsistency <- function(e1,e2,e3,testletStr) {
       il  <- list(e1, e2, e3)
       if(!all(sapply(il, ncol) == 2)) {stop("All thre data.frames must have two columns.")}
       if(!all(lapply(il, FUN = function ( x ) { all(colnames(x) == c("item", "est")) })==TRUE)) {
          il <- lapply(il, FUN = function ( x ) { colnames(x) <- c("item", "est"); return(x)})
       }
       it  <- lapply(il, FUN = function (x) {
                     stopifnot ( length(x[,"item"]) == length(unique(x[,"item"])))
                     return(x[,"item"])})
       comp<- as.data.frame(combinat::combn(length(it),2))
       chk1<- lapply(comp, FUN = function ( com ){ if ( length(intersect(it[[com[1]]], it[[com[2]]])) < 2) { stop(paste0(length(intersect(it[[com[1]]], it[[com[2]]])) ," common items between measurement occasions ",com[1], " and ", com[2], "."))} })
       if (!is.null(testletStr)) {
            if(!"data.frame" %in% class(testletStr) || "tbl" %in% class(testletStr) ) { stop("'testletStr' must be a data.frame.")}
            if ( ncol(testletStr) != 2) {stop("'testletStr' must have two columns.")}
            if(!"item" %in% colnames(testletStr)) { stop("One column in 'testletStr' must be named 'item'.")}
            if (setdiff(colnames(testletStr), "item") != "unit") {
                cat(paste0("The 'unit' column in 'testletStr' seemed to be called '",setdiff(colnames(testletStr), "item"),"'. Rename this column into 'unit'.\n"))
                colnames(testletStr)[setdiff(1:2, match("item", colnames(testletStr)))] <- "unit"
            }
            allI<- unique(unlist(lapply(il, FUN = function (x){x[,"item"]})))
            if(!all(allI %in% testletStr[,"item"])) {
                message(paste0("Following items without testlet definition: '",paste(setdiff(allI,testletStr[,"item"]) , collapse="', '"), "'.\n    Linking errors will be compute without considering testlets."))
                testletStr <- NULL
            }
            tls <- testletStr[which(testletStr[,"item"] %in% allI),]
            isna<- length(which(is.na(tls[,setdiff(colnames(tls), "item")])))
            if(isna>0) {
                message(paste0("Following ",isna," items without testlet definition: '",paste(tls[which(is.na(tls[,setdiff(colnames(tls), "item")])),"item"] , collapse="', '"), "'.\n    Linking errors will be compute without considering testlets."))
                testletStr <- NULL
            }
            return(testletStr)
       }  else  {
            return(NULL)
       } }

equa.rasch.jk <- function (pars.data, se.linkerror = FALSE, alpha1 = 0, alpha2 = 0) {
  pars.data <- as.data.frame(stats::na.omit(pars.data))
  itemunits <- unique(pars.data[, 1])
  N.units <- length(itemunits)
  N.items <- nrow(pars.data)
  mod1 <- equ.rasch(pars.data[, c(4, 2)], pars.data[, c(4, 3)])
  res1 <- data.frame(unit = itemunits, shift = 0, SD = 0, linkerror = 0)
  for (nn in 1:N.units) {
    pars.data1 <- pars.data[pars.data[, 1] != itemunits[nn], ]
    mod.nn <- equ.rasch(x = pars.data1[, c(4, 2)], y = pars.data1[, c(4, 3)])
    res1[nn, "shift"] <- mod.nn$B.est$Mean.Mean
    res1[nn, "SD"] <- mod.nn$descriptives$SD
    if (se.linkerror) {
      itemunits.nn <- itemunits[-nn]
      l1 <- NULL
      for (ii in itemunits.nn) {
        pars.data1.ii <- pars.data1[paste(pars.data1[, 1]) != ii, ]
        mod.ii <- equ.rasch(x = pars.data1.ii[,  c(4, 2)], y = pars.data1.ii[, c(4, 3)], alpha1 = alpha1, alpha2 = alpha2)
        l1 <- c(l1, mod.ii$B.est$Mean.Mean)
      }
      res1[nn, "linkerror"] <- sqrt((N.units - 2)/(N.units -  1) * sum((l1 - res1[nn, "shift"])^2))
    }
  }
  linkerror <- sqrt((N.units - 1)/N.units * sum((res1[, 2] - mod1$B.est$Mean.Mean)^2))
  se.sd <- sqrt((N.units - 1)/N.units * sum((res1[, 3] - mod1$descriptives$SD)^2))
  if (se.linkerror) {
    se.linkerror <- sqrt((N.units - 1)/N.units * sum((res1[, 4] - linkerror)^2))
  }
  else {
    se.linkerror <- NA
  }
  descriptives <- data.frame(N.items = N.items, N.units = N.units,
                             shift = mod1$B.est$Mean.Mean, SD = mod1$descriptives$SD,
                             linkerror.jackknife = linkerror, SE.SD.jackknife = se.sd,
                             se.linkerror.jackknife = se.linkerror)
  res <- list(pars.data = pars.data, itemunits = itemunits, descriptives = descriptives)
  return(res)
}

equ.rasch <- function (x, y, theta = seq(-4, 4, len = 100), alpha1 = 0, alpha2 = 0) {
  x[, 1] <- gsub(" ", "", paste(x[, 1]))
  y[, 1] <- gsub(" ", "", paste(y[, 1]))
  b.xy <- stats::na.omit(merge(x, y, by = "item"))
  colnames(b.xy) <- c("item", "Itempar.Gr1", "Itempar.Gr2")
  B.mm <- mean(b.xy[, 3]) - mean(b.xy[, 2])
  g1   <- prob_raschtype_genlogis(theta = theta, b = b.xy[, 2], alpha1 = 0, alpha2 = 0)
  opt_interval <- 10 * c(-1, 1)
  ha <- function(B) {
    fct1 <- prob_raschtype_genlogis(theta = theta, b = b.xy[, 2], alpha1 = alpha1, alpha2 = alpha2)
    fct2 <- prob_raschtype_genlogis(theta = theta, b = b.xy[, 3] - B, alpha1 = alpha1, alpha2 = alpha2)
    sum((fct1 - fct2)^2)
  }
  B.ha <- stats::optimize(f = ha, interval = opt_interval)$minimum
  sl <- function(B) {
    fct1 <- prob_raschtype_genlogis(theta = theta, b = b.xy[, 2], alpha1 = alpha1, alpha2 = alpha2)
    fct2 <- prob_raschtype_genlogis(theta = theta, b = b.xy[, 3] - B, alpha1 = alpha1, alpha2 = alpha2)
    sum((rowSums(fct1 - fct2))^2)
  }
  B.sl <- stats::optimize(f = sl, interval = opt_interval)$minimum
  B.est <- data.frame(B.mm, B.ha, B.sl)
  colnames(B.est) <- c("Mean.Mean", "Haebara","Stocking.Lord")
  b.xy$TransfItempar.Gr1 <- b.xy[, 2] + B.est[1, "Mean.Mean"]
  x[, 2] <- x[, 2] + B.est[1, "Mean.Mean"]
  transf.par <- merge(x = x, y = y, by.x = 1, by.y = 1, all = TRUE)
  colnames(transf.par) <- c("item", "TransfItempar.Gr1", "Itempar.Gr2")
  transf.par <- transf.par[order(paste(transf.par$item)), ]
  des <- data.frame(N.Items = nrow(b.xy), SD = stats::sd(b.xy$TransfItempar.Gr1 -  b.xy$Itempar.Gr2), stringsAsFactors = FALSE)
  des$Var <- des$SD^2
  des$linkerror <- sqrt(des[,"SD"]^2/des[,"N.Items"])
  res <- list(B.est = B.est, descriptives = des, anchor = b.xy[,  c(1, 2, 4, 3)], transf.par = transf.par)
  return(res)}

prob_raschtype_genlogis <- function( theta, b, alpha1, alpha2, fixed.a=1+0*b,  Qmatrix=NULL, dimensions=NULL ) {
  LT <- length(theta)
  if (is.matrix(theta)){
    LT <- nrow(theta)
  }

  if ( is.null(Qmatrix) ){
    XX <- TAM::tam_outer(x=theta, y=b, op="-")
    XX <- as.vector(XX)
    aM <- sirt::sirt_matrix2(x=fixed.a, nrow=LT)
    XX <- aM * XX
  }
  if ( ! is.null(Qmatrix) ){
    XX0 <- tcrossprod( as.matrix(theta), Qmatrix )
    XX <- XX0 - outer( rep(1,nrow(theta)), b )
    XX <- as.vector(XX)
    aM <- sirt::sirt_matrix2(x=fixed.a, nrow=LT)
    XX <- aM * XX
  }
  pm <- sirt::pgenlogis(x=XX, alpha1=alpha1, alpha2=alpha2 )
  pm <- matrix( pm, ncol=length(b))
  return(pm)}
  
checkIPD <- function(e1,e2,e3,testletStr, difBound){
       eq1  <- equat1pl(results=e2, itemF = colnames(e2)[1], valueF = colnames(e2)[2], prmNorm = e1, difBound=difBound, iterativ=TRUE)
       eq2  <- equat1pl(results=e3, itemF = colnames(e3)[1], valueF = colnames(e3)[2], prmNorm = e2, difBound=difBound, iterativ=TRUE)
       if ( "itemExcluded" %in% colnames(eq1[["items"]][["global"]][["global"]][["info"]])) {
            e1_weg <- setdiff(eq1[["items"]][["global"]][["global"]][["info"]][,"itemExcluded"], "")
       }  else  {
            e1_weg <- NULL
       }
       if ( "itemExcluded" %in% colnames(eq2[["items"]][["global"]][["global"]][["info"]])) {
            e3_weg <- setdiff(eq2[["items"]][["global"]][["global"]][["info"]][,"itemExcluded"], "")
       }  else  {
            e3_weg <- NULL
       }
       return(list(e1_weg = e1_weg, e3_weg = e3_weg))}








