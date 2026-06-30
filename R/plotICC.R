plotICC <- function ( resultsObj, defineModelObj, runModelObj = NULL, items = NULL, personPar = c("WLE", "EAP", "PV"), personsPerGroup = 30, pdfFolder = NULL, smooth = 7 ) {
           personPar  <- match.arg(arg = toupper(personPar), choices = c("WLE", "EAP", "PV"))
           if(smooth<5) {smooth <- 5}
           it  <- itemFromRes ( resultsObj )
           if(!"est" %in% colnames(it) ) { it[,"est"] <- NA }
           if(!"estOffset" %in% colnames(it) ) { it[,"estOffset"] <- NA }
           it[,"est"] <- rowSums(it[,c("est", "estOffset")], na.rm = TRUE)      ### untere Zeilen: wenn 1pl und 2pl gemeinsam im resultsobjekt auftauchen, gibt es fuer 1pl keinen
           if(!"estSlope" %in% colnames(it) ) { it[,"estSlope"] <- 1 }          ### slope parameter; die werte sind NA. Zum Plotten muessen sie daher fuer das Raschmodell auf 1 gesetzt werden
           if(length(which(is.na(it[,"estSlope"]))) > 0) { it[which(is.na(it[,"estSlope"])), "estSlope"] <- 1 }
           eapA<- eapFromRes (resultsObj)                                       ### eap fuer alle; muss wideformat haben!!!
           if(personPar == "WLE") {
              eapA <- wleFromRes(resultsObj)
              colnames(eapA) <- car::recode(colnames(eapA), "'wle_est'='EAP'")
           }
           if(personPar == "PV") {
              eapA <- pvFromRes(resultsObj, toWideFormat = TRUE)
              colnames(eapA) <- car::recode(colnames(eapA), "'pv1'='EAP'")
           }
           checkmate::assert_character(items,  null.ok = TRUE, unique = TRUE,  any.missing = FALSE)
           if((is.null(items) || length(items) > 1)  & is.null(pdfFolder)) {stop("If ICCs for more than one item should be displayed, please specify an output folder for pdf.\n")}
           if(!is.null(pdfFolder)) { grDevices::pdf(file = pdfFolder, width = 10, height = 7.5) }
           if(!is.null(items))  {
              miss <- setdiff(items, it[,"item"])
              if(length(miss)>0) {warning(paste0("Following ",length(miss), " items not included in results object: '",paste(miss,, collapse="', '"),"'."))}
              if(length(intersect(items,it[,"item"]))==0) {stop("No commons items in 'items' and results object.")}
              it <- subset(it, item %in% items)
           }
     ### plotten findet fuer jedes dichotome Item separat statt
           inf <- grep("infit", colnames(it), value=TRUE, ignore.case=TRUE)     ### infit-spalte, es darf nur eine geben
           stopifnot(length(inf) == 1)
           pl  <- do.call("rbind", by(data = it, INDICES = it[,c("model", "item")], FUN = function ( i ) {
                  if(nrow(i) == 1) {                                            ### Plotten findet fuer polytome und dichotome Items in getrennter Weise statt, deshalb muss identifiziert
                     xlm <- c(i[["est"]]+2, i[["est"]]-2)                       ### werden, welche Items polytom sind: das sind die mit mehr als einer Zeile im items-Objekt
                     anf <- -6                                                  # anf <- if ( min(xlm) < -4 ) { anf <- floor(min(xlm)) } else { anf  <- -4}
                     ende<- 6                                                   # ende<- if ( max(xlm) >  4 ) { ende<- ceiling(max(xlm)) } else { ende <- 4}
                     x   <- seq ( anf, ende, l = 400)
                     y   <- exp( i[["estSlope"]]*x - i[["est"]] ) / (1+exp( i[["estSlope"]]*x - i[["est"]] ))
                     plot (x, y, type = "l", main = paste("Item '",as.character(i[["item"]]),"'\n\n",sep=""), xlim = c(-6,6), ylim = c(0,1), xlab = "theta", ylab = "P(X=1)", col = "darkred", cex = 8, lwd = 2)
                     graphics::mtext( paste("Model = ",i[["model"]],"  |  Dimension = ",i[["dimension"]], "  |  difficulty = ",round(i[["est"]], digits = 3),"  |  Infit = ",round(i[[inf]], digits = 3),"\n",sep=""))
                     eap <- eapA[intersect ( which (eapA[,"dimension"] == i[["dimension"]]) , which (eapA[,"model"] == i[["model"]])),]
                     if(inherits(defineModelObj, "defineMultiple")) {           ### Problem: je nachdem ob modelle gesplittet wurden oder nicht, muss der Itemdatensatz woanders gesucht werden ... Hotfix
                         woIst<- which ( lapply ( defineModelObj, FUN = function ( g ) {   g[["analysis.name"]] == i[["model"]] }) == TRUE)
                         stopifnot(length(woIst) == 1)
                         dat  <-defineModelObj[[woIst]][["daten"]]
                     } else {
                         dat  <- defineModelObj[["daten"]]
                     }
     ### Hotfix: ID namen identifizieren
                     id  <- unique(resultsObj[intersect(which(resultsObj[,"type"] == "tech"), which(resultsObj[,"par"] == "ID")),"derived.par"])
                     stopifnot(length(id)==1)
                     prbs<- na.omit ( merge ( dat[,c( "ID", as.character(i[["item"]]))], eap[,c( id, "EAP")], by.x = "ID", by.y = id))
                     anz <- round ( nrow(prbs) / personsPerGroup ) + 1             ### mindestens 'personsPerGroup' Personen pro Gruppe
                     if(anz < 3 ) {anz <- 3}
                     if(anz > smooth) {anz <- round(smooth)}
                     eapQ<- quantile ( prbs[,"EAP"], probs = seq(0,1,l = anz))
                     prbs[,"gr"] <- eatTools::num.to.cat ( x = prbs[,"EAP"], cut.points = eapQ[-c(1,length(eapQ))])
                     prbs<- do.call("rbind", by ( data = prbs, INDICES = prbs[,"gr"], FUN = function ( g ) {
                             g[,"mw"] <- mean(g[,"EAP"])
                             g[,"anz"]<- length(g[,"EAP"])
                             g[,"lh"] <- mean(g[, as.character(i[["item"]]) ])
                             return(g)}))
                     matr<- prbs[!duplicated(prbs[,c("mw", "lh")]),c("mw", "lh")]
                     matr<- data.frame(matr[sort(matr[,"mw"],decreasing=FALSE,index.return=TRUE)$ix,])
                     graphics::points ( x = matr[,"mw"], y = matr[,"lh"], cex = 1, pch = 21, bg = "darkblue")
                     graphics::lines ( x = matr[,"mw"], y = matr[,"lh"], col = "blue", lty = 3, lwd = 3)
                  } else {
                     return(i)
                  }} ))
     ### jetzt die partial credit items dazu plotten
           if(!is.null(pl) && nrow(pl)>0) {
              if(inherits(defineModelObj, "defineMultiple")) {                  ### Problem: je nachdem ob modelle gesplittet wurden oder nicht, treten die items in einem oder mehreren objekten auf
                 dfm <- defineModelObj
                 rmo <- runModelObj
              } else {
                 dfm <- list(obj1 = defineModelObj)
                 rmo <- list(obj1 = runModelObj)
              }
              pl2 <- lapply(1:length(dfm), FUN = function(nr) {
                     i <- unique(intersect(colnames(dfm[[nr]][["daten"]]), pl[,"item"]))
                     i2<- eatTools::whereAre(i, rmo[[nr]]$item$item)
                     tx<- capture.output(p2 <- plot(rmo[[nr]], items = i2, type="items", export=FALSE, low=-6, high=6)) })
           }
           if(!is.null(pdfFolder)) { grDevices::dev.off() } }
