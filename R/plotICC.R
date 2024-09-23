plotICC <- function ( resultsObj, defineModelObj, item = NULL, personPar = c("WLE", "EAP", "PV"), personsPerGroup = 30, pdfFolder = NULL, smooth = 7 ) {
  personPar  <- match.arg(arg = toupper(personPar), choices = c("WLE", "EAP", "PV"))
  if (smooth<5) {smooth <- 5}
  it  <- itemFromRes ( resultsObj )
  if ( !"est" %in% colnames(it) ) { it[,"est"] <- NA }
  if ( !"estOffset" %in% colnames(it) ) { it[,"estOffset"] <- NA }
  it[,"est"] <- rowSums(it[,c("est", "estOffset")], na.rm = TRUE)
  if ( !"estSlope" %in% colnames(it) ) { it[,"estSlope"] <- 1 }
  if ( length(which(is.na(it[,"estSlope"]))) > 0) { it[which(is.na(it[,"estSlope"])), "estSlope"] <- 1 }
  eapA<- eapFromRes (resultsObj)
  if ( personPar == "WLE") {
    eapA <- wleFromRes(resultsObj)
    colnames(eapA) <- car::recode(colnames(eapA), "'wle_est'='EAP'")
  }
  if ( personPar == "PV") {
    eapA <- pvFromRes(resultsObj, toWideFormat = TRUE)
    colnames(eapA) <- car::recode(colnames(eapA), "'pv1'='EAP'")
  }
  cat("Note: To date, only 1pl/2pl dichotomous models are supported.\n"); flush.console()
  if ( is.null(item) & is.null(pdfFolder)) {stop("If ICCs for more than one item should be displayed, please specify an output folder for pdf.\n")}
  if ( !is.null(pdfFolder)) { grDevices::pdf(file = pdfFolder, width = 10, height = 7.5) }
  if ( !is.null ( item ) )  {
    if ( !item %in% it[,"item"]) { stop (paste("Item '",item,"' was not found in 'resultsObj'.\n",sep=""))}
    it <- it[which(it[,"item"] == item),]
  }
  pl  <- by ( data = it, INDICES = it[,c("model", "item")], FUN = function ( i ) {
    xlm <- c(i[["est"]]+2, i[["est"]]-2)
    anf <- -6
    ende<- 6
    x   <- seq ( anf, ende, l = 400)
    y   <- exp( i[["estSlope"]]*x - i[["est"]] ) / (1+exp( i[["estSlope"]]*x - i[["est"]] ))
    plot (x, y, type = "l", main = paste("Item '",as.character(i[["item"]]),"'\n\n",sep=""), xlim = c(-6,6), ylim = c(0,1), xlab = "theta", ylab = "P(X=1)", col = "darkred", cex = 8, lwd = 2)
    graphics::mtext( paste("Model = ",i[["model"]],"  |  Dimension = ",i[["dimension"]], "  |  difficulty = ",round(i[["est"]], digits = 3),"  |  Infit = ",round(i[["infit"]], digits = 3),"\n",sep=""))
    eap <- eapA[intersect ( which (eapA[,"dimension"] == i[["dimension"]]) , which (eapA[,"model"] == i[["model"]])),]
    if ( inherits(defineModelObj, "defineMultiple")) {
      woIst<- which ( lapply ( defineModelObj, FUN = function ( g ) {   g[["analysis.name"]] == i[["model"]] }) == TRUE)
      stopifnot(length(woIst) == 1)
      dat  <-defineModelObj[[woIst]][["daten"]]
    }  else  {
      dat  <- defineModelObj[["daten"]]
    }
    id  <- unique(resultsObj[intersect(which(resultsObj[,"type"] == "tech"), which(resultsObj[,"par"] == "ID")),"derived.par"])
    stopifnot(length(id)==1)
    prbs<- na.omit ( merge ( dat[,c( "ID", as.character(i[["item"]]))], eap[,c( id, "EAP")], by.x = "ID", by.y = id))
    anz <- round ( nrow(prbs) / personsPerGroup ) + 1
    if ( anz < 3 ) { anz <- 3 }
    if ( anz > smooth) { anz <- round(smooth)}
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
    graphics::lines ( x = matr[,"mw"], y = matr[,"lh"], col = "blue", lty = 3, lwd = 3) } )
  if ( !is.null(pdfFolder)) { grDevices::dev.off() } }
