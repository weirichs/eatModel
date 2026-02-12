### contains all functions with names like "get.X"
### mostly called by getConquest-functions or by the user


### called by getConquestPVs ---------------------------------------------------

get.plausible <- function(file, quiet = FALSE, forConquestResults = FALSE){
  checkmate::assert_file(file)
  lapply(c(quiet, forConquestResults), checkmate::assert_logical, len = 1)
  #
  input           <- scan(file,what="character",sep="\n",quiet=TRUE)
  input           <- strsplit(eatTools::crop(gsub("-"," -",input) ) ," +")
  n.spalten       <- max ( sapply(input,FUN=function(ii){ length(ii) }) )
  input           <- data.frame( matrix( t( sapply(input,FUN=function(ii){ ii[1:n.spalten] }) ),
                                         length(input), byrow = FALSE), stringsAsFactors = FALSE)
  pv.pro.person   <- sum (input[-1,1]==1:(nrow(input)-1) )
  n.person        <- nrow(input)/(pv.pro.person+3)
  weg             <- c(1, as.numeric(sapply(1:n.person, FUN=function(ii){
                        ((pv.pro.person+3)*ii-1):((pv.pro.person+3)*ii+1) })))
  cases           <- input[(1:n.person)*(pv.pro.person+3)-(pv.pro.person+2),1:2]
  input.sel       <- input[-weg,]
  n.dim <- dim(input.sel)[2]-1
  if(quiet == FALSE) {cat(paste(n.person,"persons and",n.dim,"dimensions(s) found.\n"))
    cat(paste(pv.pro.person,"plausible values were drawn for each person on each dimension.\n"))}
  ID              <- input[  (pv.pro.person + 3) *  (1:n.person) - (pv.pro.person + 2) ,2]
  colnames(input.sel) <- c("PV.Nr", paste("dim.",1:(ncol(input.sel)-1),sep=""))
  input.sel[,1]   <- gsub( " ", "0", formatC(input.sel[,1],width = max(nchar(input.sel[,1]))))
  input.sel$ID    <- rep(ID, each = pv.pro.person)
  is.na.ID        <- FALSE
  if(is.na(input.sel$ID[1])) {
    is.na.ID        <- TRUE
    input.sel$ID    <- rep( 1: n.person, each = pv.pro.person)
  }
  input.melt      <- reshape2::melt(input.sel, id.vars = c("ID", "PV.Nr") , stringsAsFactors = FALSE)
  input.melt[,"value"] <- as.numeric(input.melt[,"value"])
  input.wide      <- data.frame( case = gsub(" ", "0",formatC(as.character(1:n.person),
                                                              width = nchar(n.person))),
                                 reshape2::dcast(input.melt, ... ~ variable + PV.Nr),
                                 stringsAsFactors = FALSE)
  colnames(input.wide)[-c(1:2)] <- paste("pv.", paste(rep(1:pv.pro.person,n.dim),
                                                      rep(1:n.dim, each = pv.pro.person),
                                                      sep = "."), sep = "")
  weg.eap      <- (1:n.person)*(pv.pro.person+3) - (pv.pro.person+2)
  input.eap    <- input[setdiff(weg,weg.eap),]
  input.eap    <- na.omit(input.eap[,-ncol(input.eap),drop=FALSE])
  stopifnot(ncol(input.eap) ==  n.dim)
  input.eap    <- lapply(1:n.dim, FUN=function(ii) {matrix(unlist(as.numeric(input.eap[,ii])),
                                                           ncol=2,byrow = TRUE)})
  input.eap    <- do.call("data.frame",input.eap)
  colnames(input.eap) <- paste(rep(c("eap","se.eap"),n.dim),
                               rep(paste("Dim",1:n.dim,sep="."),each=2),sep="_")
  PV           <- data.frame(input.wide,input.eap, stringsAsFactors = FALSE)
  numericColumns <- grep("pv.|eap_|case",colnames(PV))
  if(is.na.ID == TRUE) {PV$ID <- NA}
  for (ii in numericColumns) {PV[,ii] <- as.numeric(as.character(PV[,ii]))  }
  if(  forConquestResults == TRUE ) {
    return(list ( pvWide = PV, pvLong = input.melt, eap = input.eap))
  }  else {
    return(PV)}}

### called by getConquestWles --------------------------------------------------

get.wle <- function(file){
  checkmate::assert_file(file)
  #
  input <- eatTools::crop(scan(file, what = "character", sep = "\n", quiet = TRUE))
  input <- strsplit(input," +")
  n.spalten <- max ( sapply(input,FUN=function(ii){ length(ii) }) )
  n.wle <- floor((n.spalten-1) / 4)
  input <- suppressWarnings(eatTools::asNumericIfPossible(data.frame(
    matrix(t(sapply(input,FUN=function(ii){ ii[1:n.spalten] }) ),length(input),
           byrow = FALSE), stringsAsFactors = FALSE), force.string = FALSE))
  valid <- na.omit(input)
  cat(paste("Found valid WLEs of ", nrow(valid)," person(s) for ", n.wle, " dimension(s).\n",sep=""))
  if(nrow(valid) != nrow(input)){cat(paste("    ", nrow(input)-nrow(valid),
                                           " persons with missings on at least one latent dimension.\n",
                                           sep="")) }
  namen1<- c(rep ( x = c("n.solved", "n.total"), times = n.wle),
             rep(x = c("wle", "std.wle"), times = n.wle))
  namen2<- rep(rep ( paste(".", 1:n.wle, sep=""), each = 2),2)
  colnames(valid)[(ncol(valid)-length(namen2)):1] <- c("ID","case")[1:(ncol(valid)-length(namen2))]
  colnames(valid)[(ncol(valid)-length(namen2)+1):ncol(valid)] <- paste(namen1,namen2,
                                                                       sep="")
  return(valid)}

### called by getConquestResults -----------------------------------------------

### liest Conquest-Outputfiles (*.shw) als R-Objekte ein (siehe auch P:\Aufgabenentwicklung\Grundschule\Daten\Misc\R-Routinen\R2Conquest.R)
### "dif.term" definiert dabei die Variable, nach der ggf. DIF-Analysen ausgelesen werden sollen, z.B. "item*sex". Wird "dif.term" nicht
### spezifiziert, werden keine DIF-Daten ausgegeben.
get.shw <- function(file, dif.term, split.dif = TRUE, abs.dif.bound = 0.6, sig.dif.bound = 0.3, p.value = 0.9) {
            all.output<- list();   all.terms <- NULL                            ### "dif.term" muss nur angegeben werden, wenn DIF-Analysen geschehen sollen.
            input.all <- scan(file,what="character",sep="\n",quiet=TRUE)
            rowToFind <- c("Final Deviance","Total number of estimated parameters")
            rowToFind <- sapply(rowToFind, FUN = function(ii) {                 ### Find the rows indicated in "rowToFind"
                         row.ii <- grep(ii,input.all)                           ### get the parameter of desired rows
                         stopifnot(length(row.ii) == 1)
                         row.ii <- as.numeric(unlist(lapply (strsplit(input.all[row.ii], " +"), FUN=function(ll) {ll[length(ll)]}) ))
                         return(row.ii)})
            ind       <- grep("TERM",input.all)                                 ### Wieviele Tabellen gibt es einzulesen?
            grenzen   <- grep("An asterisk",input.all)
            if(length(ind)==0) {stop(paste("No TERM-statement found in file ",file,".\n",sep=""))}
            for (i in 1:length(ind)) {
                 term <- eatTools::crop(unlist(strsplit(input.all[ind[i]], ":"))[2])
                 cat(paste0("Found TERM ",i,": '",term,"' \n"))
                 all.terms <- c(all.terms,term)                                 ### Dies dient nur dazu, hinterher die Liste mit ausgelesenen Tabellen beschriften zu koennen.
                 bereich <- (ind[i]+6) : (grenzen[i] -2)                        ### Dies der Bereich, der ausgewaehlt werden muss
                 namen   <- gsub("\\^","",c("No.", strsplit(input.all[bereich[1]-2]," +")[[1]][-1]))
                 namen   <- rep(namen, car::recode(namen, "'CI'=2; else=1"))    ### Wenn ein "CI" als Spaltenname erscheint, muessen daraus im R-Dataframe zwei Spalten werden!
                 inp.sel <- gsub("\\(|)|,"," ",eatTools::crop(input.all[bereich]))# Textfile wird reduziert, und voranstehende und abschliessende Leerzeichen werden entfernt. entferne ausserdem Klammern und Kommas (wenns welche gibt)
                 inp.sel <- gsub("\\*    ", "  NA", inp.sel)                    ### hier: gefaehrlich: wenn mittendrin Werte fehlen, wuerde stringsplit eine unterschiedliche Anzahl Elemente je Zeile finden und die fehlenden Elemente stets ans Ende setzen. Fatal!
                 foo        <- strsplit(inp.sel," +")
                 maxColumns <- max(sapply(foo, length))                         ### Gefahr 2: '*' bezeichnet fixierte Parameter, die keinen Standardfehler haben. Manchmal steht aber trotzdem einer da (z.B. in DIF). Ersetzung soll nur stattfinden, wenn mehr als vier Leerzeichen hinterher
                 nDifferentColumns <- length( table(sapply(foo, length)))
                 maxColumns <- which( sapply(foo, FUN=function(ii){ length(ii) == maxColumns  }) )[1]
    ### untere Zeile: WICHTIG! wo stehen in der Zeile mit den meisten nicht fehlenden Werten Leerzeichen?
                 foo.2      <- stringr::str_locate_all(inp.sel[maxColumns], " ")[[1]][,1]
                 foo.3      <- diff(foo.2)                                      ### zeige die Position des letzten Leerzeichens vor einem Nicht-Leerzeichen
                 foo.3      <- foo.2[foo.3 !=1]                                 ### suche nun in jeder Zeile von input.sel: ist das Zeichen zwei Stellen nach foo.3 ein Leerzeichen? Wenn ja: NA!
                 ESTIMATE   <- stringr::str_locate(input.all[ind[i] + 4], "ESTIMATE")[1,1]
                 foo.3      <- foo.3[foo.3>(ESTIMATE-3)]                        ### Achtung: das alles soll aber nur fuer Spalten beginnen, die hinter ESTIMATE stehen! (missraet sonst fuer Produktterme, z.B. "item*sex")
                 if(nDifferentColumns>1) {
                    if(length(foo.3)>0) {                                       ### Und nochmal: das soll NUR geschehen, wenn es in mindestens einer Zeile nicht die vollstaendige (=maximale) Anzahl von Elementen gibt!
                       for (ii in 1:length(inp.sel)) {                          ### also wenn nDifferentColumns groesser als EINS ist (kleiner darf es nicht sein)
                            for (iii in 1:length(foo.3)) {
                                 if(substr( inp.sel[ii], foo.3[iii] + 2 , foo.3[iii] + 2 ) == " ") {inp.sel[ii] <- paste(substr(inp.sel[ii],1,foo.3[iii]), "NA", substring(inp.sel[ii],foo.3[iii]+3) , sep="")}
                            }
                       }
                    }
                    if(length(foo.3)==0) {cat(paste("There seem to be no values in any columns behind 'ESTIMATE'. Check outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep=""))}
                 }
                 inp.sel    <- strsplit(inp.sel," +")
                 if(length(inp.sel[[1]]) == 0 ) {
                    cli::cli_warn(c(paste("There seem to be no valid values associated with term '",all.terms[length(all.terms)],"' in file: '",file,"'.",sep=""), "i"=paste("Term '",all.terms[length(all.terms)], "' will be ignored.", sep="")))
                    all.terms <- all.terms[-i]
                 }
                 if(length(inp.sel[[1]]) > 0 ) {
                    referenzlaenge <- max (sapply( inp.sel, length) )
                    if(referenzlaenge < length(namen) ) {
                       cat(paste("Several columns seem to be empty for term '",all.terms[length(all.terms)],"' in file: '",file,"'.\n",sep=""))
    ### bloeder spezialfall: wenn dif-Analyse mit 'compute.fit=FALSE' gemacht wurde, fehlen die infit-Spalten ... die eigentlich notwendige zusaetzliche spalte 'add.column' wird dann nicht eingefuegt. Finde raus, ob das der Fall ist
                       head <- eatTools::crop(input.all[bereich[1]-2])          ### Ueberschrift
                       leerz<- gregexpr(" ", head)[[1]]                         ### suche: wo beginnt der zweite Block aufeinanderfolgender leerzeichen?
                       leerd<- which ( diff ( leerz) > 1 )[2]
                       vgl  <- length(strsplit ( eatTools::crop(substr(input.all[bereich[1]], 1, leerd)), split = " +")[[1]])
                       if ( vgl == 4 ) {
                            namen <- c(namen[1:2], "add.column1", namen[3:(referenzlaenge-1)])
                       }  else  {
                            referenzlaenge <- length(namen)
                       }
                    }
    ### Ende spezialfall
                    if(referenzlaenge > length(namen) ) {
                       if(referenzlaenge == length(namen) + 1) {
                          cat(paste("There seem to be one more column than columns names. Expect missing column name before 'ESTIMATE'. \nCheck outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep=""))
                          ind.name <- which(namen == "ESTIMATE")
                          namen    <- append(namen, "add.column", after = 3)
                       }
                       if(referenzlaenge >  length(namen) + 1) {
                          cli::cli_warn(c("There seem to be more columns than names for it.", "i"=paste("Check outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"'. \n",sep="")))
                          namen<- c(namen, rep("add.column",referenzlaenge-length(namen) ))
                       }
                    }
                    inp.sel  <- t(sapply(inp.sel, FUN=function(ii){ c(ii, rep(NA,referenzlaenge-length(ii))) }))
                    colnames(inp.sel) <- namen                                  ### untere Zeile: entferne eventuelle Sternchen und wandle in Dataframe um!
                    inp.sel  <- suppressWarnings(eatTools::asNumericIfPossible(data.frame( gsub("NNA", NA, gsub("NA", NA, gsub("\\*","",inp.sel))), stringsAsFactors = FALSE), force.string = FALSE))
    ### estimate und error muessen auf jeden fall numerisch sein
                    results.sel<- data.frame(inp.sel,filename=as.character(file),stringsAsFactors = FALSE)
                    sapply(intersect(c("ESTIMATE", "ERROR"), colnames(results.sel)), FUN = function (colx) {if(any(is.na(as.numeric(results.sel[,colx])))) {cat(paste0("'",colx,"' column in Outputfile for term '",all.terms[length(all.terms)],"' in file: '",file,"' does not seem to be a numeric value. Please check!\n"))}})
                    if(!missing(dif.term)) {                                    ### Der absolute DIF-Wert ist 2 * "Betrag des Gruppenunterschieds". Fuer DIF muessen ZWEI Kriterien erfuellt sein:
                       if(all.terms[length(all.terms)] == dif.term) {           ### Der absolute DIF-Wert muss groesser als 'abs.dif.bound' (z.B. 0.6) und zugleich signifikant groesser als 'sig.dif.bound' (z.B. 0.3) sein
                          cat(paste("Treat '",all.terms[length(all.terms)],"' as DIF TERM.\n",sep=""))
                          results.sel <- data.frame(results.sel,abs.dif = 2*results.sel$ESTIMATE,stringsAsFactors=FALSE)
                          konfNiveau  <- round(100*p.value)                     ### Das bedeutet, fuer Werte groesser 0.6 darf 0.3 NICHT im 90 bzw. 95%-Konfidenzintervall liegen. Nur dann haben wir DIF!
                          results.sel[,paste("KI.",konfNiveau,".u",sep="")] <- results.sel$abs.dif-2*abs(qnorm(0.5*(1-p.value)))*results.sel$ERROR
                          results.sel[,paste("KI.",konfNiveau,".o",sep="")] <- results.sel$abs.dif+2*abs(qnorm(0.5*(1-p.value)))*results.sel$ERROR
                          results.sel[,paste("sig.",konfNiveau,sep="")] <- ifelse(abs(results.sel[,"abs.dif"])>abs.dif.bound & abs(results.sel[,paste("KI.",konfNiveau,".u",sep="")])>sig.dif.bound & abs(results.sel[,paste("KI.",konfNiveau,".o",sep="")])>sig.dif.bound,1,0)
                          results.sel$filename <- file
                          if(split.dif==TRUE) {
                             results.sel <- results.sel[1:(dim(results.sel)[1]/2),]
                             if(dim(results.sel)[1]!=dim(results.sel)[1]) {warning("missing variables in DIF table.")}
                          }
                       }
                    }
                    all.output[[i]] <- results.sel
                 }
              }
              if(!missing(dif.term)) {if(sum(all.terms==dif.term)==0) {cat(paste("Term declarated as DIF: '",dif.term,"' was not found in file: '",file,"'. \n",sep=""))  }}
              names(all.output) <- all.terms
    ### ggf. Regressionsparameter einlesen!
            	regrStart <- grep("REGRESSION COEFFICIENTS", input.all) + 2
              isRegression <- length(regrStart) > 0
            	if ( isRegression)   {
                  regrEnd <- grep("An asterisk next", input.all)
              		regrEnd <- regrEnd[which(regrEnd > regrStart)][1] - 2
              		if ( is.na(regrEnd)) {
              		     warning(paste0("Regression coefficients seems to be missing or corrupted in '",file,"'."))
              		} else {
                    	 regrInput <- eatTools::crop(input.all[regrStart:regrEnd])
                  		 zeileDimensions <- grep("Regression Variable",input.all)
                       stopifnot(length(zeileDimensions) ==1)
                  		 nameDimensions  <- unlist(strsplit(input.all[zeileDimensions], "  +"))[-1]
                  		 regrRows <- grep("CONSTANT",input.all)
                       regrRows <- regrRows[regrRows<=regrEnd][1]
                  		 regrNamen <- unlist(lapply(strsplit(input.all[regrRows:regrEnd],"  +"), FUN=function(ii) {unlist(ii)[1]} ))
                       regrInputSel <- eatTools::crop(input.all[regrRows:regrEnd])
                       regrInputSel <- gsub("\\*","  NA",gsub("\\(|)","",regrInputSel))
                  		 regrInputSel <- unlist( strsplit(regrInputSel," +") )
                  		 nDimensions  <- (length(  regrInputSel ) / length(regrNamen) - 1 )/2
                       cat(paste("Found ",nDimensions," dimension(s): ",paste(nameDimensions,collapse=", "),"\n",sep=""))
                       cat(paste("Found ",length(regrNamen)-1," regressor(s).\n",sep=""))
                       regrInputSel <- data.frame(matrix(regrInputSel, ncol=2*nDimensions+1, byrow=TRUE),stringsAsFactors=FALSE)
                       for (ii in 2:ncol(regrInputSel))  {regrInputSel[,ii] <- as.numeric(regrInputSel[,ii])}
                       colnames(regrInputSel) <- c("reg.var", paste(rep(c("coef","error"),nDimensions), rep(nameDimensions,each=2),sep="_") )
                       regrInputSel$filename <- file
                  		 all.output$regression <- regrInputSel
              		}
              }
    ### Kovarianz-/ Korrelationsmatrix einlesen: schwierig, also Trennen nach ein- vs. mehrdimensional. Eindimensional: zweimal "-----" zwischen Beginn und Ende des COVARIANCE-Statements
              korStart <- grep("COVARIANCE/CORRELATION MATRIX", input.all)
              korEnd   <- grep("An asterisk next", input.all)
              korEnd   <- min(korEnd[korEnd > korStart])
              korStriche <- grep("-----",input.all)
              korStriche <- korStriche[korStriche > korStart & korStriche < korEnd]
              if(length(korStriche) == 2) {                                     ### eindimensional!
                 varRow    <- grep("Variance", input.all)
                 variance  <- as.numeric( unlist( lapply(strsplit(input.all[varRow]," +"), FUN=function(ll) {ll[length(ll)]}) ) )
                 names(variance) <- "variance"
                 all.output$cov.structure <- variance
              }
              if(length(korStriche) > 2) {                                      ### mehrdimensional!
                 bereich     <- input.all[ (min(korStriche) + 1) : (max(korStriche) - 1 ) ]
                 bereich     <- bereich[ -grep("----",bereich)]
                 bereich     <- strsplit(eatTools::crop(bereich),"  +")
                 for (ii in 2:(length(bereich)-1) )  {
                     if(ii <= length(bereich[[ii]]) )  {
                        bereich[[ii]] <- c(bereich[[ii]][1:(ii-1)], NA, bereich[[ii]][ii:length(bereich[[ii]])])
                     }
                     if(ii > length(bereich[[ii]]) )  {
                        bereich[[ii]] <- c(bereich[[ii]][1:(ii-1)], NA)
                     }
                 }
                 bereich.data.frame <- suppressWarnings(eatTools::asNumericIfPossible(data.frame(do.call("rbind", bereich[-1]),stringsAsFactors=FALSE), force.string = FALSE))
                 colnames(bereich.data.frame) <- bereich[[1]]
                 all.output$cov.structure <- bereich.data.frame
              }
            all.output$final.deviance <- rowToFind
    ### Reliabilitaetsindices einlesen
            i1   <- grep("Dimension: \\(Dimension", input.all)
            all.output$reliability <- do.call("rbind", lapply(i1, FUN = function (z ) {
                    stopifnot(substr(input.all[z+3], 2,35) == "WLE Person separation RELIABILITY:")
                    stopifnot(substr(input.all[z+4], 2,20) == "EAP/PV RELIABILITY:")
                    return(data.frame ( dim = eatTools::crop(eatTools::crop(substring(input.all[z], 13)), ")"),
                           wle.rel = as.numeric(eatTools::crop(substring(input.all[z+3], 36))),
                           eap.rel = as.numeric(eatTools::crop(substring(input.all[z+4], 36))),
                           stringsAsFactors = FALSE))}))
            return(all.output)}


### not called -----------------------------------------------------------------

get.prm <- function(file){
  checkmate::assert_file(file)
  #
  input <- scan(file,what="character",sep="\n",quiet=TRUE)
  input <- strsplit( gsub("\\\t"," ",eatTools::crop(input)), "/\\*")
  ret   <- data.frame(do.call("rbind", strsplit(eatTools::crop(unlist(lapply(input, FUN = function ( l ) {l[1]}))), " +")),
                      stringsAsFactors = FALSE)
  nameI <- eatTools::crop(eatTools::removePattern(eatTools::crop(eatTools::crop(unlist(lapply(input, FUN = function ( l ) {l[length(l)]}))),
                                                                 char = "item"), pattern = "\\*/"))
  ret   <- data.frame ( Case= as.numeric(ret[,1]), item = nameI,
                        parameter = as.numeric(ret[,2]), stringsAsFactors = FALSE)
  return(ret)}

### called by getConquestItn ---------------------------------------------------

get.itn <- function(file){
  checkmate::assert_file(file)
  #
  input <- scan(file, what = "character", sep="\n", quiet = TRUE)
  ind.1 <- grep("==========",input)
  items <- grep( "item:", input )
  diff.last <- ind.1[length(ind.1)-1] - items[length(items)] + 4
  items <- cbind(1:length(items),items,c(diff(items),diff.last))
  ind.2 <- gregexpr(":", input[items[,2]])
  ind.3 <- unlist(ind.2)
  ind.3 <- matrix(ind.3,length(ind.2),byrow=T)
  item.namen <- substr(input[items[,2]], ind.3[,dim(ind.3)[2]]+1+nchar(as.character(items[,1])),100)
  item.namen <- gsub(" ","",item.namen)
  item.namen <- gsub("\\)","",item.namen); item.namen <- gsub("\\(","",item.namen)
  if(dim(ind.3)[2]>1)
  {stopifnot(length(table(ind.3[,1]))==1)
    dif.name <- rep(substr(input[items[,2]], 1, ind.3[,1]-1),(items[,3]-11))
    dif.value <- rep(as.numeric(substr(input[items[,2]], ind.3[,1]+1, ind.3[,1]+1)),(items[,3]-11))}
  zeilen <- list(); reihe <- NULL
  for (i in 1:dim(items)[1])
  {zeilen[[i]] <- (items[i,2]+7) : (items[i,2]+ (items[i,3]-5) )
  cases       <- gsub("NA ","NA",input[zeilen[[i]]])
  cases <- gsub("_BIG_ ","NA",cases)
  cases <- gsub("_BIG_","NA",cases)
  if(length(table(sapply(1:length(cases),FUN=function(ii){length(unlist(strsplit(cases[ii]," +"))) }) ) )>1 )
  {cases <- gsub("          ","    NA    ",cases)}
  cases       <- data.frame( matrix ( unlist( strsplit(eatTools::crop(gsub(" +"," ", cases))," ") ),
                                      nrow=length(zeilen[[i]]),byrow=T ) , stringsAsFactors=F)
  ind         <- grep("\\)",cases[1,]); cases[,ind] <- gsub("\\)","",cases[,ind] )
  cases       <- data.frame(cases[,1:(ind-1)],matrix(unlist(strsplit(cases[,6],"\\(")),
                                                     nrow=length(zeilen[[i]]),byrow=T),
                            cases[,-c(1:ind)],stringsAsFactors=F)
  for(jj in 1:ncol(cases)) {cases[,jj] <- as.numeric(cases[,jj])}
  colnames(cases) <- c("Label","Score","Abs.Freq","Rel.Freq","pt.bis","t.value",
                       "p.value",paste(rep(c("PV1.Avg.","PV1.SD."), ((ncol(cases)-7)/2)),
                                       rep(1:((ncol(cases)-7)/2),each=2),sep=""))
  threshold.zeile   <- input[items[i,2]+2]; threshold <- NULL; delta <- NULL
  bereich <- ifelse( (items[i,3]-12)<1,1,(items[i,3]-12))
  if((items[i,3]-12)<1) {cat(paste("Item",i,"hat nur eine Antwortkategorie.\n"))}
  for (j in 1: bereich )
  {threshold  <- c(threshold ,as.numeric(substr(threshold.zeile,  6*j+16,6*j+21)))
  delta      <- c(delta,     as.numeric(substr(input[items[i,2]+3],6*j+13,6*j+18)))}
  while(length(threshold) < nrow(cases)) {threshold <- c(threshold,NA)}
  while(length(delta) < nrow(cases)) {delta <- c(delta,NA)}
  item.p <- NA
  valid.p <- which(is.na(cases$Score))
  if(length(valid.p) == 0)
  {item.p <- cases[which(cases$Score == max(cases$Score)),"Abs.Freq"] / sum(cases$Abs.Freq)}
  sub.reihe   <- data.frame(item.nr=i, item.name=item.namen[i], cases[,1:2],
                            n.valid = sum(cases$Abs.Freq), cases[,3:4], item.p = item.p,
                            diskrim=as.numeric(substr(input[items[i,2]+1],45,55)),
                            cases[,-c(1:4)], threshold, delta, stringsAsFactors=F)
  reihe <- rbind(reihe,sub.reihe)}
  if(dim(ind.3)[2]>1)
  {reihe <- data.frame(dif.name,dif.value,reihe,stringsAsFactors=FALSE)}
  return(reihe)}

### not called -----------------------------------------------------------------

get.dsc <- function(file) {
  checkmate::assert_file(file)
  input  <- scan(file,what="character",sep="\n",quiet=TRUE)
  n.grp  <- grep("Group: ",input)
  grpNams<- unlist( lapply( strsplit(input[n.grp]," ") , function(ll) {paste(ll[-1],collapse=" ")} ) )
  trenn1 <- grep(paste(rep("-", 18), collapse=""), input)
  trenn2 <- grep(paste(rep("\\.", 18), collapse=""), input)
  stopifnot(length(trenn1) == length(trenn2))
  daten  <- lapply(1:(length(trenn1)/2), FUN=function(ii) {
    dat <- data.frame ( do.call("rbind", strsplit(input[(trenn1[2*ii]+1):(trenn2[2*ii-1]-1)]," +")))
    dat <- data.frame(group.name = grpNams[ii], tidyr::unite(dat, col = "X1", colnames(dat)[1:(ncol(dat)-4)], sep=" ")) |> eatTools::asNumericIfPossible(force.string = FALSE) |> suppressWarnings()
    colnames(dat) <- c("group.name","dimension","N","mean","std.dev","variance")
    desc <- data.frame ( do.call("rbind", strsplit(input[(trenn2[2*ii-1]+1):(trenn2[2*ii]-1)]," +")))
    desc <- tidyr::unite(desc, col = "X1", colnames(desc)[1:(ncol(desc)-3)], sep=" ") |> eatTools::asNumericIfPossible(force.string = FALSE) |> suppressWarnings()
    colnames(desc) <- c("dimension","mean","std.dev","variance")
    return(list( single.values=dat, aggregates=desc, ndim = length(grep("Error", desc[,"dimension"])))) } )
  names(daten) <- grpNams
  ndim   <- unique(unlist(lapply(daten, FUN = function(d){d[["ndim"]]})))
  stopifnot(length(ndim)==1)
  names(grpNams) <- rep("i", length(grpNams))
  cli::cli_inform(c(paste0("Found ",length(n.grp)," group(s) and ",ndim," dimension(s) in '",file,"'"), grpNams))
  return(daten)}


### not called -----------------------------------------------------------------

get.equ <- function(file){
  checkmate::assert_file(file)
  #
  input       <- scan(file,what="character",sep="\n",quiet = TRUE)
  dimensionen <- grep("Equivalence Table for",input)
  cat(paste("Found ",length(dimensionen), " dimension(s).\n",sep=""))
  ende        <- grep("================",input)
  ende        <- sapply(dimensionen, FUN=function(ii) {ende[ende>ii][1]})
  tabellen    <- lapply(1:length(dimensionen), FUN=function(ii)
  {part <- eatTools::crop(input[(dimensionen[ii]+6):(ende[ii]-1)])
  part <- data.frame(matrix(as.numeric(unlist(strsplit(part," +"))),ncol=3,
                            byrow=T),stringsAsFactors=F)
  colnames(part) <- c("Score","Estimate","std.error")
  return(part)})
  regr.model  <- grep("The regression model",input)
  item.model  <- grep("The item model",input)
  stopifnot(length(regr.model) == length(item.model))
  name.dimensionen <- unlist( lapply(dimensionen,FUN=function(ii) {
    unlist(lapply(strsplit(input[ii], "\\(|)"),FUN=function(iii){iii[length(iii)]}))}) )
  model       <- lapply(1:length(regr.model), FUN=function(ii) {
    rbind(eatTools::crop(gsub("The regression model:","",input[regr.model[ii]])),
          eatTools::crop(gsub("The item model:","",input[item.model[ii]])) ) })
  model       <- do.call("data.frame",args=list(model,row.names=c("regression.model","item.model"),
                                                stringsAsFactors=F))
  colnames(model) <- name.dimensionen
  tabellen$model.specs <- model
  names(tabellen)[1:length(dimensionen)] <- name.dimensionen
  return(tabellen)}
