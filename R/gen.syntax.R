### called by defineModel()

gen.syntax     <- function(Name,daten, all.Names, namen.all.hg = NULL, all.hg.char = NULL, var.char, model = NULL, anchored, constraints=c("cases","none","items"), pfad=NULL, Title=NULL,n.plausible=5,std.err=c("quick","full","none"), compute.fit ,
           distribution=c("normal","discrete"), method=c("gauss", "quadrature", "montecarlo"), n.iterations=200, nodes=NULL, p.nodes=2000, f.nodes=2000, converge=0.001,deviancechange=0.0001, equivalence.table=c("wle","mle","NULL"), use.letters=use.letters, model.statement=model.statement, conquest.folder = NULL, allowAllScoresEverywhere,
           seed , export = list(logfile = TRUE, systemfile = FALSE, history = TRUE, covariance = TRUE, reg_coefficients = TRUE, designmatrix = FALSE) )  {
   if(is.null(anchored)) {anchored <- FALSE} else {anchored <- TRUE}
   export.default <- list(logfile = TRUE, systemfile = FALSE, history = TRUE, covariance = TRUE, reg_coefficients = TRUE, designmatrix = FALSE)
   mustersyntax <- c("title = ####hier.title.einfuegen####;",
   "export logfile >> ####hier.name.einfuegen####.log;",
   "datafile ####hier.name.einfuegen####.dat;",
   "Format pid ####hier.id.einfuegen####",
   "group",
   "codes ####hier.erlaubte.codes.einfuegen####;",
   "labels  << ####hier.name.einfuegen####.lab;",
   "import anchor_parameters << ####hier.name.einfuegen####.ank;",
   "caseweight",
   "set constraints=####hier.constraints.einfuegen####;",
   "set warnings=no,update=yes,n_plausible=####hier.anzahl.pv.einfuegen####,p_nodes=####hier.anzahl.p.nodes.einfuegen####,f_nodes=####hier.anzahl.f.nodes.einfuegen####;",
   "set seed=####hier.seed.einfuegen####;",
   "export par    >> ####hier.name.einfuegen####.prm;",
   "regression",
   "model ####hier.model.statement.einfuegen####;",
   "estimate ! fit=####hier.fitberechnen.einfuegen####,method=####hier.method.einfuegen####,iter=####hier.anzahl.iterations.einfuegen####,nodes=####hier.anzahl.nodes.einfuegen####,converge=####hier.converge.einfuegen####,deviancechange=####hier.deviancechange.einfuegen####,stderr=####hier.std.err.einfuegen####,distribution=####hier.distribution.einfuegen####;",
   "Itanal >> ####hier.name.einfuegen####.itn;",
   "show cases! estimates=latent >> ####hier.name.einfuegen####.pvl;",
   "show cases! estimate=wle >> ####hier.name.einfuegen####.wle;",
   "equivalence ####hier.equivalence.table.einfuegen#### >> ####hier.name.einfuegen####.equ;",
   "show >> ####hier.name.einfuegen####.shw;",
   "export history >> ####hier.name.einfuegen####.his;",
					 "export covariance >> ####hier.name.einfuegen####.cov;",
					 "export reg_coefficients >> ####hier.name.einfuegen####.reg;",
					 "export designmatrix >> ####hier.name.einfuegen####.mat;",
   "put >> ####hier.name.einfuegen####.cqs;  /* export systemfile */",
   "descriptives !estimates=pv >> ####hier.name.einfuegen####_pvl.dsc;",
   "descriptives !estimates=wle >> ####hier.name.einfuegen####_wle.dsc;",
   "quit;")
   if(is.null(Title))   {                                       ### wenn kein Titel gesetzt, erstelle ihn aus Sys.getenv()
      all.inf  <- Sys.getenv()
      Title    <- paste("Analysis name: ",Name, ", User: ",all.inf["USERNAME"],", Computername: ",all.inf["COMPUTERNAME"],", ", R.version$version.string , ", Time: ",date(),sep="")}
   converge <- paste("0",substring(as.character(converge+1),2),sep="")
   deviancechange <- paste("0",substring(as.character(deviancechange+1),2),sep="")
   syntax    <- gsub("####hier.title.einfuegen####",Title,mustersyntax)
   if(is.null(n.plausible))   {n.plausible <- 0}  ; if(is.na(n.plausible))     {n.plausible <- 0}
   if(n.plausible == 0 )     {                                  ### wenn Anzahl PVs = 0 oder NULL, loesche Statement; andernfalls: setze Anzahl zu ziehender PVs ein!
      syntax    <- gsub("n_plausible=####hier.anzahl.pv.einfuegen####,","",syntax) } else {
      syntax    <- gsub("####hier.anzahl.pv.einfuegen####",n.plausible,syntax)
   }
   syntax    <- gsub("####hier.name.einfuegen####",Name,syntax)
   ID.char   <- max(as.numeric(names(table(nchar(daten[,"ID"])))))
   syntax    <- gsub("####hier.id.einfuegen####",paste("1-",as.character(ID.char)," ",sep="" ) ,syntax)
   syntax    <- gsub("####hier.anzahl.iterations.einfuegen####",n.iterations,syntax)
   syntax    <- gsub("####hier.anzahl.p.nodes.einfuegen####",p.nodes,syntax)
   syntax    <- gsub("####hier.anzahl.f.nodes.einfuegen####",f.nodes,syntax)
   syntax    <- gsub("####hier.converge.einfuegen####",converge,syntax)
   syntax    <- gsub("####hier.deviancechange.einfuegen####",deviancechange,syntax)
   if(!is.null(seed)) {syntax    <- gsub("####hier.seed.einfuegen####",seed,syntax)}
   syntax    <- gsub("####hier.constraints.einfuegen####",match.arg(constraints),syntax)
   compute.fit  <- if(compute.fit == TRUE ) compute.fit <- "yes" else compute.fit <- "no"
   syntax    <- gsub("####hier.fitberechnen.einfuegen####",compute.fit,syntax)
   syntax    <- gsub("####hier.anzahl.nodes.einfuegen####",nodes,syntax)
   syntax    <- gsub("####hier.std.err.einfuegen####",match.arg(std.err),syntax)
   syntax    <- gsub("####hier.distribution.einfuegen####",match.arg(distribution),syntax)
   syntax    <- gsub("####hier.equivalence.table.einfuegen####",match.arg(equivalence.table),syntax)
   syntax    <- gsub("####hier.model.statement.einfuegen####",tolower(model.statement),syntax)
   erlaubte.codes <- paste(gsub("_","",sort(gsub(" ","_",formatC(names(eatTools::tableUnlist(daten[, all.Names[["variablen"]], drop=FALSE ])),width=var.char)),decreasing=TRUE)),collapse=",")
   syntax    <- gsub("####hier.erlaubte.codes.einfuegen####",erlaubte.codes, syntax )
   ind       <- grep("Format pid",syntax)
   beginn    <- NULL                                            ### setze "beginn" auf NULL. Wenn DIF-Variablen spezifiziert sind, wird "beginn" bereits
   if(length(namen.all.hg)>0)    {                              ### untere Zeile: wieviele "character" haben Hintergrundvariablen?
     all.hg.char.kontroll <- all.hg.char
     all.hg.char <- sapply(namen.all.hg, FUN=function(ii) {max(nchar(as.character(na.omit(daten[,ii]))))})
     stopifnot(all(all.hg.char == all.hg.char.kontroll))        ### Trage nun die Spalten in das Format-Statement ein: Fuer ALLE expliziten Variablen
     for (ii in 1:length(namen.all.hg))  {
          if(is.null(beginn)) {beginn <- ID.char+1}
          ende   <- beginn-1+all.hg.char[ii]
          if (beginn != ende) {syntax[ind] <- paste(syntax[ind],namen.all.hg[ii], " ", beginn,"-",ende," ",sep="")}
          if (beginn == ende) {syntax[ind] <- paste(syntax[ind],namen.all.hg[ii], " ", beginn," ",sep="")}
          beginn  <- ende+1 }
   }
   if(length(all.Names[["DIF.var"]])>0)   {                     ### in folgender Schleife ueberschrieben und dann in der Schleife "if(!is.null(HG.var))" ergaenzt, nicht neu geschrieben
      if(model.statement != "item") {
        cat(paste("Caution! DIF variable was specified. Expected model statement is: 'item - ",tolower(all.Names[["DIF.var"]])," + item*",tolower(all.Names[["DIF.var"]]),"'.\n",sep=""))
        cat(paste("However, '",tolower(model.statement),"' will used as 'model statement' to accomplish your will.\n",sep=""))
      }
      if(model.statement == "item") {
         ind.model <- grep("model item", syntax)                ### Aendere model statement
         stopifnot(length(ind.model)==1)
         syntax[ind.model] <- paste("model item - ",paste(tolower(all.Names[["DIF.var"]]),collapse=" - ") ," + ", paste("item*",tolower(all.Names[["DIF.var"]]),collapse=" + "), ";",sep="")
      }
   }
   if(length(all.Names[["HG.var"]])>0)  {
      ind.2   <- grep("^regression$",syntax)
      syntax[ind.2] <- paste(eatTools::crop(paste( c(syntax[ind.2], tolower(all.Names[["HG.var"]])), collapse=" ")),";",sep="")
      if(method == "gauss") {
         cat(paste0("Warning: Gaussian quadrature is only available for models without latent regressors. Use 'Bock-Aiken quadrature' for estimation.\n"))
         method <- "quadrature"                                 ### method muss "quadrature" oder "montecarlo" sein
      }
   }
   syntax    <- gsub("####hier.method.einfuegen####",method,syntax)
   if(length(all.Names[["weight.var"]])>0)  {                   ### Method wird erst hier gesetzt, weil sie davon abhaengt, ob es ein HG-Modell gibt
      ind.4   <- grep("caseweight",syntax)
      syntax[ind.4] <- paste( syntax[ind.4], " ", tolower(all.Names[["weight.var"]]),";",sep="") }
   if(length(all.Names[["group.var"]])>0) {
       ind.3   <- grep("^group$",syntax)
       stopifnot(length(ind.3) == 1)
       syntax[ind.3] <- paste(eatTools::crop(paste( c(syntax[ind.3], tolower(all.Names[["group.var"]])), collapse=" ")),";",sep="")
       ### gebe gruppenspezifische Descriptives
       add.syntax.pv  <- as.vector(sapply(all.Names[["group.var"]], FUN=function(ii) {paste("descriptives !estimates=pv, group=",tolower(ii)," >> ", Name,"_",tolower(ii),"_pvl.dsc;",sep="")} ))
       add.syntax.wle <- as.vector(sapply(all.Names[["group.var"]], FUN=function(ii) {paste("descriptives !estimates=wle, group=",tolower(ii)," >> ", Name,"_",tolower(ii),"_wle.dsc;",sep="")} ))
       ind.3    <- grep("quit",syntax)
       stopifnot(length(ind.3)==1)
       syntax   <- c(syntax[1:(ind.3-1)],add.syntax.pv, add.syntax.wle, syntax[ind.3:length(syntax)]) }
   if(is.null(beginn)) {beginn <- ID.char+1}
   syntax[ind] <- paste(syntax[ind], "responses ",beginn,"-",beginn-1+var.char*ncol(data.frame(daten[,all.Names[["variablen"]]],stringsAsFactors = FALSE)),";",sep="")
   if(var.char>1)  {                                            ### Items haben mehr als eine Spalte Stelligkeit (Conquest-Handbuch, S.177)
      syntax[ind] <- paste(gsub(";","",syntax[ind]), " (a",var.char,");",sep="")}
   score.statement <- .writeScoreStatementMultidim (data=daten, itemCols=all.Names[["variablen"]], qmatrix=model, columnItemNames = 1 ,use.letters=use.letters, allowAllScoresEverywhere = allowAllScoresEverywhere )
   expected.nodes  <- nodes^(ncol(model)-1)
   if(expected.nodes>3500 & method != "montecarlo") {cat(paste("Specified model probably will use ",expected.nodes," nodes. Choosen method ",method," may not appropriate. Recommend to use 'montecarlo' instead.\n",sep=""))}
   ind <- grep("labels ",syntax)
   stopifnot(length(ind)==1)
   syntax <- c(syntax[1:ind],score.statement,syntax[(ind+1):length(syntax)])
   if(length(all.Names[["HG.var"]])==0) {                       ### wenn kein HG-model, loesche entsprechende Syntaxzeilen
      ind.2 <- grep("^regression$",syntax)
      stopifnot(length(ind.2)==1)
      syntax <- syntax[-ind.2]
      ind.3 <- grep("export reg_coefficients",syntax)
      stopifnot(length(ind.3)==1)
      syntax <- syntax[-ind.3] }
   if(length(all.Names[["group.var"]]) ==0) {                   ### wenn keine Gruppen definiert, loesche Statement
      ind.3 <- grep("^group$",syntax)
      stopifnot(length(ind.3)==1)
      syntax <- syntax[-ind.3]}
   if(length(all.Names[["weight.var"]]) ==0) {                  ### wenn keine Gewichte definiert, loesche Statement
      ind.4 <- grep("^caseweight$",syntax)
      stopifnot(length(ind.4)==1)
      syntax <- syntax[-ind.4]}
   if(match.arg(equivalence.table) == "NULL") {                 ### wenn keine Equivalence-Statement definiert, loesche Zeile
      ind.5   <- grep("^equivalence",syntax)
      stopifnot(length(ind.5)==1)
      syntax <- syntax[-ind.5]}
   if(is.null(seed)) {                                          ### wenn keine seed-Statement definiert, loesche Zeile
      ind.7   <- grep("^set seed",syntax)
      stopifnot(length(ind.7)==1)
      syntax <- syntax[-ind.7]}
   if(n.plausible == 0)     {                                   ### wenn Anzahl PVs = 0 oder NULL, loesche Statement
      ind.6   <- grep("^show cases! estimates=latent", syntax)
      stopifnot(length(ind.6) == 1)
      syntax  <- syntax[-ind.6]}
   if(anchored == FALSE) {ind.2 <- grep("anchor_parameter",syntax)# wenn keine ANKER gesetzt, loesche entsprechende Syntaxzeile
                        syntax <- syntax[-ind.2]}
   if(anchored == TRUE)  {ind.2 <- grep("^set constraints",syntax)# wenn ANKER gesetzt, setze constraints auf "none"
                        if(match.arg(constraints) != "none") { cat("Anchorparameter were defined. Set constraints to 'none'.\n")}
                        syntax[ind.2]  <- "set constraints=none;"}
   if(!all(sapply(export, inherits, what="logical"))) {stop("All list elements of argument 'export' have to be of class 'logical'.")}
   export <- as.list(userSpecifiedList ( l = export, l.default = export.default ))
   weg <- names(export[which(export == FALSE)])
   if(length(weg)>0)    {                                       ### hier wird, was nicht exportiert werden soll, aus Syntax geloescht.
      for (ii in seq(along=weg) ) {
           ind.x <- grep(paste("export ", weg[ii], sep=""), syntax)
           stopifnot(length(ind.x) == 1)
           syntax <- syntax[-ind.x]}}
   write(syntax,file.path(pfad,paste(Name,".cqc",sep="")),sep="\n")}


### gen.syntax help functions --------------------------------------------------

.writeScoreStatementMultidim <- function(data, itemCols, qmatrix, columnItemNames = 1 ,columnsDimensions = -1, use.letters=use.letters , allowAllScoresEverywhere) {
  n.dim      <- (1:ncol(qmatrix) )[-columnItemNames]
  stopifnot(length( which( rowSums(qmatrix[,n.dim,drop = FALSE]) == 0))==0)
  if(length(setdiff(names(eatTools::tableUnlist(qmatrix[,-1, drop = FALSE])), c("0","1"))) > 0 )  {
    cat("Found unequal factor loadings for at least one dimension. This will result in a 2PL model.\n")
    for (u in 2:ncol(qmatrix)) {qmatrix[,u] <- as.character(round(qmatrix[,u], digits = 3))}
  }
  stopifnot(all(qmatrix[,1] == itemCols))
  cat(paste("Q matrix specifies ",length(n.dim)," dimension(s).\n",sep=""))
  stopifnot(length(setdiff(colnames(data[,itemCols]),  qmatrix[,columnItemNames]) )==0)
  unique.patter <- qmatrix[which(!duplicated(do.call("paste", qmatrix[,-1, drop = FALSE] ))), -1, drop = FALSE]
  colnames(unique.patter) <- paste("Var",1:ncol(unique.patter), sep="")
  score.matrix  <- data.frame(score=1, unique.patter, matrix(NA, nrow= nrow(unique.patter), ncol=length(itemCols), dimnames=list(NULL, paste("X",1:length(itemCols),sep=""))),stringsAsFactors = FALSE)
  scoreColumns  <- grep("^Var",colnames(score.matrix))
  for (i in 1:length(itemCols))  {
    qmatrix.i    <- qmatrix[qmatrix[,columnItemNames] == itemCols[i],]
    matchRow     <- which(sapply ( 1:nrow(score.matrix) , function(ii) {all ( as.numeric(qmatrix.i[,n.dim]) == as.numeric(score.matrix[ii,scoreColumns])) }))
    stopifnot(length(matchRow) == 1)
    matchColumn  <- min(which(is.na(score.matrix[matchRow,])))
    stopifnot(length(matchColumn) == 1)
    score.matrix[matchRow,matchColumn] <- i
  }
  rowsToDelete <- which(is.na(score.matrix[, max(scoreColumns) + 1]))
  if(length(rowsToDelete)>0) {score.matrix <- score.matrix[-rowsToDelete, ]}
  for (ii in 1:nrow(score.matrix)) {score.matrix[,ii] <- as.character(score.matrix[,ii])}
  score.matrix <- fromMinToMax(dat = data[,itemCols, drop = FALSE], score.matrix = score.matrix, qmatrix = qmatrix, allowAllScoresEverywhere = allowAllScoresEverywhere, use.letters = use.letters)
  kollapse <- lapply(1:nrow(score.matrix), FUN=function(ii) {na.omit(as.numeric(score.matrix[ii,-c(1,scoreColumns)]))})
  kollapse.diff   <- lapply(kollapse,FUN=function(ii) {c(diff(ii),1000)})
  kollapse.ascend <- lapply(kollapse.diff, FUN=function(ii) {unique(c(0, which(ii!=1)))})
  kollapse.string <- list()
  for (a in 1:length(kollapse.ascend))  {
    string   <- list()
    for (i in 2:length(kollapse.ascend[[a]]))   {
      string.i <- unique( c(kollapse[[a]][kollapse.ascend[[a]][i-1]+1], kollapse[[a]][kollapse.ascend[[a]][i]]))
      string.i <- ifelse(length(string.i) == 2,paste(string.i[1],"-",string.i[2],sep=""),as.character(string.i))
      string[[i]] <- string.i
    }
    string <- paste(unlist(string),collapse=", ")
    kollapse.string[[a]] <- string
  }
  control <- lapply(kollapse.string,FUN=function(ii) {eval(parse(text=paste("c(",gsub("-",":",ii),")",sep="")))})
  if (!all(unlist(lapply(1:length(control), FUN=function(ii) {all(kollapse[[ii]] == control[[ii]])})))) {
    cat("Error in creating score statement.\n")
  }
  score.matrix <- data.frame(prefix="score",score.matrix[,c(1,scoreColumns)],items="! items(",kollapse.string=unlist(kollapse.string),suffix=");",stringsAsFactors=F)
  score.statement <- sapply(1:nrow(score.matrix), FUN=function(ii) { paste(score.matrix[ii,],collapse=" ")})
  return(score.statement) }

### ----------------------------------------------------------------------------

fromMinToMax <- function(dat, score.matrix, qmatrix, allowAllScoresEverywhere, use.letters)    {
  all.values <- plyr::alply(as.matrix(score.matrix), .margins = 1, .fun = function(ii) {sort(names(eatTools::tableUnlist(dat[,na.omit(as.numeric(ii[grep("^X", names(ii))])), drop = FALSE])) ) })
  if ( length(all.values) > 1) {
    if ( all ( outer ( all.values, all.values, Vectorize(identical))) == FALSE ) {
      cat(paste("Found different values for dimensions: \n",sep=""))
      for ( u in 1:length(all.values)) {
        cat(paste0("   Dimension ", u, ": values '",paste(all.values[[u]], collapse= "', '"), "' \n"))
      }
      if ( allowAllScoresEverywhere == TRUE ) {
        all.values <- lapply(all.values, FUN = function ( ii ) { sort(unique( unlist ( all.values ) ))})
        cat(paste("Following value definition was done according to 'allowAllScoresEverywhere == TRUE': \n",sep=""))
        for ( u in 1:length(all.values)) {
          cat(paste0("   Dimension ", u, ": values '",paste(all.values[[u]], collapse= "', '"), "' \n"))
        }
      }
    }
  }
  if(use.letters == TRUE )  {minMaxRawdata  <- unlist ( lapply( all.values, FUN = function (ii) {paste("(",paste(LETTERS[which(LETTERS == ii[1]) : which(LETTERS == ii[length(ii)])], collapse=" "),")") } ) ) }
  if(use.letters == FALSE ) {minMaxRawdata  <- unlist ( lapply( all.values, FUN = function (ii) {paste("(",paste(ii[1] : ii[length(ii)],collapse = " "),")")  } ) ) }
  scoring <- unlist( lapply( minMaxRawdata , FUN = function(ii) { paste("(", paste( 0 : (length(unlist(strsplit(ii, " ")))-3), collapse = " "),")")}) )
  stopifnot(length(scoring) == length( minMaxRawdata ), length(scoring) == nrow(score.matrix )  )
  for (i in 1:nrow(score.matrix))    {
    score.matrix$score[i] <- minMaxRawdata[i]
    targetColumns         <- suppressWarnings(intersect ( grep("Var",colnames(score.matrix)), which(as.numeric(score.matrix[i,]) != 0 ) ))
    stopifnot(length(targetColumns) > 0 )
    score.matrix[i,targetColumns]  <- suppressWarnings(unlist(lapply(score.matrix[i,targetColumns], FUN = function ( y ) {paste( "(", paste(as.numeric(y) * na.omit(as.numeric(unlist(strsplit(scoring[i]," ")))), collapse = " "), ")")})))
    nonTargetColumns      <- suppressWarnings(intersect ( grep("Var",colnames(score.matrix)), which(as.numeric(score.matrix[i,]) == 0 ) ))
    if ( length ( nonTargetColumns ) > 0 )    {
      score.matrix[i,nonTargetColumns]  <- "()"
    }
  }
  return(score.matrix)}

### ----------------------------------------------------------------------------

userSpecifiedList <- function ( l, l.default ) {
  if ( !is.null ( names ( l ) ) ) {
    names ( l ) <- match.arg ( names(l) , names(l.default) , several.ok = TRUE )
  } else {
    if(length(l) > length(l.default) )  {
      stop("Length of user-specified list with more elements than default list.\n")
    }
    names ( l ) <- names ( l.default )[seq(along=l)]
  }
  if ( length(l) < length(l.default) ) {
    l <- c ( l , l.default )
    l <- l[!duplicated(names(l))]
    l <- l[match ( names (l) , names(l.default) )]
  }
  return(l)}
