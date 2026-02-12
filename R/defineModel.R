defineModel <- function(dat, items, id, splittedModels = NULL, irtmodel = c("1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", "GPCM.groups", "2PL.groups", "GPCM.design", "3PL"),
               qMatrix=NULL, DIF.var=NULL, HG.var=NULL, group.var=NULL, weight.var=NULL, anchor = NULL, domainCol=NULL, itemCol=NULL, valueCol=NULL,catCol = NULL, check.for.linking = TRUE,
               minNperItem = 50, removeMinNperItem = FALSE, boundary = 6, remove.boundary = FALSE, remove.no.answers = TRUE, remove.no.answersHG = TRUE, remove.missing.items = TRUE, remove.constant.items = TRUE,
               remove.failures = FALSE, remove.vars.DIF.missing = TRUE, remove.vars.DIF.constant = TRUE, verbose=TRUE, software = c("conquest","tam", "mirt"), dir = NULL,
               analysis.name, schooltype.var = NULL, model.statement = "item",  compute.fit = TRUE, pvMethod = c("regular", "bayesian"), fitTamMmlForBayesian = TRUE, n.plausible=5, seed = NULL, conquest.folder=NULL,
               constraints=c("cases","none","items"),std.err=c("quick","full","none"), distribution=c("normal","discrete"), method=c("gauss", "quadrature", "montecarlo", "quasiMontecarlo"),
               n.iterations=2000,nodes=NULL, p.nodes=2000, f.nodes=2000,converge=0.001,deviancechange=0.0001, equivalence.table=c("wle","mle","NULL"), use.letters=FALSE,
               allowAllScoresEverywhere = TRUE, guessMat = NULL, est.slopegroups = NULL, fixSlopeMat = NULL, slopeMatDomainCol=NULL, slopeMatItemCol=NULL, slopeMatValueCol=NULL,
               progress = NULL, Msteps = NULL, increment.factor=1 , fac.oldxsi=0, export = list(logfile = TRUE, systemfile = FALSE, history = TRUE, covariance = TRUE, reg_coefficients = TRUE, designmatrix = FALSE) )   {
       options(warn=1)
       argL <- mget(ls())                                                       ### Argumentenliste erzeugen
     ### soll der Aufruf fuer mehrere Modelle stattfinden? dann ist splittedModels NICHT null.
     ### Sektion 'multiple models handling': jedes Modell einzeln von 'defineModel' aufbereiten lassen
       if(!is.null(splittedModels)) {                                           ### erster Schritt: Hier wird jetzt erstmal nur die bescheuerte Liste aus 'splitModels' aufbereitet (wenn der Nutzer sie verhunzt hat)
         cleared<- cleanifySplittedModels(lst=splittedModels, argL = argL)
     ### jetzt wird argL modifiziert: es gibt so viele Listenelemente, wie modelle per multiplit gerechnet werden sollen
         argL2  <- lapply(cleared[["models.splitted"]], FUN = function (x) {doAufb(x, argL = argL)})
     ### konsoleninformationen erzeugen (noch nicht printen)
         infos  <- lapply(1:length(argL2), FUN = function (model.nr) { generateConsoleInfo(argL = argL2[[model.nr]], x = cleared[["models.splitted"]][[model.nr]])})
     ### alle Modelle abarbeiten, wahrweise single core oder multicore
         if(is.null(splittedModels[["nCores"]]) || splittedModels[["nCores"]] == 1) {
             resAll <- lapply(1:length(argL2), FUN = function (m) {             ### single core
                       cat(infos[[m]], sep="\n")
                       r1 <- defineModelSingle(a = argL2[[m]])
                       return(r1)})
         }  else  {                                                             ### multicore
             doIt   <- function (laufnummer,  argL2 ) {
                       if(!exists("getResults"))  { library(eatModel) }
                       txt <- capture.output (r1 <- defineModelSingle(a = argL2[[laufnummer]]))
                       return(list(txt=txt, r1=r1))}
             beg    <- Sys.time()
             if(splittedModels[["mcPackage"]] == "parallel") {
                cl  <- makeCluster(splittedModels[["nCores"]], type = "SOCK")
             }  else  {
                cl  <- future::makeClusterPSOCK(splittedModels[["nCores"]], verbose=FALSE)
             }                                                                  ### fuer check der Funktion ueber single core:
             resList<- clusterApply(cl = cl, x = 1:length(argL2), argL2=argL2, fun = doIt)
             stopCluster(cl)                                                    ### resList<- lapply(1:length(argL2), FUN = doIt, argL2=argL2)
             resAll <- lapply(resList, FUN = function (l) {l[["r1"]]})
     ### Konsoleninformationen anzeigen
             for(i in 1:length(resList)) {
                cat(infos[[i]], sep="\n")
                cat(resList[[i]][["txt"]], sep="\n")
             }
             cat(paste ( length(argL2), " models were prepared for estimation: ", sep="")); print( Sys.time() - beg, digits = 3)
         }
         attr(resAll, "split") <- splittedModels
         class(resAll) <- c("defineMultiple", "list")
       }  else  {
     ### kein model split
         resAll <- defineModelSingle(a=argL)
       }
       options(warn=0)
       return(resAll) }

defineModelSingle <- function (a) {
     ### assertions
       lapply(a[c("minNperItem", "boundary", "n.iterations", "p.nodes", "f.nodes","converge","deviancechange", "increment.factor" , "fac.oldxsi", "n.plausible")],checkmate::assert_numeric, lower = 0, len = 1)
       checkmate::assert_numeric(a[["nodes"]], lower = 1, null.ok = TRUE, len = 1)
       lapply(a[c("check.for.linking", "removeMinNperItem", "remove.boundary", "remove.no.answers", "remove.no.answersHG", "remove.missing.items", "remove.vars.DIF.missing", "remove.vars.DIF.constant", "verbose", "compute.fit", "fitTamMmlForBayesian", "use.letters", "allowAllScoresEverywhere")],checkmate::assert_logical, len = 1)
       for ( i in names(a)) { assign(i, a[[i]]) }                               ### alle Objekte in a auf den NAMESPACE exportieren
       checkmate::assert_logical(progress, null.ok = TRUE, len = 1)
       if(is.null(progress)) {progress <- software != "tam"}
       dat  <- eatTools::makeDataFrame(dat, name = "dat")
     ### software checken
       software <- match.arg(arg = tolower(software), choices = eval(formals(defineModel)[["software"]]))
       if(software == "conquest" && is.null(a[["conquest.folder"]]) ) {conquest.folder <- identifyConquestFolder() }
       if(software %in% c("conquest", "tam")) {
          irtmodel <- match.arg(irtmodel, choices = eval(formals(defineModel)[["irtmodel"]]))
          if(is.null(Msteps) ) {                                                ### den Default fuer Msteps so setzen wie in TAM
             if ( irtmodel == "3PL" ) { Msteps <- 10 } else { Msteps <- 4 }
          }
       }
       method   <- match.arg(method, choices = eval(formals(defineModel)[["method"]]))
       pvMethod <- match.arg(pvMethod, choices = eval(formals(defineModel)[["pvMethod"]]))
       if(software == "conquest") {
          original.options <- options("scipen")                                 ### lese Option fuer Anzahl der Nachkommastellen
          options(scipen = 20)                                                  ### setze Option fuer Anzahl der Nachkommastellen
          if(is.symbol(a[["analysis.name"]])) {stop("Please specify 'analysis.name' or use 'software = \"tam\"'\n")}
       }  else  {
          if(is.symbol(a[["analysis.name"]])) {analysis.name <- "not_specified"}
       }
       checkmate::assert_character(model.statement, len = 1)
       if(length(items) == 0 ) {stop("Argument 'items' has no elements.\n",sep="")}
       if(length(items) != length(unique(items)) ) {
          cat(paste0("Warning: Item identifier is not unique. Only ",length(unique(items))," unique item identifiers will be used.\n"))
          items <- unique(items)
       }
       if(length(id) != 1 ) {stop("Argument 'id' must be of length 1.\n",sep="")}
     ### wenn model.statement != "item", muessen die zusaetzlichen Variablen im Datensatz sein (ausser z.B. 'step'. Deshalb hier nur eine Warnung, keine Fehlermeldung)
       if(model.statement != "item") {
          vars <- setdiff(eatTools::crop(unlist(strsplit(model.statement, "\\+|-|\\*"))), "item")
          mis  <- which(!vars %in% colnames(dat))
          if(length(mis)>0) {
             cat(paste0("Model statement '",model.statement,"': Variable(s) '",paste(vars[mis], collapse="', '"), "' from 'model.statement' not found in data.\n"))
             vars <- setdiff(vars,vars[mis])                                    ### ueberschreibt vars objekt und loescht die items aus dem model statement, die es nicht im datensatz gibt, raus
          }
       }  else  {
          vars <- NULL
       }
       allVars     <- list(ID = id, variablen=items, DIF.var=DIF.var, HG.var=HG.var, group.var=group.var, weight.var=weight.var, schooltype.var = schooltype.var, add.vars = vars)
       all.Names   <- lapply(allVars, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = dat, variable=ii)})
     ### wenn software = conquest, duerfen variablennamen nicht mehr als 11 Zeichen haben! das heisst, falls doch, werden die items hier umbenannt
       renam       <- NULL                                                      ### initialisieren
       if(software == "conquest") {
          if(max(nchar(all.Names[["variablen"]]))>10) {
             neu   <- paste0("a", as.numeric(as.factor(all.Names[["variablen"]])))
             i     <- 1
             while(length(intersect(neu, colnames(dat))) > 0) {
                i  <- i+1;
                neu<- paste0(letters[i], as.numeric(as.factor(all.Names[["variablen"]])))
             }
             renam <- data.frame(old = all.Names[["variablen"]], new = neu, stringsAsFactors = FALSE)
             all.Names[["variablen"]] <- eatTools::recodeLookup(all.Names[["variablen"]], renam)
             colnames(dat) <- eatTools::recodeLookup(colnames(dat), renam)
          }
       } else {
          if(!is.null(all.Names[["group.var"]])) {
             message("Specifying group.vars is only allowed for software = 'conquest'. 'group.vars' will be ignored.")
             all.Names[["group.var"]] <- NULL
          }
       }
     ### ID-Variable pruefen und ggf. aendern
       dat <- checkID_consistency(dat=dat, allNam=all.Names, software=software)
     ### Verzeichnis ('dir') pruefen oder erzeugen
       dir <- checkDir(dir=dir, software=software)
     ### pruefen, ob es Personen gibt, die weniger als <boundary> items gesehen haben (muss VOR den Konsistenzpruefungen geschehen)
       dat <- checkBoundary(dat=dat, allNam=all.Names, boundary=boundary, remove.boundary=remove.boundary)
     ### Sektion 'explizite Variablennamen ggf. aendern' ###
       subsNam <- .substituteSigns(dat=dat, variable=unlist(all.Names[-c(1:2)]), all.Names = all.Names)
       if(software == "conquest" || !is.null(all.Names[["DIF.var"]])) {
          if(!all(subsNam$old == subsNam$new)) {                                ### Conquest erlaubt keine gross geschriebenen und expliziten Variablennamen, die ein "." oder "_" enthalten
             sn     <- subsNam[which( subsNam$old != subsNam$new),]
             if(nrow(sn) > 4) {toadd <- " (truncated)"} else {toadd <- ""}
             message("'.', '-', and '_' nor upper case letters are allowed in explicit variable names and numbers in DIF variable name. Delete signs from variables names for explicit and DIF variables",toadd,": \n\n", eatTools::print_and_capture (head(sn, n=4), spaces = 5), "\n")
             colnames(dat) <- eatTools::recodeLookup(colnames(dat), sn[,c("old", "new")])
             all.Names     <- lapply(all.Names, FUN = function ( y ) {eatTools::recodeLookup(y, sn[,c("old", "new")]) })
             if(model.statement != "item") {
                cat("    Remove deleted signs from variables names for explicit variables also in the model statement. Please check afterwards for consistency!\n")
                for ( uu in 1:nrow(sn))  {model.statement <- gsub(sn[uu,"old"], sn[uu,"new"], model.statement)}
             }
          }
          if("item" %in% unlist(all.Names[-c(1:2)])) { stop("Conquest does not allow labelling explicit variable(s) with 'Item' or 'item'.\n") }
       }                                                                        ### untere Zeilen: Dif-Variablen und Testitems duerfen sich nicht ueberschneiden
       if(length(intersect(all.Names$DIF.var, all.Names$variablen))>0)    {stop("Test items and DIF variable have to be mutually exclusive.\n")}
       if(length(intersect(all.Names$weight.var, all.Names$variablen))>0) {stop("Test items and weighting variable have to be mutually exclusive.\n")}
       if(length(intersect(all.Names$HG.var, all.Names$variablen))>0)     {stop("Test items and HG variable have to be mutually exclusive.\n")}
       if(length(intersect(all.Names$group.var, all.Names$variablen))>0)  {stop("Test items and group variable have to be mutually exclusive.\n")}
     ### Sektion 'Q matrix ggf. erstellen und auf Konsistenz zu sich selbst und zu den Daten pruefen' ###
       if(is.null(a[["qMatrix"]])) {
          qMatrix <- data.frame ( item = all.Names$variablen, Dim1 = 1, stringsAsFactors = FALSE)
       } else {
          qMatrix <- checkQmatrixConsistency(qMatrix)                           ### pruefe Konsistenz der q-matrix
          notInDat<- setdiff(qMatrix[,1], all.Names$variablen)
          notInQ  <- setdiff( all.Names$variablen , qMatrix[,1])
          if(length(notInDat)>0) {
             cat(paste("Following ", length(notInDat)," item(s) missed in data frame will be removed from Q matrix: \n    ",paste(notInDat,collapse=", "),"\n",sep=""))
             qMatrix <- checkQmatrixConsistency(qMatrix[-match(notInDat, qMatrix[,1]),])
          }
          if(length(notInQ)>0) {
             cat(paste("Following ", length(notInQ)," item(s) missed in Q matrix will be removed from data: \n    ",paste(notInQ,collapse=", "),"\n",sep=""))
          }                                                                     ### Wichtig! Sicherstellen, dass Reihenfolge der Items in Q-Matrix mit Reihenfolge der Items im Data.frame uebereinstimmt!
          all.Names[["variablen"]] <- qMatrix[,1]
       }
       flush.console()
     ### Sektion 'Alle Items auf einfache Konsistenz pruefen'
       cic <- checkItemConsistency(dat=dat, allNam = all.Names, remove.missing.items=remove.missing.items, verbose=verbose, removeMinNperItem=removeMinNperItem, minNperItem=minNperItem, remove.constant.items=remove.constant.items, model.statement=model.statement, software=software)
     ### Sektion 'Hintergrundvariablen auf Konsistenz zu sich selbst und zu den Itemdaten pruefen'. Ausserdem Stelligkeit (Anzahl der benoetigten character) fuer jede Variable herausfinden
       cbc <- checkBGV(allNam = cic[["allNam"]], dat=cic[["dat"]], software=software, remove.no.answersHG=remove.no.answersHG, remove.vars.DIF.missing=remove.vars.DIF.missing, namen.items.weg=cic[["namen.items.weg"]], remove.vars.DIF.constant=remove.vars.DIF.constant)
     ### Sektion 'Itemdatensatz zusammenbauen' (fuer Conquest ggf. mit Buchstaben statt Ziffern)
       if(length(cbc[["namen.items.weg"]])>0)  {
          cat(paste("Remove ",length(unique(cbc[["namen.items.weg"]]))," test item(s) overall.\n",sep=""))
          cbc[["allNam"]]$variablen <- setdiff(cbc[["allNam"]]$variablen, unique(cbc[["namen.items.weg"]]) )
          qMatrix             <- qMatrix[match(cbc[["allNam"]]$variablen, qMatrix[,1]),]
       }
     ### Sektion 'Personen ohne gueltige Werte identifizieren und ggf. loeschen'. Gibt dat, perNA, datL zurueck
       pwvv<- personWithoutValidValues(dat=cbc[["dat"]], allNam=cbc[["allNam"]], remove.no.answers=remove.no.answers)
     ### Sektion 'Summenscores fuer Personen pruefen'
       cpsc<- checkPersonSumScores(datL = pwvv[["datL"]], allNam = cbc[["allNam"]], dat=pwvv[["dat"]], remove.failures=remove.failures)
     ### Sektion 'Verlinkung pruefen'
       if(check.for.linking == TRUE) {                                          ### Dies geschieht auf dem nutzerspezifisch reduzierten/selektierten Datensatz
          linkNaKeep <- checkLink(dataFrame = cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], remove.non.responser = FALSE, verbose = FALSE )
          linkNaOmit <- checkLink(dataFrame = cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], remove.non.responser = TRUE, verbose = FALSE )
          if(linkNaKeep == FALSE & linkNaOmit == TRUE )  {cat("Note: Dataset is not completely linked. This is probably only due to missings on all cases.\n")}
          if(linkNaKeep == TRUE )                        {cat("Dataset is completely linked.\n")}
       }
     ### Sektion 'Anpassung der Methode (gauss, monte carlo) und der nodes'
       met <- adaptMethod(method=method, software=software, nodes=nodes)
     ### Sektion 'Datensaetze softwarespezifisch aufbereiten: Conquest' ###
       if(length(cbc[["namen.all.hg"]])>0) {all.hg.char <- sapply(cbc[["namen.all.hg"]], FUN=function(ii) {max(nchar(as.character(na.omit(cpsc[["dat"]][,ii]))))})} else {all.hg.char <- NULL}
       if(software == "conquest" )   {                                          ### untere Zeile: wieviele character muss ich fuer jedes Item reservieren?
          var.char  <- sapply(cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], FUN=function(ii) {max(nchar(as.character(na.omit(ii))))})
          no.number <- setdiff(1:length(var.char), grep("[[:digit:]]",var.char))
          if(length(no.number)>0) {var.char[no.number] <- 1}                    ### -Inf steht dort, wo nur missings sind, hier soll die Characterbreite auf 1 gesetzt sein
          if(use.letters == TRUE)   {                                           ### sollen Buchstaben statt Ziffern benutzt werden? Dann erfolgt hier Recodierung.
             rec.statement <- paste(0:25,"='",LETTERS,"'",sep="",collapse="; ")
             daten.temp <- cpsc[["dat"]]                                        ### temporaer numerischen Datensatz erstellen, denn wenn letters = TRUE, stehen hier schon Buchstaben drin, und dann kann man keine part-whole correlation mehr bestimmen
             for (i in cbc[["allNam"]][["variablen"]])  {                       ### Warum erst hier? Weil Pruefungen (auf Dichotomitaet etc. vorher stattfinden sollen)
                cpsc[["dat"]][,i] <- car::recode(cpsc[["dat"]][,i], rec.statement)
             }
             var.char <- rep(1,length(cbc[["allNam"]][["variablen"]]))          ### var.char muss nun neu geschrieben werden, da nun alles wieder einstellig ist!
          }
       }
     ### Sektion 'deskriptive Ergebnisse berechnen und durchschleifen' ###
       daten   <- data.frame(ID=as.character(cpsc[["dat"]][,cbc[["allNam"]][["ID"]]]), cpsc[["dat"]][,cbc[["namen.all.hg"]], drop = FALSE], cpsc[["dat"]][,cbc[["allNam"]][["variablen"]], drop = FALSE], stringsAsFactors = FALSE)
       deskRes <- desk.irt(daten = daten, itemspalten = match(cbc[["allNam"]][["variablen"]], colnames(daten)))
       crit    <- which (deskRes[,"valid"] < minNperItem)
       if(length(crit)>0) {
          cat ( paste ( "Following ",length(crit), " items with less than ",minNperItem," item responses:\n",sep=""))
          options(width=1000)
          deskRes[crit,] |> dplyr::select(-tidyselect::any_of(c("item.nr", "Label", "KB", "Codes", "Abs.Freq", "category"))) |> unique() |> print(digits = 3)
       }                                                                        ### diskriminierung kann nicht bestimmt werden, wenn letters = TRUE, weil hier dann bereits buchstaben drin stehen
       if(inherits(try(discrim <- item.diskrim(daten,match(cbc[["allNam"]][["variablen"]], colnames(daten)))  ),"try-error"))  {
          discrim <- item.diskrim(daten.temp,match(cbc[["allNam"]][["variablen"]], colnames(daten.temp)))
          rm(daten.temp)
       }
       if(length ( cbc[["allNam"]][["schooltype.var"]] ) > 0 ) {                ### jetzt ggf. noch schulformspezifische p-Werte, falls gewuenscht
          deskS <- by ( data = cpsc[["dat"]], INDICES = cpsc[["dat"]][, cbc[["allNam"]][["schooltype.var"]] ], FUN = function ( st ) {
                   drst <- desk.irt(daten = st, itemspalten = match(cbc[["allNam"]][["variablen"]], colnames(st)))
                   colnames(drst) <- car::recode (colnames(drst) , paste0("'item.p'='item.p.",st[1,cbc[["allNam"]][["schooltype.var"]]],"'") )
                   return(drst)})
          for(uu in 1:length( deskS) ) {
             matchU <- match(c("item.nr","Label", "KB", "cases", "Missing", "valid", "Codes" , "Abs.Freq"), colnames(deskS[[uu]]))
             stopifnot ( length (which(is.na(matchU))) == 0 , ncol(deskS[[uu]]) - length ( matchU) == 2)
             deskRes <- merge ( deskRes, deskS[[uu]][,-matchU], by = "item.name", all = TRUE)
          }
       }
       lab <- data.frame(itemNr = 1:length(cbc[["allNam"]][["variablen"]]), item = cbc[["allNam"]][["variablen"]], stringsAsFactors = FALSE)
       if(!is.null(a[["anchor"]]))  {
          ankFrame <- anker (lab = lab, prm = anchor, qMatrix = qMatrix, domainCol=domainCol, itemCol=itemCol, valueCol=valueCol, catCol=catCol)
       } else {
          ankFrame <- NULL
          if(fitTamMmlForBayesian == FALSE ) {
             cat("   Note: 'anchor' is necessary if 'fitTamMmlForBayesian' is FALSE. Because 'anchor' is NULL, 'fitTamMmlForBayesian' is set to be TRUE now.\n")
             fitTamMmlForBayesian <- TRUE
          }
       }
       if(software == "conquest" )   {
          daten$ID <- gsub ( " ", "0", formatC(daten$ID, width=max(as.numeric(names(table(nchar(daten$ID)))))) )
          fixed.width <- c(as.numeric(names(table(nchar(daten[,"ID"])))), all.hg.char, rep(max(var.char),length(var.char)))
     ### erstmal testen, ob die Characterzahl wirklich einheitlich ist ... datensatz wird dazu nicht auf festplatte geschrieben
          txt  <-  capture.output ( gdata::write.fwf(daten , colnames = FALSE,rownames = FALSE, sep="",quote = FALSE,na=".", width=fixed.width))
          stopifnot(length(table(nchar(txt)))==1)                               ### Check: hat der Resultdatensatz eine einheitliche Spaltenanzahl? Muss unbedingt sein!
          rm(txt)                                                               ### Speicher sparen
          gdata::write.fwf(daten , file.path(dir,paste(analysis.name,".dat",sep="")), colnames = FALSE,rownames = FALSE, sep="",quote = FALSE,na=".", width=fixed.width)
          colnames(lab) <- c("===>","item")                                     ### schreibe Labels!
          write.table(lab,file.path(dir,paste(analysis.name,".lab",sep="")),col.names = TRUE,row.names = FALSE, dec = ",", sep = " ", quote = FALSE)
          batch <- paste( normalize.path(conquest.folder),paste(analysis.name,".cqc",sep=""), sep=" ")
          write(batch, file.path(dir,paste(analysis.name,".bat",sep="")))
          foo <- gen.syntax(Name=analysis.name, daten=daten, all.Names = cbc[["allNam"]], namen.all.hg = cbc[["namen.all.hg"]], all.hg.char = all.hg.char, var.char= max(var.char), model=qMatrix, anchored=anchor, pfad=dir, n.plausible=n.plausible, compute.fit = compute.fit,
                 constraints=constraints, std.err=std.err, distribution=distribution, method=met[["method"]], n.iterations=n.iterations, nodes=met[["nodes"]], p.nodes=p.nodes, f.nodes=f.nodes, converge=converge,deviancechange=deviancechange, equivalence.table=equivalence.table, use.letters=use.letters, model.statement=model.statement, conquest.folder = conquest.folder, allowAllScoresEverywhere = allowAllScoresEverywhere, seed = seed, export = export)
          if(!is.null(anchor))  {
             write.table(ankFrame[["resConquest"]], file.path(dir,paste(analysis.name,".ank",sep="")) ,sep=" ", col.names = FALSE, row.names = FALSE, quote = FALSE)
          }
     ### wenn Conquest gewaehlt, dann ggf. Logfile umbenennen, falls es bereits (unter demselben namen) existiert
          if(file.exists( file.path ( dir,  paste(analysis.name,".log",sep=""))) )  {
             cat(paste("Found existing log file '",paste(analysis.name,".log",sep=""), "' in folder '",dir,"'\nConquest analysis will overwrite log file. Original log file will be saved as '",paste(analysis.name,"_old.log'\n",sep=""),sep=""))
             do <- file.rename(from = file.path(dir, paste(analysis.name,".log",sep="")), to = file.path(dir, paste(analysis.name,"_old.log",sep="")))
          }
     ### Sektion 'Rueckgabeobjekt bauen', hier fuer Conquest                    ### setze Optionen wieder in Ausgangszustand
          options(scipen = unlist(original.options)); flush.console()           ### Achtung: setze Konsolenpfade in Hochkommas, da andernfalls keine Leerzeichen in den Ordner- bzw. Dateinamen erlaubt sind!
          ret <- list ( software = software, input = paste("\"", file.path(dir, paste(analysis.name,"cqc",sep=".")), "\"", sep=""), conquest.folder = paste("\"", conquest.folder, "\"", sep=""), dir=dir, analysis.name=analysis.name, model.name = analysis.name, qMatrix=qMatrix, all.Names=cbc[["allNam"]], deskRes = deskRes, discrim = discrim, perNA=pwvv[["perNA"]], per0=cpsc[["per0"]], perA = cpsc[["perA"]], perExHG = cbc[["perExHG"]], itemsExcluded = cbc[["namen.items.weg"]], daten=daten, method=met[["method"]], nodes=met[["nodes"]], p.nodes=p.nodes, f.nodes=f.nodes, renam=renam)
          class(ret) <-  c("defineConquest", "list")
       }
     ### Sektion 'Rueckgabeobjekt fuer tam'
       if(software %in% c("tam", "mirt") )   {
          cat(paste("Q matrix specifies ",ncol(qMatrix)-1," dimension(s).\n",sep=""))
          est.slopegroups <- prepEstSlopegroupsTAM(esg = est.slopegroups, allNam = cbc[["allNam"]])
     ### Ladungsmatrix fuer 2pl und gpcm
          fixSlopeMatPrep <- prepFixSlopeMatTAM(fsm = fixSlopeMat, allNam = cbc[["allNam"]], qma =  qMatrix, slopeMatDomainCol=slopeMatDomainCol, slopeMatItemCol=slopeMatItemCol, slopeMatValueCol=slopeMatValueCol, dat=daten, irtmodel=irtmodel)
          fixSlopeMatPrep[["ori"]] <- fixSlopeMat
          guessMat        <- prepGuessMat(guessMat, allNam = fixSlopeMatPrep[["allNam"]])
       }
       if(software == "tam" )   {
          control         <- list ( snodes = met[["snodes"]] , QMC=met[["QMC"]], convD = deviancechange ,conv = converge , convM = .0001 , Msteps = Msteps , maxiter = n.iterations, max.increment = 1 ,
                                  min.variance = .001 , progress = progress , ridge=0 , seed = seed , xsi.start0=FALSE,  increment.factor=increment.factor , fac.oldxsi= fac.oldxsi)
          if ( !is.null(met[["nodes"]])) { control$nodes <- met[["nodes"]] }
          ret     <- list ( software = software, constraint = match.arg(constraints, choices = eval(formals(defineModel)[["constraints"]])) , qMatrix=qMatrix, anchor=list(ank = ankFrame[["resTam"]], allNam = cbc[["allNam"]]),
                            all.Names=fixSlopeMatPrep[["allNam"]], daten=daten, irtmodel=fixSlopeMatPrep[["irtmodel"]], est.slopegroups = est.slopegroups[["esgNm"]], guessMat=guessMat, control = control,
                            n.plausible=n.plausible, dir = dir, analysis.name=analysis.name, deskRes = deskRes, discrim = discrim, perNA=pwvv[["perNA"]], per0=cpsc[["per0"]], perA = cpsc[["perA"]],
                            perExHG = cbc[["perExHG"]], itemsExcluded = cbc[["namen.items.weg"]], fixSlopeMat = fixSlopeMatPrep[["slopMat"]], estVar = fixSlopeMatPrep[["estVar"]], pvMethod = pvMethod,  fitTamMmlForBayesian=fitTamMmlForBayesian)
          class(ret) <-  c("defineTam", "list")
       }
       if(software == "mirt" )   {
          irtmodel <- prepItemTypeMirt(irtmodel = irtmodel, allNam = cbc[["allNam"]], qMatrix=qMatrix)
          ret      <- list ( software = software, qMatrix=qMatrix, allNam = fixSlopeMatPrep[["allNam"]], daten=daten, irtmodel=irtmodel, anchor=list(ank = ankFrame[["resTam"]], allNam = cbc[["allNam"]]), fixSlopeMat = fixSlopeMatPrep,
                            n.plausible=n.plausible, dir = dir, analysis.name=analysis.name, deskRes = deskRes, discrim = discrim, perNA=pwvv[["perNA"]], per0=cpsc[["per0"]], perA = cpsc[["perA"]],
                            perExHG = cbc[["perExHG"]], itemsExcluded = cbc[["namen.items.weg"]], est.slopegroups = est.slopegroups, progress=progress)
          class(ret) <-  c("defineMirt", "list")
       }
       return(ret)}

