getTrafo <- function(dataBase = "I:/Methoden/10_sonstige Materialien/trafo.rda", mode=c("paper","pc"), grade=c("primary", "secondary"), subject = c("math", "deu", "eng", "frz", "bio", "che", "phy"),
            domain = c("all", "GL", "ZO", "RF", "MS", "GM", "DHW", "ZA", "ME", "FZ", "DZ", "lesen", "hoeren", "ortho", "sg", "CE", "CF", "PE", "PF", "BE", "BF"),
            study = c("bt", "vera")) {
    if(inherits(dataBase, "character")) {load(dataBase)} else {trafo <- dataBase}
    mode   <- match.arg(arg=mode, choices=eval(formals(getTrafo)[["mode"]]))
    grade  <- match.arg(arg=grade, choices=eval(formals(getTrafo)[["grade"]]))
    subject<- match.arg(arg=subject, choices=eval(formals(getTrafo)[["subject"]]), several.ok = TRUE)
    domain <- match.arg(arg=domain, choices=eval(formals(getTrafo)[["domain"]]), several.ok = TRUE)
    if("all" %in% domain) {
       domain <- setdiff(eval(formals(getTrafo)[["domain"]]), "all")
       suppr  <- TRUE
    } else {
       suppr  <- FALSE
    }
    study  <- match.arg(arg=study, choices=eval(formals(getTrafo)[["study"]]))
    target1<- trafo[[mode]][[grade]]
    subj1  <- names(target1)
    mis1   <- setdiff(subject, subj1)
    if(length(mis1)>0) {message(paste0("subject(s) '",paste(mis1, collapse="', '"), "' not included in mode '",mode,"', grade '",grade,"'. These subjects will be ignored."))}
    subj2  <- intersect(subject, subj1)
    if(length(subj2) > 0 ) {
       dom1  <- unique(unlist(lapply(subj2, FUN = function(su) {names(trafo[[mode]][[grade]][[su]][[study]])})))
       dom2  <- setdiff(domain, dom1)
       if(length(dom2)>0 && isFALSE(suppr)) {message(paste0("domain(s) '",paste(dom2, collapse="', '"), "' not included in mode '",mode,"', grade '",grade,"', study '",study,"'. These domains will be ignored."))}
       dom3  <- intersect(domain, dom1)
       if(length(dom3)>0) {
          ret    <- lapply(subj2, FUN = function(su) {
                    extr <- intersect(dom3, names(trafo[[mode]][[grade]][[su]][[study]]))
                    if(length(extr) > 0) {
                       rp <- do.call("rbind", lapply(extr, FUN = function(e) {trafo[[mode]][[grade]][[su]][[study]][[e]][["refPop"]]}))
                       cts<- do.call("c", lapply(extr, FUN = function(e) {trafo[[mode]][[grade]][[su]][[study]][[e]][["cuts"]]}))
                       anc<- do.call(plyr::rbind.fill, lapply(extr, FUN = function(e) {trafo[[mode]][[grade]][[su]][[study]][[e]][["anchor"]]}))
                       inf<- do.call("c", lapply(extr, FUN = function(e) {trafo[[mode]][[grade]][[su]][[study]][[e]][["info"]]}))
                    }
                    return(list(refPop=rp, cuts = cts, anchor = anc, info=inf))})
          refPop <- do.call("rbind", lapply(ret, FUN = function(r) {r[["refPop"]]}))
          cuts   <- do.call("c", lapply(ret, FUN = function(r) {r[["cuts"]]}))
          anchor <- do.call("rbind", lapply(ret, FUN = function(r) {r[["anchor"]]}))
          info   <- do.call("c", lapply(ret, FUN = function(r) {r[["info"]]}))
          return(list(anchor=anchor, refPop = refPop, cuts=cuts, info=info))
       }
    }}
# mode <- "paper"; grade <- "secondary"; subject <- c("bio", "che"); domain = c("BF", "BE", "CF", "PF"); study = "bt"

