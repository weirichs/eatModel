
.onLoad <- function (libname, pkgname) {
     root <- system.file(package = pkgname)
     if ( !file.exists(file.path(root, "exec"))) {dir.create(file.path(root, "exec"))}
     if ( !file.exists( system.file("exec", "console_Feb2007.exe", package = pkgname) )) {
           if ( !file.exists("i:/Methoden/00_conquest_console/console_Feb2007.exe") ) {
               packageStartupMessage("Cannot find conquest executable file. Please choose manually.")
               fname <- file.choose()
           }  else  {
               fname <- "i:/Methoden/00_conquest_console/console_Feb2007.exe"
           }
           foo <- file.copy(from = fname, to = file.path(root, "exec", "console_Feb2007.exe") )
     }
}

