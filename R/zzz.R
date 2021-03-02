
.onLoad <- function (libname, pkgname) {
     print(libname)
     print(pkgname)
     root <- system.file(package = "eatModel")
     if ( !file.exists(file.path(root, "exec"))) {dir.create(file.path(root, "exec"))}
     if ( !file.exists( system.file("exec", "console_Feb2007.exe", package = "eatModel") )) {
           if ( !file.exists("i:/Methoden/00_conquest_console/console_Feb2007.exe") ) {
               packageStartupMessage("Cannot find conquest executable file. Please choose manually.")
               fname <- file.choose()
           }  else  {
               fname <- system.file("exec", "console_Feb2007.exe", package = "eatModel")
           }
           foo <- file.copy(from = fname, to = file.path(root, "exec", "console_Feb2007.exe") )
     }
}

