chooseCores <- function(cores = NULL, GBcore = NULL, max.cores = NULL ) {
             # if(!exists("detectCores"))   {library(parallel)}
             n.cores <- parallel::detectCores()
             if(is.na(n.cores))     {cat("Cannot detect cores. Computer does not seem to be suited for multicore analysis.\n")}
             ram     <- ps::ps_system_memory()[["total"]]
             if(!is.null(cores)) {
                use.cores <- as.integer(cores)
                if(use.cores > n.cores) {warning(paste0(use.cores," cores were desired. Your machine only provides ",n.cores," cores."))}
             } else {
                use.cores <- n.cores
             }
             if(!is.null(GBcore)) {
                ramGB <- ram/(1024^3)
				        if(GBcore > ramGB ) {
                   warning(paste0("Not able to use ",GBcore, " giga byte per core. Only ",round(ramGB,digits = 2)," giga byte found altogether."))
                   GBcore <- ramGB
                }
                use.cores.new <- (ramGB) / GBcore
                if(use.cores.new < use.cores ) {use.cores <- use.cores.new}
             }
			 if (!is.null(max.cores)) {
					use.cores <- min( use.cores, max.cores )
			 } 
			 return( use.cores ) }
