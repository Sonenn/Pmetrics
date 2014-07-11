#' \code{PMBuild} will ensure all dependent packages are installed and compile
#' Fortran source code for permanent Pmetrics modules
#'
#' @title Build Pmetrics
#' @author Michael Neely

PMbuild <- function(){
  
  #load necessary packages
  packages <- packageDescription("Pmetrics")$Suggests
  packages <- gsub("\n","",packages)
  packages <- unlist(strsplit(packages,","))
  cat("\nChecking for required packages...\n")
  for(i in packages){
    if(system.file(package=i)==""){
      #tempf? fix for MTSKNN
      if(i=="MTSKNN") {
        install.packages(i,repos="http://www.lapk.org/software/Pmetrics/MTSKNN",quiet=T)
        next
      }
      if(getOption("repos")[1]=="") {setRepositories()}
      install.packages(i,repos=getOption("repos"),dependencies=T)
    }
    
  }
  currwd <- getwd()
  OS <- getOS()
  compiler <- PMFortranConfig()
  sourcedir <- system.file("code",package="Pmetrics")
  destdir <- switch(OS,"~/.config/Pmetrics/compiledFortran",
                    paste(Sys.getenv("APPDATA"),"\\Pmetrics\\compiledFortran",sep=""),
                    "~/.config/Pmetrics/compiledFortran")
  #remove old files if present
  oldfiles <- Sys.glob(paste(destdir,"*.o",sep="/"))
  if(length(oldfiles)>0) {file.remove(oldfiles)}
  #compile new files
  setwd(sourcedir)
  if(!file.exists(destdir)) dir.create(destdir,showWarnings=F)
  PMfiles <- data.frame(filename=as.character(c("NPprep","NPeng","ITprep","ITeng","ITerr","SIMeng")))
  PMfiles$path <- sapply(PMfiles$filename,function(x) 
    shQuote(list.files(getwd(),pattern=as.character(x))))
  
  
  for(i in 1:nrow(PMfiles)){
    cat(paste("\nCompiling ",i," of ",nrow(PMfiles),": ",PMfiles$filename[i],"...",sep=""))
    flush.console()
    command <- sub("<exec>",paste(PMfiles$filename[i],".o -c",sep=""),compiler)
    command <- sub("<files>",PMfiles$path[i],command)
    fortstatus <- suppressWarnings(system(command,intern=T,ignore.stderr=F))
    if(!is.null(attr(fortstatus,"status"))){
      unlink(switch(OS,"~/.config/Pmetrics",
                    paste(Sys.getenv("APPDATA"),"\\Pmetrics",sep=""),
                    "~/.config/Pmetrics"),recursive=T)
      stop(paste("\nThere was an error compiling ",PMfiles$filename[i],".\nDid you select the right fortran compiler?  If yes, try reinstalling fortran.\nFor gfortran, log into www.lapk.org and access system-specific tips on the Pmetrics installation page (step 5).\n",sep=""))
    }

  }
  cat("\nAll packages installed and permanent Fortran modules compiled.\n")
  flush.console()
  invisible(file.copy(from=Sys.glob("*.o"),to=destdir))
  invisible(file.remove(Sys.glob("*.o")))
  fort <- paste(system.file("config",package="Pmetrics"),"newFort.txt",sep="/")
  writeLines("0",fort) #reset to zero
  setwd(currwd)
}

