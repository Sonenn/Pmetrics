Dopt <- function(run,data,clean=T){
  
  #copy output file if available
  if(!missing(run)){
    if(!file.exists(as.character(run))) {stop(paste(run," not found in the current working directory.\n",sep=""))}
    if(!file.exists(paste(run,"/dopt",sep=""))) dir.create(paste(run,"/dopt",sep=""))
    
    #get output file
    
    outfile <- basename(Sys.glob(paste(run,"/outputs/OUT[0-9]*",sep="")))
    if(length(outfile)==0){ #no output file found
      outfile <- readline(paste("Default output file not found.\nEnter another filename or 'end' to quit: \n")) 
      if (tolower(outfile)=="end") {stop(paste("No output file specified.\n"))}
    }
    if(length(outfile)>1){ #multiple output files found
      cat("The following output files were found.\n")
      for(i in 1:length(outfile)){
        cat(paste(i,". ",outfile[i],"\n",sep=""))
      }
      outfileIndex <- readline(paste("Enter the number of the file you want to use  or 'end' to quit: \n")) 
      if (tolower(outfileIndex)=="end") {stop(paste("No output file specified.\n"))} else {outfile <- outfile[as.numeric(outfileIndex)]}
      
    }
    
    #copy this output file to new /dopt folder
    invisible(file.copy(from=paste(run,"/outputs/",outfile,sep=""),to=paste(run,"/dopt",sep="")))
    
    #now move the data file in dopt folder if supplied  
    if(!missing(data)){
      if(file.exists(data)){
        file.copy(from=data,to=paste(run,"/dopt",sep=""),overwrite=T)
        file.remove(data)
      } else {stop(paste(data,"is not in current directory.\n"))}
    } else {data <- "PATQZPX.001" } #default file name
    
  } else {stop("Please supply a run number.\n")}
  
  
  OS <- getOS()
  #read or define the Fortran compiler
  fortSource <- switch(OS,"~/.config/Pmetrics/compiledFortran",
                       paste(Sys.getenv("APPDATA"),"\\Pmetrics\\compiledFortran",sep=""),
                       "~/.config/Pmetrics/compiledFortran")
  if(!file.exists(fortSource)){
    PMbuild()
  }
  compiler <- PMFortranConfig()
  #choose serial compiliation
  if(length(compiler)==2){
    compiler <- compiler[1]
  }
  if(is.null(compiler)) {cat("\nExecute DOPT after fortran is installed.\n");return(invisible(NULL))}
  
  
  
  enginefiles <- shQuote(normalizePath(list.files(fortSource,pattern="DOeng",full.names=T)))
  enginecompile <- sub("<exec>","do_run",compiler)
  enginecompile <- sub("<files>",paste(enginefiles,"modelqzpx.for"),enginecompile,fixed=T)
  
  
  #make the control file
  if(data=="PATQZPX.001"){
    ControlFile <- c(outfile,   #output file name file
                     "1",       #don't change density file
                     "1",       #don't change  data file
                     "1")       #don't change model file
    
  } else {
    ControlFile <- c(outfile,   #output file name file
                     "1",       #don't change density file
                     "0",       #change  data file
                     data,      #enter data filename
                     "go",      #continue
                     "1")       #don't change model file
  }
  
  
  f <- file(paste(run,"/dopt/DOcontrol",sep=""),"w")
  writeLines(ControlFile,f)
  close(f)
  
  
  wd <- getwd()
  setwd(paste(run,"/dopt",sep=""))
  
  
  #run prep program
  cat("Reading files....\n")
  flush.console()
  invisible(file.copy(from=paste(fortSource,"/DOprep.exe",sep=""),to=getwd()))
  if(OS==1 | OS==3){system("./DOprep.exe MacOSX < DOcontrol",ignore.stdout=T)} 
  if(OS==2) {shell("DOprep.exe DOS < DOcontrol")}
  
  
  
  # RUN  engine
  cat("Calculating requested D-optimal times....\n")
  flush.console()
  
  if(OS==1 | OS==3){
    if(!file.exists("extnum")){
      system("echo 1 > extnum")
    }
    system(enginecompile,ignore.stdout=T)
    system("./do_run",ignore.stdout=T)
  } else {
    if(!file.exists("extnum")){
      system("echo 1 > extnum")
    }
    shell(enginecompile)
    shell("do_run")
  }
  
  #get dopt calc number
  doptNum <- scan("extnum",quiet=T)
  doptNum <- doptNum - 1 
  
  #read results
  doptRes <- readLines(paste("OUTDOPT",sprintf("%04d", doptNum),sep=""))
  lineNum <- grep("DESIGN OF ",doptRes)
  count <- 0
  while(doptRes[lineNum+count]!=""){
    count <- count+1
  }
  doptNum <- count - 1
  gridPts <- grep("FOR GRID PT.",doptRes)
  probPos <- regexpr("[[:digit:]]+\\.*[[:digit:]]*E-*[[:digit:]]*",doptRes[gridPts])
  probs <- sapply(1:length(gridPts),function(x) substr(doptRes[gridPts[x]],start=probPos[x],stop=probPos[x]+attr(probPos,"match.length")[x]))
  eachDopt <- doptRes[unlist(lapply(gridPts,function(x) x+2+(1:doptNum)))]
  allDopt <- data.frame(time=as.numeric(gsub(" ","",eachDopt)))
  allDopt$gridpoint <- rep(1:length(gridPts),each=doptNum)
  allDopt$prob <- as.numeric(rep(gsub(" ","",probs),each=doptNum))
  meanDoptLine <- grep(" THE WEIGHTED MAX.",doptRes)
  meanDopt <- as.numeric(gsub(" ","",doptRes[(meanDoptLine+3):(meanDoptLine+2+doptNum)]))
  results <- list(all=allDopt,mean=meanDopt)
  class(results) <- c("PMdopt","list")
  
  #clean up run
  if(clean){
    toDelete <- Sys.glob(c("*QZPX*","DOprep*","DOcontrol","fort*","do_run*"))
    file.remove(toDelete)
  }
  setwd(wd)
  return(results)
} #end function

