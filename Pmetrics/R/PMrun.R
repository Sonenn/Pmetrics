
.PMrun <- function(type,model,data,run,
                   include,exclude,ode,tol,salt,cycles,
                   indpts,icen,aucint,
                   idelta,prior,xdev,search,
                   auto,intern,silent,overwrite,nocheck){
  
  #check for files
  #copy model file if available
  if(is.numeric(model)){
    modelfile <- suppressWarnings(tryCatch(scan(paste(model,"etc/instr.inx",sep="/"),what="character",quiet=T,skip=4,n=1),error=function(e) NULL))
    modelfile <- Sys.glob(paste(model,"/inputs/",strsplit(modelfile,"\\.")[[1]][1],"*",sep=""))
    if(length(modelfile)>0){
      file.copy(from=modelfile,to=getwd())
      model <- basename(modelfile)
    }
  }
  #make sure model file name is <=8 characters
  if(!FileNameOK(model)) {endNicely(paste("Model file name must be 8 characters or fewer.\n"),model=-99,data)}
  
  while(!file.exists(model)) {
    model <- readline(paste("The model file",shQuote(paste(getwd(),model)),"does not exist.\nEnter another filename or 'end' to quit: \n")) 
    if (tolower(model)=="end") {endNicely(paste("No model file specified.\n"),model=-99,data); break}
  }
  
  #copy wrk files if available
  wrkFlag <- F
  if(is.numeric(data)){
    wrkfiles <- Sys.glob(paste(data,"wrkcopy/*.ZMQ",sep="/"))
    if(length(wrkfiles)>0){
      file.copy(from=wrkfiles,to=getwd())
      wrkFlag <- T
    } else {endNicely(paste("\nNo working copy files found in ",data,"/wrkcopy folder.\n",sep=""),model,data=-99)}
    RFfile <- suppressWarnings(tryCatch(readLines(Sys.glob(paste(data,"outputs/??_RF0001.TXT",sep="/"))),error=function(e) NULL))
    if(length(RFfile)>0){
      datafileName <- tail(RFfile,1)
      file.copy(from=paste(data,"inputs",datafileName,sep="/"),to=getwd())
      data <- datafileName
    } else {endNicely(paste("\nNo RF file found in ",data,"/outputs folder to extract matrix data filename.\n",sep=""),model,data=-99)}
  }
  
  #make sure data file name is <=8 characters
  if(!FileNameOK(data)) {endNicely(paste("Data file name must be 8 characters or fewer.\n"),model,data=-99)}
  
  #ok look for a data file
  while(!file.exists(data)) {
    data <- readline(paste("The data file",shQuote(paste(getwd(),data)),"does not exist.\nEnter another filename or 'end' to quit: \n")) 
    if (tolower(data)=="end") {endNicely(paste("No data file specified.\n"),model,data=-99); break}
  }
  
  #get information from datafile
  dataFile <- PMreadMatrix(data,quiet=T)
  
  #check for errors in data if nocheck=T
  if(nocheck){
    err <- PMcheck(dataFile,quiet=T)
    if(attr(err,"error")==-1){
      endNicely("\nThere are errors in your data file.  Run PMcheck.\n",model,data)
    }
  }
  
  dataoffset <- 2*as.numeric("addl" %in% names(dataFile))
  ncov <- ncol(dataFile)-(12+dataoffset)
  if(ncov>0) {covnames <- names(dataFile)[(13+dataoffset):ncol(dataFile)]} else {covnames <- NA}
  numeqt <- max(dataFile$outeq,na.rm=T)
  id <- unique(dataFile$id)
  nsubtot <- length(id)
  
  if(is.null(include)) include <- id
  if(!is.null(exclude)) include <- id[!id %in% exclude]
  
  nsub <- length(include)
  activesub <- c(which(id %in% include),0)
  ndrug <- max(dataFile$input,na.rm=T)
  if(is.null(salt)) {
    salt <- rep(1,ndrug)
  } else {
    if(length(salt)!=length(ndrug)) endNicely("\nSalt fraction length must equal number of drugs.\n",model,data)
  }
  
  #AUC interval, index of gridpts, prior, and MIC for NPAG
  if(type=="NPAG"){
    #AUC interval
    if(is.null(aucint)){
      maxTime <- max(dataFile$time,na.rm=T)
      if(maxTime>24*48) {aucint <- ceiling(maxTime/48)} else {aucint <- 24}
    }
    
    #gridpts
    if(is.null(indpts)) indpts <- -99
    
    #check if prior and if so, get name for instruction file
    if(is.null(prior)) {  #prior not specified
      prior <- -99
      priorString <- 1
    } else { #prior specified, so choose how
      if(inherits(prior,"NPAG")){priorString <- c(0,"prior.txt")} #prior is an NPdata object
      if(is.character(prior)) {priorString <- c(0,prior)}  #prior is the name of a file
      if(is.numeric(prior)){  #prior is a run number
        priorDEN <- Sys.glob(paste(prior,"outputs/DEN*",sep="/"))[1]
        if(length(priorDEN)>0){
          file.copy(from=priorDEN,to=paste(getwd(),"/prior.txt",sep=""))
          prior <- "prior.txt"
          priorString <- c(0,prior)
          
        }
      }
    }
    #MIC
    xmic <- 1
  }
  
  #set fraction of ab range for initial parameter SD for IT2B and ERR
  if(type=="IT2B" | type=="ERR"){
    xsig=0.5
  }
  
  
  #attempt to translate model and data files into separate fortran model  and instruction files
  
  if(type=="NPAG"){
    engine <- list(alg="NP",nsubtot=nsubtot,nsub=nsub,activesub=activesub,ncov=ncov,covnames=covnames,ndrug=ndrug,tol=tol,
                   salt=salt,numeqt=numeqt,cycles=cycles,icen=icen,indpts=indpts,aucint=aucint,idelta=idelta,xmic=xmic,
                   ode=ode,limits=NA,priorString=priorString,wrkFlag=wrkFlag)
  }
  if(type=="IT2B"){
    engine <- list(alg="IT",nsubtot=nsubtot,nsub=nsub,activesub=activesub,ncov=ncov,covnames=covnames,ndrug=ndrug,
                   salt=salt,numeqt=numeqt,cycles=cycles,xsig=xsig,tol=tol,xdev=xdev,indpts=-99,
                   ode=ode,limits=NA,wrkFlag=wrkFlag)
    
  }
  if(type=="ERR"){
    engine <- list(alg="ERR",nsubtot=nsubtot,nsub=nsub,activesub=activesub,ncov=ncov,covnames=covnames,ndrug=ndrug,
                   salt=salt,numeqt=numeqt,cycles=cycles,xsig=xsig,tol=tol,xdev=xdev,indpts=-99,
                   ode=ode,limits=NA,wrkFlag=wrkFlag)
  }
  
  trans <- makeModel(model=model,data=data,engine=engine,write=T,silent=silent)
  
  if(trans$status==-1) endNicely(trans$msg,model,data) #error
  if(trans$status==0){ 
    useOldFortran <- T  #old fortran file model
    auto <- F #disable automatic running
  } 
  if(trans$status==1){ #new model and instruction files made
    useOldFortran <- F
    modelFor <- trans$modelFor
    ptype <- trans$ptype
    if(!wrkFlag) {ctype <- trans$ctype} else {ctype <- 0}
    instr <- "instr.inx"
  }
  
  OS <- getOS()
  
  fortSource <- switch(OS,"~/.config/Pmetrics/compiledFortran",
                       paste(Sys.getenv("APPDATA"),"\\Pmetrics\\compiledFortran",sep=""),
                       "~/.config/Pmetrics/compiledFortran")
  if(!file.exists(fortSource)){
    PMbuild()
  }
  compiler <- PMFortranConfig()
  if(is.null(compiler)) endNicely(paste("\nExecute ",type," run after gfortran is installed.\n",sep=""),
                                  model,data)
  
  #make new ouput directory
  if(is.null(run)){
    olddir <- list.dirs(recursive=F)
    olddir <- olddir[grep("^\\./[[:digit:]]+",olddir)]
    olddir <- sub("^\\./","",olddir)
    if(length(olddir)>0){
      newdir <- as.character(max(as.numeric(olddir))+1)
    } else {newdir <- "1"}
  } else {
    if(!is.numeric(run)) {endNicely("'run' must be numeric.\n",model,data)} else {newdir <- as.character(run)}
  }
  
  if(file.exists(newdir)){
    if(overwrite) {unlink(newdir,recursive=T)} else {endNicely(paste("\n",newdir," exists already.  Set overwrite=T to overwrite.\n",model,data))}
  }
  
  #substitution string for directory separator according to OS  
  rep <- c("/","\\\\","/")[OS]
  
  #change units of ODE tolerance to linear from log
  ode <- c(0,10**ode)
  
  #generate the names of the permanent modules  
  prepfiles <- shQuote(normalizePath(list.files(fortSource,
                                                pattern=switch(type,NPAG="NPprep",IT2B="ITprep",ERR="ITprep"),full.names=T)))
  enginefiles <- shQuote(normalizePath(list.files(fortSource,
                                                  pattern=switch(type,NPAG="NPeng",IT2B="ITeng",ERR="ITerr"),full.names=T)))
  
  #generate names of files that will be created
  prepFileName <- switch(type,NPAG="np_prep",IT2B="it_prep",ERR="err_prep")
  runFileName <- switch(type,NPAG="np_run",IT2B="it_run",ERR="err_run")
  drivFileName <- switch(type,NPAG="npagdriv.f",IT2B="it2bdriv.f",ERR="assdriv.f")
  scriptFileName <- switch(type,NPAG="npscript",IT2B="itscript",ERR="errscript")
  #list of output files
  if(type=="NPAG") {outlist <- c("DEN*","OUT0*","OUTT*","PRTB*","ILOG*","NP_RF*")}
  if(type=="IT2B") {outlist <- c("DENF*","FROM*","LAST*","OUFF*","OUTF0*","IT_RF*")}
  if(type=="ERR") {outlist <- "ASS*"}
  
  #generate the compile statements  
  prepcompile <- sub("<exec>",prepFileName,compiler)
  prepcompile <- sub("<files>",prepfiles,prepcompile,fixed=T)
  enginecompile <- sub("<exec>",runFileName,compiler)
  enginecompile <- sub("<files>",enginefiles,enginecompile,fixed=T)
  #now substitute the file separator for inclusion in the batch file
  prepfiles <- gsub("/",rep,prepfiles)
  enginefiles <- gsub("/",rep,enginefiles)
  
  #initiation statement  
  if (OS==1 | OS==3){  #Mac and Linux
    system(paste("echo '",type," run initiated at",format(Sys.time(),"%H:%M on %Y %b %d."),
                 "\nPREP FILES:",prepfiles,"\nENGINE FILES:",enginefiles,"\n\n' > log.txt"))
    system(prepcompile,ignore.stderr=F)
  } 
  
  if (OS==2){  #Windows
    shell(paste("echo '",type," run initiated at",format(Sys.time(),"%H:%M on %Y %b %d."),
                "\nPREP FILES:",prepfiles,"\nENGINE FILES:",enginefiles,"\n\n'"))
    shell(prepcompile,ignore.stderr=F)        
  }
  
  if(!intern | !auto){ #run as batch file script
    #build the batch file script
    PMscript <- vector("character")
    #format working directory for batch file
    workdir <- gsub("/",rep,getwd())
    #change working directory
    PMscript[getNext(PMscript)] <- paste("cd",shQuote(workdir))
    if (!auto){
      #manual run of prep program
      PMscript[getNext(PMscript)] <- c(paste("./",prepFileName," MacOSX",sep=""),
                                       paste(prepFileName," DOS",sep=""),
                                       paste("./",prepFileName," MacOSX",sep=""))[OS]
      
    } else {
      if(type=="NPAG"){
        #handle the prior for NPAG runs
        if (inherits(prior,"NPAG")){ 
          nvar <- trans$nvar
          prior$ab <- as.matrix(trans$ab)
          prior$popPoints <- makeFinal(prior)$popPoints
          for(i in 1:nvar){
            if(prior$ab[i,1] > min(prior$popPoints[,i])) endNicely(paste("You have changed ",prior$par[i]," so that the minimum range of ",prior$ab[i,1]," is greater than the minimum prior point value of ",min(prior$popPoints[,i]),"\nThis will cause NPAG to crash.\n",sep=""),model,data)
            if(prior$ab[i,2] < max(prior$popPoints[,i])) endNicely(paste("You have changed ",prior$par[i]," so that the maximum range of ",prior$ab[i,2],"is less than the maximum prior point value of ",max(prior$popPoints[,i]),"\nThis will cause NPAG to crash.\n",sep=""),model,data)
          }
          err <- makeDen(prior,F)
          if(err==-1) stop("\nYour NPdata prior object is older and does not contain the number of dimensions in your model.\nRe-run your NPAG analysis with Pmetrics 0.25 or later before bootstrapping.\n")
          prior <- c(0,"prior.txt")
        } else {
          if(prior != -99){
            prior <- c(0,prior) } else {prior <- 1}
        }
      } #end prior handling block for NPAG runs
      
      #make the control file to execute run with instructions
      ControlFile <- c("1",   #using instruction file
                       "instr.inx", #name of instruction file 
                       "1") #extra in case of error    
      f <- file("PMcontrol","w")
      writeLines(ControlFile,f)
      close(f)
      
      #run prep program
      PMscript[getNext(PMscript)] <- c(paste("./",prepFileName," MacOSX < PMcontrol",sep=""),
                                       paste(prepFileName," DOS < PMcontrol",sep=""),
                                       paste("./",prepFileName," MacOSX < PMcontrol",sep=""))[OS]
    }  
    PMscript[getNext(PMscript)] <- "echo 1 > extnum"
    PMscript[getNext(PMscript)] <- "echo go > go"
    PMscript[getNext(PMscript)] <- paste(enginecompile,drivFileName,sep=" ")
    PMscript[getNext(PMscript)] <- c(paste("./",runFileName," < go",sep=""),
                                     paste(runFileName," <go",sep=""),
                                     paste("./",runFileName," < go",sep=""))[OS]
    PMscript[getNext(PMscript)] <- c("echo;echo Cleaning up....;echo","echo. & echo Cleaning up.... & echo.","echo;echo Cleaning up....;echo")[OS]
    PMscript[getNext(PMscript)] <- c("stty -echo","echo off","stty -echo")[OS]
    PMscript[getNext(PMscript)] <- paste("mkdir ",newdir,sep="")
    PMscript[getNext(PMscript)] <- paste("mkdir ",newdir,c("/inputs","\\inputs","/inputs")[OS],sep="")
    PMscript[getNext(PMscript)] <- paste("mkdir ",newdir,c("/outputs","\\outputs","/outputs")[OS],sep="")
    PMscript[getNext(PMscript)] <- paste("mkdir ",newdir,c("/wrkcopy","\\wrkcopy","/wrkcopy")[OS],sep="")
    PMscript[getNext(PMscript)] <- paste("mkdir ",newdir,c("/etc","\\etc","/etc")[OS],sep="")
    
    
    
    #add the name of the data file to the end of the output for NPAG or IT2B
    if(type=="NPAG" | type=="IT2B"){
      PMscript[getNext(PMscript)] <- paste("echo ",data,
                                           " >> ",substr(type,1,2),"_RF0001.TXT",sep="")   
    }
    
    #check to make sure run completed
    if(type=="NPAG" | type=="IT2B"){
      PMscript[getNext(PMscript)] <- c(paste("if [ ! -f ",substr(type,1,2),"_RF0001.TXT ]; then error=true; else error=false; fi",sep=""),
                                       paste("if not exist ",substr(type,1,2),"_RF0001.TXT (set error=1) ELSE (set error=0)",sep=""),
                                       paste("if [ ! -f ",substr(type,1,2),"_RF0001.TXT ]; then error=true; else error=false; fi",sep=""))[OS]
    }
    if(type=="ERR"){
      PMscript[getNext(PMscript)] <- c("if [ ! -f ASS0001 ]; then error=true; else error=false; fi",
                                       "if not exist ASS0001 (set error=1) ELSE (set error=0)",
                                       "if [ ! -f ASS0001 ]; then error=true; else error=false; fi")[OS]
    }
    
    #move output files    
    for (i in 1:6){
      PMscript[getNext(PMscript)] <- paste(c(paste("if [ -f ",outlist[i]," ]; then mv ",sep=""),
                                             paste("if exist ",outlist[i]," move ",sep=""),
                                             paste("if [ -f ",outlist[i]," ]; then mv ",sep=""))[OS],
                                           outlist[i]," ",newdir,c("/outputs; fi","\\outputs","/outputs; fi")[OS],sep="")
    }
    #if error file exists
    PMscript[getNext(PMscript)] <- c(paste("if [ -f ERROR* ]; then mv ERROR* ",newdir,"/outputs; fi",sep=""),
                                     paste("if exist ERROR* move ERROR* ",newdir,"\\outputs",sep=""),
                                     paste("if [ -f ERROR* ]; then mv ERROR* ",newdir,"/outputs; fi",sep=""))[OS] 
    
    if (instr!=-99){
      PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],instr," ",newdir,c("/etc","\\etc","/etc")[OS],sep="")
      PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],"log.txt ",newdir,c("/outputs","\\outputs","/outputs")[OS],sep="")
      PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],"PMcontrol ",newdir,c("/etc","\\etc","/etc")[OS],sep="")
      
    }
    if(!useOldFortran){  #we are using the new model template
      PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],modelFor," ",newdir,c("/etc/","\\etc\\","/etc/")[OS],modelFor,sep="")  #move fortran file to etc
      PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],model," ",newdir,c("/inputs/","\\inputs\\","/inputs/")[OS],model,sep="")  #move template file to inputs
    } else {PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],model," ",newdir,c("/inputs","\\inputs","/inputs")[OS],sep="")}  #using fortran file directly, so move to inputs
    
    PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],"XQZPJ*.ZMQ ",newdir,c("/wrkcopy","\\wrkcopy","/wrkcopy")[OS],sep="")
    PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],"extnum ",newdir,c("/etc","\\etc","/etc")[OS],sep="")
    if(type=="NPAG"){
      PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],"npag*.* ",newdir,c("/etc","\\etc","/etc")[OS],sep="")
      PMscript[getNext(PMscript)] <- paste(c("rm ","erase ","rm ")[OS],"CHMAX*.*",sep="")
      PMscript[getNext(PMscript)] <- c(paste("if [ -f FROM0001 ]; then mv FROM0001 ",newdir,"/inputs; fi",sep=""),
                                       paste("if exist FROM0001 move FROM0001 ",newdir,"\\inputs",sep=""),
                                       paste("if [ -f FROM0001 ]; then mv FROM0001 ",newdir,"/inputs; fi",sep=""))[OS] 
    }
    if(type=="IT2B" | type=="ERR"){
      PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],"it2b*.* ",newdir,c("/etc","\\etc","/etc")[OS],sep="")
      PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],"itas*.* ",newdir,c("/etc","\\etc","/etc")[OS],sep="")
    }
    if(type=="ERR"){
      PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],"assdriv*.* ",newdir,c("/etc","\\etc","/etc")[OS],sep="")
    }
    PMscript[getNext(PMscript)] <- paste(c("rm ","erase ","rm ")[OS],"fort.*",sep="")
    PMscript[getNext(PMscript)] <- paste(c("rm ","erase ","rm ")[OS],"go",sep="")
    PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],prepFileName,"* ",newdir,c("/etc","\\etc","/etc")[OS],sep="")
    PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],runFileName,"* ",newdir,c("/etc","\\etc","/etc")[OS],sep="")
    PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],data," ",newdir,c("/inputs","\\inputs","/inputs")[OS],sep="")
    if(type=="NPAG" && prior[1]==0) PMscript[getNext(PMscript)] <- paste(c("mv ","move ","mv ")[OS],basename(prior[2])," ",newdir,c("/inputs","\\inputs","/inputs")[OS],sep="")
    
    #make report
    reportscript <- paste(normalizePath(Sys.getenv("PmetricsPath"),winslash="/"),"/Pmetrics/report/",
                          switch(type,NPAG="NP",IT2B="IT",ERR="ERR"),"repScript.R",sep="")
    outpath <- c(paste(workdir,"/",newdir,"/outputs",sep=""),
                 paste(workdir,"\\",newdir,"\\outputs",sep=""),
                 paste(workdir,"/",newdir,"/outputs",sep=""))[OS]
    
    #close the error loop    
    PMscript[getNext(PMscript)] <- c("if ! $error ; then ", "if %error% == 0 (","if ! $error ; then ")[OS]
    
    #call report script and then open HTML file
    PMscript[getNext(PMscript)] <- c(paste(normalizePath(R.home("bin"),winslash="/"),"/Rscript ",shQuote(reportscript)," ",shQuote(outpath)," ",icen,sep=""),
                                     paste(shQuote(paste(gsub("/",rep,normalizePath(R.home("bin"),winslash="/")),"\\Rscript",sep=""))," ",shQuote(reportscript)," ",shQuote(outpath)," ",icen,sep=""),
                                     paste(normalizePath(R.home("bin"),winslash="/"),"/Rscript ",shQuote(reportscript)," ",shQuote(outpath)," ",icen,sep=""))[OS]
    PMscript[getNext(PMscript)] <- c(paste("open ",shQuote(paste(gsub("/",rep,outpath),"/",type,"report.html",sep=""))," ; fi",sep=""),
                                     paste("start ",shQuote(paste(type,"Report"))," ",shQuote(paste(gsub("/",rep,outpath),"\\",type,"report.html",sep="")),")",sep=""),
                                     paste("xdg-open ",shQuote(paste(gsub("/",rep,outpath),"/",type,"report.html",sep=""))," ; fi",sep=""))[OS]
    
    #final clean up
    if(OS==1 | OS==3){ #for Mac or Linux
      PMscript[getNext(PMscript)] <- paste("mv ",scriptFileName," ",newdir,"/etc",sep="")
      
    } else {  #for Windows
      PMscript[getNext(PMscript)] <- paste("copy ",scriptFileName,".bat ",newdir,"\\etc",sep="")
      PMscript[getNext(PMscript)] <- paste("echo. & echo Press any key to complete run and close this window... & echo.")
      PMscript[getNext(PMscript)] <- paste("pause > nul")
      PMscript[getNext(PMscript)] <- paste("erase ",scriptFileName,".bat ",sep="")
    }
    
    PMscript <- PMscript[!is.na(PMscript)]
    
    
    f <- file(c(scriptFileName,paste(scriptFileName,".bat",sep=""),scriptFileName)[OS],"w")
    writeLines(PMscript,f)
    close(f)
    
    if (OS==1){ #Mac
      system(paste("chmod +x ",scriptFileName))
      system(paste("open -a Terminal.app ",shQuote(paste(getwd(),"/",scriptFileName,sep="")),sep=""))
    } 
    if (OS==2){ #Windows
      cat(paste("Launch ",scriptFileName,".bat in your working directory to execute the NPAG run.\n",sep=""))
      outpath <- gsub("\\\\","/",outpath)
    }
    if (OS==3){ #Linux
      system(paste("chmod +x ",scriptFileName))
      system(paste("openvt ",shQuote(paste(getwd(),"./",scriptFileName,sep="")),sep=""))
    } 
    
    return(outpath)
    
    
  } else { #run internally, also for servers
    workdir <- getwd()
    if (!auto){  #allow users to see questions
      if(OS==1 | OS==3){system(paste("./",prepFileName," MacOSX",sep=""))} 
    } else { #we do have instructions for prep
      if(type=="NPAG"){  #handle prior for NPAG
        if (inherits(prior,"NPAG")){ 
          nvar <- trans$nvar
          prior$ab <- as.matrix(trans$ab.df)
          prior$popPoints <- makeFinal(prior)$popPoints
          for(i in 1:nvar){
            if(prior$ab[i,1] > min(prior$popPoints[,i])) endNicely(paste("You have changed ",prior$par[i]," so that the minimum range of ",prior$ab[i,1]," is greater than the minimum prior point value of ",min(prior$popPoints[,i]),"\nThis will cause NPAG to crash.\n",sep=""),model,data)
            if(prior$ab[i,2] < max(prior$popPoints[,i])) endNicely(paste("You have changed ",prior$par[i]," so that the maximum range of ",prior$ab[i,2],"is less than the maximum prior point value of ",max(prior$popPoints[,i]),"\nThis will cause NPAG to crash.\n",sep=""),model,data)
          }
          err <- makeDen(prior,F)
          if(err==-1) stop("\nYour NPdata prior object is older and does not contain the number of dimensions in your model.\nRe-run your NPAG analysis with Pmetrics 0.25 or later before bootstrapping.\n")
          prior <- c(0,"prior.txt")
        } else {
          if(prior != -99){
            prior <- c(0,prior) } else {prior <- 1}
        }
      } #end prior handling block
      
      ControlFile <- c("1",   #we have an instruction file
                       "instr.inx")   #name of the instruction file
      
      f <- file("PMcontrol","w")
      writeLines(ControlFile,f)
      close(f)
      
      if(OS==1 | OS==3){system(paste("./",prepFileName," MacOSX < PMcontrol",sep=""))} 
      if(OS==2) {shell(paste(prepFileName," DOS < PMcontrol",sep=""))}
      
    }  
    
    # RUN  engine
    if(OS==1 | OS==3){
      system("echo 1 > extnum")
      system("echo go > go")
      system(paste(enginecompile,drivFileName,sep=" "))
      system(paste("./",runFileName," < go",sep=""))
    } else {
      shell("echo 1 > extnum")
      shell("echo go > go")
      shell(paste(enginecompile,drivFileName,sep=" "))
      shell(paste(runFileName," < go",sep=""))
    }
    
    #CLEAN UP
    dir.create(newdir)
    dir.create(paste(newdir,"/inputs",sep=""))
    dir.create(paste(newdir,"/outputs",sep=""))
    dir.create(paste(newdir,"/wrkcopy",sep=""))
    dir.create(paste(newdir,"/etc",sep=""))
    
    #write data file name to end of NPAG/IT2B output
    if(type=="NPAG" | type=="IT2B"){
      if(length(Sys.glob("??_RF*"))>0){
        pmfile <- file(Sys.glob("??_RF*"),open="a")
        writeLines(data,pmfile)
        close(pmfile)
      }
    }
    
    #move output files
    file.copy(from=Sys.glob(outlist),to=paste(newdir,"/outputs",sep=""))    
    file.remove(Sys.glob(outlist))
    
    if(auto){
      file.copy(from=instr,to=paste(newdir,"/etc",sep=""))
      file.copy(from="log.txt",to=paste(newdir,"/outputs",sep=""))
      file.copy(from="PMcontrol",to=paste(newdir,"/etc",sep="")) 
      file.copy(from=data,to=paste(newdir,"/inputs",sep=""))
      file.remove(instr)
      file.remove("log.txt")
      file.remove("PMcontrol")
      file.remove(data)
    }
    if(type=="NPAG" && prior[1]==0){
      file.copy(from=prior[2],to=paste(newdir,"/inputs",sep=""))
      file.remove(prior[2])
    }
    
    if(!useOldFortran){  #we are using the new model template
      file.copy(from=modelFor,to=paste(newdir,"/etc/",modelFor,sep=""))  #move fortran file to etc
      file.copy(from=model,to=paste(newdir,"/inputs/",model,sep="")) #move template file to inputs
      file.remove(modelFor)         
    } else {
      file.copy(from=model,to=paste(newdir,"/inputs",sep="")) #using fortran file directly, so move to inputs
      file.remove(model)   
    }  
    
    if(length(Sys.glob("CHMAX*.*"))>0){
      file.remove(Sys.glob("CHMAX*.*"))
    }
    
    if(length(Sys.glob("FROM*"))>0) {
      file.copy(from=Sys.glob("FROM*"),to=paste(newdir,"/inputs",sep=""))
      file.remove(Sys.glob("FROM*"))
    }
    if(length(Sys.glob("ERROR*"))>0) {
      file.copy(from=Sys.glob("ERROR*"),to=paste(newdir,"/outputs",sep=""))
      file.remove(Sys.glob("ERROR*"))
    }
    
    file.remove(Sys.glob("fort.*"))
    file.remove("go")
    
    if(type=="NPAG"){
      file.copy(from=Sys.glob("npag*"),to=paste(newdir,"/etc",sep=""))
      file.remove(Sys.glob("npag*"))
    }
    
    if(type=="IT2B" | type=="ERR"){
      file.copy(from=Sys.glob("it2b*"),to=paste(newdir,"/etc",sep=""))
      file.copy(from=Sys.glob("itas*"),to=paste(newdir,"/etc",sep=""))
      file.remove(Sys.glob("it2b*"))
      file.remove(Sys.glob("itas*"))
    }
    
    if(type=="ERR"){
      file.copy(from="assdriv.f",to=paste(newdir,"/etc",sep=""))   
      file.remove("assdriv.f")
    }
    
    
    file.copy(from=Sys.glob("XQZPJ*.ZMQ"),to=paste(newdir,"/wrkcopy",sep=""))
    file.copy(from="extnum",to=paste(newdir,"/etc",sep=""))
    file.copy(from=Sys.glob("*_prep*"),to=paste(newdir,"/etc",sep=""))
    file.copy(from=Sys.glob("*_run*"),to=paste(newdir,"/etc",sep=""))
    file.remove(Sys.glob("XQZPJ*.ZMQ"))
    file.remove("extnum")
    file.remove(Sys.glob("*_prep*"))
    file.remove(Sys.glob("*_run*"))
    
    #make report
    if(type=="NPAG" | type=="IT2B") {PMreport(paste(workdir,newdir,"outputs",sep="/"),icen=icen,type=type)}
    if(type=="ERR") {ERRreport(paste(workdir,newdir,"outputs",sep="/"),icen=icen,type=type)}
    
    #final clean up
    setwd(workdir)
    file.copy(from=Sys.glob("*.*"),to=paste(newdir,"/inputs",sep=""))
    file.remove(Sys.glob("*.*"))
    outpath <- paste(workdir,newdir,"outputs",sep="/")
    
    
    return(outpath)    
    
  }
  
}