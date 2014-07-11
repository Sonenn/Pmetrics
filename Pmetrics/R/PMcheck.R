#' This function will check a .csv file or a data frame containing a
#' previously loaded .csv file (the output of \code{\link{PMreadMatrix}} for errors
#' which would cause the analysis to fail.  If a model file is provided, and the data
#' file has no errors, it will also check the model file for errors.
#'
#' Either a filename or a data object in memory are accepted as \code{data}.
#' The format of the .csv matrix file is fairly rigid.
#' It must have the following features.  Text is case-sensitive.
#' \itemize{
#'  \item A header in row 1 with the appropriate version, currently \dQuote{POPDATA DEC_11}
#'	\item Column headers in row 2.  These headers are: #ID, EVID, TIME, DUR, DOSE, ADDL, II, INPUT, OUT, OUTEQ,
#' C0, C1, C2, C3.
#'  \item No cell should be empty.  It should either contain a value or \dQuote{.} as a placeholder.
#'	\item Columns after OUTEQ are interpreted as covariates.
#'	\item All subject records must begin with a dose event (EVID=1).
#'	\item All subject records must begin with TIME=0.
#'  \item All dose events (EVID=1) must have entries in ID, EVID, TIME, DUR, DOSE and INPUT.  ADDL and II are optional, but if ADDL is not 0 or
#' missing, then II is mandatory.
#'  \item All observation events (EVID=0) must have entries in ID, EVID, TIME, OUT, OUTEQ.
#'  If an observation is missing, use \emph{-99}; otherwise use a \dQuote{.} as a placeholder
#'  in cells that are not required (e.g. INPUT for an observation event).
#'  \item If covariates are present in the data, there must be an entry for every covariate at time 0 fore each subject.
#'  \item All covariates must be numeric.
#'  \item All times within a subject ID must be monotonically increasing.
#'  \item All subject IDs must be contiguous.
#'  \item All rows must have EVID and TIME values.
#'  \item EVID, TIME, DUR, DOSE, ADDL, II, INPUT, OUT, OUTEQ,
#' C0, C1, C2, C3 must all have numeric entries.
#'  }
#'  
#' To use this function, see the example below.  As another example, assume that 
#' there is a PMdata object called \dQuote{mdata} loaded in memory.
#' Also assume that on rows 1, 10, 150, there are incomplete dose records.
#' 
#' Use \code{err <- PMcheck(mdata)} to run the check.
#' To see the rows in mdata with problems, read the report on the console.
#' In this example, you would use \code{mdata[err$doseComp$results,]}.
#'  
#' You could then try to fix the problem(s) with \code{mdata2 <- PMcheck(mdata,fix=T)}.  Note that we are now returning
#' a PMmatrix data object (hopefully cleaned of errors) rather than the PMerr object returned when \code{fix=FALSE}.
#' Pmetrics handles each of the errors in the following ways.
#' \itemize{
#'  \item If the columns are simply out of order, they will be reordered.  If some are missing, the fix must
#'  be done by the user, i.e. manually.
#'  \item All id and covariate values are truncated to 11 characters.
#'  \item Time=0 observations are deleted.
#'  \item Missing observations are set to -99 (not \dQuote{.}).
#'  \item Incomplete dose records are flagged for the user to fix manually.
#'  \item Incomplete observation records are flagged for the user to fix manually.
#'  \item Subjects without an EVID=1 as first event are flagged for the user to fix manually.
#'  \item Subjects with TIME != 0 as first event have dummy dose=0 events inserted at time 0.
#'  \item Subjects with a missing covariate at time 0 are flagged for the user to fix manually.
#'  \item Non-numeric covariates are converted to numeric (via \link{factor}).
#'  \item Non-ordered times are sorted within a subject if there are no EVID=4 events; otherwise the
#'  user must fix manually.
#'  \item Non-contiguous subject ID rows are combined and sorted if there are no EVID=4 events; otherwise the
#'  user must fix manually.
#'  \item Rows missing an EVID are assigned a value of 0 if DOSE is  missing, 1 otherwise.
#'  \item Rows missing a TIME value are flagged for the user to fix manually.
#'  \item Columns that are non-numeric which must be numeric are flagged for the user to fix manually.
#'  These are EVID, TIME, DUR, DOSE, ADDL, II, INPUT, OUT, OUTEQ, C0, C1, C2, and C3.  
#'  Covariate columns are fixed separately (see above).
#' }
#' 
#' @title Check Pmetrics Inputs for Errors
#' @param data The name of a Pmetrics .csv matrix file in the current working directory,
#' the full path to one not in the current working directory, or a data.frame containing 
#' the output of a previous \code{\link{PMreadMatrix}} command.
#' @param model The filename of a Pmetrics model file in the current working directory.  This parameter is optional.
#' If specified, and the data object has no errors, the model file will be evaluated.
#' @param fix Boolean operator; if \code{TRUE}, Pmetrics will attempt to fix errors in the data file.
#' Default is \code{FALSE}.
#' @param quiet Boolean operator to suppress printed output.  Default is false.
#' @return If \code{fix=TRUE}, then \code{PMcheck} returns a PMmatrix data object which has been
#' cleaned of errors as much as possible, displaying a report on the console.  
#' If \code{fix=FALSE}, then \code{PMcheck} returns a list of objects of class \emph{PMerr}.  Each object is itself a list whose 
#' first object (\code{$msg}) is a character vector with \dQuote{OK} plus a brief description if there is no error, or the error.  
#' The second object (\code{$results}) is a vector of the row numbers that contain that error.
#'  \item{colorder}{The first 14 columns must be named id, evid, time, dur, dose, addl, ii, input, out, outeq, c0, c1, c2, and c3 in that order.} 
#'  \item{maxchar}{All id and covariate values should be less than or equal to 11 characters.}
#'  \item{obsT0}{Time=0 events should not be observations.}
#'  \item{obsMiss}{Missing observations should be -99 (not \dQuote{.}, which is simply a placeholder).}
#'  \item{doseComp}{Make sure all dose records are complete, i.e. contain id, time, evid=1 or 4, duration (0 for bolus), dose, input number.}
#'  \item{obsComp}{Make sure all observation records are complete, i.e. contain id, time, output, and outeq number.}
#'  \item{evid1}{Ensure that each subject's first record is an evid=1.}
#'  \item{T0}{Make sure each subject's first time=0.}
#'  \item{covT0}{Make sure that there is an non-missing entry for each covariate at time=0 for each subject.}
#'  \item{covNumeric}{Ensure that all covariate entries are numeric.}
#'  \item{timeOrder}{Ensure that all times within a subject ID are monotonically increasing.}
#'  \item{contigID}{Ensure that all subject IDs are contiguous.}
#'  \item{missEVID}{Ensure that all rows have an EVID value.}
#'  \item{missTIME}{Ensure that all rows have a TIME value.}
#'  \item{nonNum}{Ensure that all columns which must be numeric are numeric.  These are EVID,
#'  TIME,DUR,DOSE,ADDL,II,INPUT,OUT,OUTEQ,C0,C1,C2,C3.  Covariate columns are checked separately
#'  (see above) but also must be numeric.}



#' @author Michael Neely
#' @seealso \code{\link{PMwriteMatrix}}, \code{\link{PMreadMatrix}}
#' @examples
#' data(PMex3)
#' err <- PMcheck(badData)
#' badData[err$obsLast$results,]
#' badData[err$covT0$results,]
#' goodData <- PMcheck(badData,fix=T)
#' PMcheck(goodData)
#' #you have to fix manually problems which require data entry



PMcheck <- function(data,model,fix=F,quiet=F){
  
  
  #here's the subfunction to check for errors
  errcheck <- function(data2,model,quiet=quiet){
    err <- list(colorder=list(msg="OK - The first 14 columns are appropriately named and ordered.",results=NA),
                maxchar=list(msg="OK - All columns contain entries of 11 or fewer characters.",results=NA),
                obsT0=list(msg="OK - No observations are at time 0.",results=NA),
                obsMiss=list(msg="OK - No missing observations.",results=NA),
                doseComp=list(msg="OK - All dose records are complete.",results=NA),
                obsComp=list(msg="OK - All observation records are complete.",results=NA),
                evid1=list(msg="OK - All subjects have evid=1 as first record.",results=NA),
                T0=list(msg="OK - All subjects have time=0 as first record.",results=NA),
                covT0=list(msg="OK - There are no covariates in the dataset.",results=NA),
                covNumeric=list(msg="OK - there are no covariates in the dataset.",results=NA),
                timeOrder=list(msg="OK - All times are increasing within a subject, given any EVID=4.",results=NA),
                contigID=list(msg="OK - All subject IDs are contiguous.",results=NA),
                missEVID=list(msg="OK - All rows have an EVID value.",results=NA),
                missTIME=list(msg="OK - All rows have a TIME value.",results=NA),
                nonNum=list(msg="OK - All columns that must be numeric are numeric.",results=NA)
    )
    #set initial attribute to 0 for no error
    attr(err,"error") <- 0
    
    # 1 - check to make sure first 14 columns are correct
    t <- tolower(names(data2))
    if(any(!c("id","time","evid") %in% t)){
      #must at least have id, evid, and time columns to proceed with the check
      return(-1)
    }
    if(length(t)<14 | any(!c("id","evid","time","dur","dose","addl","ii","input","out","outeq","c0","c1","c2","c3") %in% t)) {
      err$colorder$msg <- "FAIL - The first 14 columns must be named id, evid, time, dur, dose, addl, ii, input, out, outeq, c0, c1, c2, and c3 in that order"
      attr(err,"error") <- -1
      } else {if(!identical(t[1:14],c("id","evid","time","dur","dose","addl","ii","input","out","outeq","c0","c1","c2","c3"))){
      err$colorder$msg <- "FAIL - The first 14 columns must be named id, evid, time, dur, dose, addl, ii, input, out, outeq, c0, c1, c2, and c3 in that order."
      attr(err,"error") <- -1}}
   
    # 2 - check to make sure ids and cols are 11 char or less
    t <- which(nchar(as.character(data2$id))>11) | nchar(names(data2)>11)
    if(length(t)>0) {
      err$maxchar$msg <- "FAIL - The following row numbers have columns that contain entries >11 characters:"
      err$maxchar$results <- t
      attr(err,"error") <- -1
    } 
    # 3 - check for observations at time 0
    t <- which(data2$evid==0 & data2$time==0)
    if(length(t)>0) {
      err$obsT0$msg <-"FAIL - The following row numbers have observations at time 0:"
      err$obsT0$results <- t
      attr(err,"error") <- -1
    } 
    # 4 - check for NA observations (should be -99)
    t <- which(is.na(data2$out) & data2$evid==0)
    if(length(t)>0) {
      err$obsMiss$msg <- "FAIL - The following row numbers are observations (evid=0) with no entry (should be -99 if missing):"
      err$obsMiss$results <- t
      attr(err,"error") <- -1
    } 
    # 5 - check for complete dose records
    t <- which(data2$evid!=0 & (is.na(data2$time) | is.na(data2$dur) | is.na(data2$dose) | is.na(data2$input)))
    if(length(t)>0) {
      err$doseComp$msg <- "FAIL - The following row numbers have incomplete dose records (unused addl or ii should have '.' placeholders):"
      err$doseComp$results <- t
      attr(err,"error") <- -1
    } 
    # 6 - check for complete observation records
    t <- which(data2$evid==0 & (is.na(data2$time) | is.na(data2$out) | is.na(data2$outeq)))
    if(length(t)>0) {
      err$obsComp$msg <- "FAIL - The following row numbers have incomplete observation records:"
      err$obsComp$results <- t
      attr(err,"error") <- -1
    } 
    # 7 - check for evid=1 as first record for each subject
    t <- which(tapply(data2$evid,data2$id,function(x) x[1])!=1)
    t2 <- match(names(t),data2$id)
    if(length(t)>0){
      err$evid1$msg <- "FAIL - The following row numbers do not have a dose (evid=1) as the first record:"
      err$evid1$results <- t2
      attr(err,"error") <- -1
    } 
    # 8 - check for time=0 for each subject as first record
    t <- which(tapply(data2$time,data2$id,function(x) x[1])!=0)
    t2 <- match(names(t),data2$id)
    if(length(t)>0){
      err$T0$msg <- "FAIL - The following row numbers do not have time=0 as first record:"
      err$T0$results <- t2
      attr(err,"error") <- -1
    } 
    #9 was deleted (former requirement to have last event be observation)
    # 10,11 - covariate checks
    numcol <- ncol(data2)
    if(numcol>14){
      # 10 - check for missing covariates at time 0
      time0 <- which(data2$time==0 & data2$evid==1)
      if(length(time0)>1) {t <- apply(as.matrix(data2[time0,15:numcol],ncol=numcol-14),1,function(x) any(is.na(x)))} else {t <- is.na(time0)}
      if(length(time0[t])>0){
        err$covT0$msg <- "FAIL - The following row numbers have missing covariate data at time 0."
        err$covT0$results <- time0[t]
        attr(err,"error") <- -1
      } else {err$covT0$msg <- "OK - all subjects have covariate data at time 0."}
      
      # 11 - check for non-numeric covariates
      t <- which(sapply(data2[,15:numcol],function(x) !is.numeric(x)))
      if(length(t)>0){
        err$covNumeric$msg <- "FAIL - The following covariates are non-numeric."
        err$covNumeric$results <- names(data2)[14+t]
        attr(err,"error") <- -1
      } else {err$covNumeric$msg <- "OK - all covariates are numeric."}
    }
    
    #12 - check that all times within a given ID block are monotonically increasing
    misorder <- NA
    for(i in 2:nrow(data2)){
      if((data2$time[i]-data2$time[i-1]<0) & data2$id[i]==data2$id[i-1] & data2$evid[i]!=4) misorder <- c(misorder,i) 
    }
    if(length(misorder)>1){
      err$timeOrder$msg <- "FAIL - The following rows are from subject IDs with unsorted time entries."
      err$timeOrder$results <- misorder[-1]
      attr(err,"error") <- -1
    }
    #13 - check that all records for a given subject ID are grouped
    temp <- data.frame(row=as.numeric(rownames(data2)),id=data2$id)
    t <- tapply(temp$row,temp$id,function(x) any(diff(x)>1))
    if(any(t)) {t2 <- which(data2$id %in% unique(data2$id)[t])} else {t2 <- NULL}
    if(length(t2)>0){
      err$contigID$msg <- "FAIL - The following rows are from subject IDs that are not contiguous."
      err$contigID$results <- t2
      attr(err,"error") <- -1
    }
    #14 - check that all records have an EVID value
    t <- which(is.na(data2$evid))
    if(length(t)>0) {
      err$missEVID$msg <-"FAIL - The following row numbers have missing EVID values:"
      err$missEVID$results <- t
      attr(err,"error") <- -1
    } 
    #15 - check that all records have an EVID value
    t <- which(is.na(data2$time))
    if(length(t)>0) {
      err$missTIME$msg <-"FAIL - The following row numbers have missing TIME values:"
      err$missTIME$results <- t
      attr(err,"error") <- -1
    } 
    #16 - check that all non-missing colmuns that have to be numeric are numeric
    allMiss <- which(apply(data2[,2:15],2,function(x) all(is.na(x))))
    nonNumeric <- which(sapply(data2[,2:15],function(x) !is.numeric(x)))
    if(length(allMiss)>0){
      nonNumeric <- nonNumeric[!nonNumeric %in% allMiss]
    } 
    if(length(nonNumeric)>0){
      err$nonNum$msg <-"FAIL - The following columns must be all numeric:"
      err$nonNum$results <- toupper(c("evid","time","dur","dose","addl","ii","input","out","outeq","c0","c1","c2","c3")[nonNumeric])
      attr(err,"error") <- -1
    }
    
    
    
    class(err) <- c("PMerr","list")
    if(!quiet) cat("\nDATA FILE REPORT:\n")
    if(!quiet) {print(err);flush.console()}
    
  
    
    #if no errors in data, and model is specified, check it for errors too
    if(all(unlist(sapply(err,function(x) is.na(x$results)))) & !is.na(model)){
      #get information from data      
      dataoffset <- 2*as.numeric("addl" %in% names(data2))
      ncov <- ncol(data2)-(12+dataoffset)
      if(ncov>0) {covnames <- names(data2)[(13+dataoffset):ncol(data2)]} else {covnames <- NA}
      numeqt <- max(data2$outeq,na.rm=T)  
      
      modeltxt <- model
      #attempt to translate model file into separate fortran model file and instruction files
      engine <- list(alg="NP",ncov=ncov,covnames=covnames,numeqt=numeqt,indpts=-99,limits=NA)
      
      if(!quiet) cat("\nMODEL REPORT:\n")
      trans <- makeModel(model=model,data=data,engine=engine,write=F,silent=quiet)
      
      if(trans$status==-1){
        if("mQRZZZ.txt" %in% Sys.glob("*",T)){
          file.remove(model)
          file.rename("mQRZZZ.txt",model)
        }
        if(!quiet) cat("There are errors in your model file.\n")
        cat(trans$msg)
      }
      if(trans$status==0) {
        if(!inherits(data,"PMmatrix") & !quiet) cat(paste("Associated data file: ",data,sep=""))
        if(!quiet) cat("\nYou are using a Fortran model file rather than the new text format.\nThis is permitted, but see www.lapk.org/ModelTemp.php for details.\n")
      }
      if(trans$status==1){
        if(!inherits(data,"PMmatrix") & !quiet) cat(paste("Associated data file: ",data,sep=""))
        if(!quiet) cat("\nExcellent - there were no errors found in your model file.\n")
        if(trans$model %in% Sys.glob("*",T)){
          file.remove(trans$model)
        }
      }
    }
    if(!quiet) flush.console()
    return(err)
  }

#here's the function to try and fix errors in the data file
errfix <- function(data2,model,quiet){
  report <- NA
  err <- errcheck(data=data2,model=model,quiet=quiet)
  # 1 - fix first 14 columns
  if(length(grep("FAIL",err$colorder$msg))>0){
    fixedcols <- c("id","evid","time","dur","dose","addl","ii","input","out","outeq","c0","c1","c2","c3")
    t <- tolower(names(data2))
    PMcols <- match(fixedcols,t)
    if(any(is.na(PMcols))) {
      misscols <- fixedcols[is.na(PMcols)]
      report <- c(report,paste("Cannot fix columns; the following are missing: ",paste(misscols,collapse="'', '"),".",sep=""))
    } else {
      covcols <- (1:ncol(data2))[!(1:ncol(data2)) %in% PMcols]
      data2 <- data2[,c(PMcols,covcols)]
      report <- c(report,paste("Columns are now ordered appropriately."))  
    }
  }
  # 2 - Make sure ids and cols are 11 char or less
  if(length(grep("FAIL",err$maxchar$msg))>0){
    names(data2) <- substr(names(data2),1,11)
    report <- c(report,paste("Column names are all 11 characters or fewer."))    
  }
  # 3 - remove observations at time 0
  if(length(grep("FAIL",err$obsT0$msg))>0){
    data2 <- data2[-err$obsT0$results,]
    report <- c(report,paste("Observations at time 0 have been removed."))    
    err <- errcheck(data=data2,model=model,quiet=quiet)
  }
  # 4 - check for NA observations (should be -99)
  if(length(grep("FAIL",err$obsMiss$msg))>0){
    data2 <- data2[err$obsMiss$results,"out"] < -99
    report <- c(report,paste("Missing observations for evid=0 have been replaced with -99."))    
    err <- errcheck(data=data2,model=model,quiet=quiet)
  }
  # 5 - check for complete dose records
  if(length(grep("FAIL",err$doseComp$msg))>0){
    report <- c(report,paste("Dose records (evid=1 or evid=4) must have time, duration, dose and input; addl and ii should be '.' if not needed.  Run PMcheck and fix manually."))    
  }
  # 6 - check for complete observation records
  if(length(grep("FAIL",err$obsComp$msg))>0){
    report <- c(report,paste("Observation records (evid=0) must have time, out, and outeq. Run PMcheck and fix manually."))    
  }
  # 7 - flag evid!=1 as first event
  if(length(grep("FAIL",err$evid1$msg))>0){
    report <- c(report,paste("The first event for every subject must be a dose (evid=1). Run PMcheck and fix manually."))    
  }
  # 8 - insert dummy doses of 0 for those missing time=0 first events
  if(length(grep("FAIL",err$T0$msg))>0){
    T0 <- data2[err$T0$results,]
    T0$time <- 0; T0$evid <- 1; T0$dose <- 0
    data2 <- rbind(data2,T0)
    data2 <- data2[order(data2$id,data2$time),]
    report <- c(report,paste("Subjects with first time > 0 have had a dummy dose of 0 inserted at time 0."))    
    err <- errcheck(data=data2,model=model,quiet=quiet)
  }

  #10 - alert for missing covariate data2
  if(length(grep("FAIL",err$covT0$msg))>0){
    report <- c(report,paste("All covariates must have values for each subject's first event.  Run PMcheck and fix manually."))    
  }
  #11 - change non-numeric to numeric covariates
  if(length(grep("FAIL",err$covNumeric$msg))>0){
    covcols <- which(names(data2) %in% err$covNumeric$results)
    data2[,covcols] <- lapply(covcols,function(x) as.numeric(factor(data2[,x])))
    report <- c(report,paste("Non-numeric covariates have been converted to numeric."))    
  }
  
  #12 - reorder times
  if(length(grep("FAIL",err$timeOrder$msg))>0){
    if(any(data2$evid==4)){
      report <- c(report,paste("Your dataset has EVID=4 events. Unable to sort times automatically."))    
    } else {
      data2 <- data2[order(data2$id,data2$time),]
      report <- c(report,paste("Times within each subject have been ordered."))}    
  }
  #13 - reorder IDs
  if(length(grep("FAIL",err$contigID$msg))>0){
    if(any(data2$evid==4)){
      report <- c(report,paste("Your dataset has EVID=4 events. Unable to sort subjecst and times automatically."))    
    } else {
      data2 <- data2[order(data2$id,data2$time),]
      report <- c(report,paste("Subjects have been grouped and ordered."))}    
  }
  #14 - fix missing EVID
  if(length(grep("FAIL",err$missEVID$msg))>0){
    data2$evid[err$missEVID$results] <- ifelse(is.na(data2$dose[err$missEVID$results]),0,1)
    report <- c(report,paste("EVID for events with doses changed to 1, otherwise 0."))    
  }
  
  #15 - report missing TIME
  if(length(grep("FAIL",err$missTIME$msg))>0){
    report <- c(report,paste("Your dataset has missing times.  Please fix manually."))    
  }
  
  #16 - report non-numeric columns
  if(length(grep("FAIL",err$nonNum$msg))>0){
    report <- c(report,paste("Your dataset has non-numeric columns.  Please fix manually."))    
  }
  
  if(!quiet) cat("\nFIX DATA REPORT:\n")
  print(report[-1])
  flush.console()
  row.names(data2) <- 1:nrow(data2)
  return(data2)
}

#get the data
if(is.character(data)) {
  data2 <- tryCatch(suppressWarnings(PMreadMatrix(data,quiet=T)),error = function(e) return(invisible(e)))
} else {data2 <- data}
if(missing(model)) model <- NA

#check for errors
err <- errcheck(data2,model=model,quiet=quiet)
if(length(err)==1){
  cat("You must at least have id, evid, and time columns to proceed with the check.\n")
  flush.console()
  return(invisible(NULL))
}
maxTime <- max(data2$time,na.rm=T)
if(maxTime>24*48 & !quiet) cat(paste("Warning: The maximum number of AUC intervals in NPAG is 48.\nYour longest event horizon is ",maxTime," hours.\nPmetrics will automatically choose an AUC interval of at least ",ceiling(maxTime/48)," hours during an NPAG run.\nYou can calculate AUCs for other intervals after the run using makeAUC().\n\n",sep=""))
flush.console()


#try to fix errors if asked
if(fix){
  if( attr(err,"error")==0){
    cat("\nFIX DATA REPORT:\n\nThere were no errors found in your data file.\n")
    return(invisible(err))
  } else {
    newdata <- errfix(data=data2,model=model,quiet=quiet)
    flush.console()
    return(invisible(newdata))
  }    
} else { #didn't ask to fix errors so return error object
  return(invisible(err))
}

}










