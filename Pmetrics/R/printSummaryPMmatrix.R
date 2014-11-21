#' Print the summary of a Pmetrics PMmatrix object
#'
#' Summarize the raw data used for a Pmetrics run.
#'
#' @title Summarize Covariates and Bayesian Posterior Parameter Values
#' @method print summary.PMmatrix
#' @param x A summary.PMmatrix object made by \code{\link{summary.PMmatrix}} 
#' @return A formatted printing of a \emph{summary.PMmatrix} object
#' @author Michael Neely
#' @seealso \code{\link{summary.PMmatrix}}
#' @export

print.summary.PMmatrix <- function(x){
  # order of objects
  #   nsub
  #   ndrug
  #   numeqt
  #   nobsXouteq
  #   missObsXouteq
  #   ncov
  #   ndoseXid
  #   nobsXid
  #   doseXid 
  #   obsXid 
  #   formula 
  attach(x)
  cat(paste("\nNumber of subjects:",nsub,"\n"))
  cat(paste("Number of inputs:",ndrug,"\n"))
  cat(paste("Number of outputs:",numeqt,"\n"))
  for(i in 1:numeqt){
    cat(paste("Total number of observations (outeq ",i,"): ",nobsXouteq[i],", with ",missObsXouteq[i]," (",sprintf("%.3f",missObsXouteq[i]/nobsXouteq[i]),"%) missing\n",sep=""))
    
  }
  cat(paste("Number of covariates:",ncov,"\n"))
  cat(paste("\nTHE FOLLOWING ARE MEAN (SD), MIN TO MAX\n",paste(rep("-",75),collapse=""),"\n",sep=""))
  cat("\nINPUTS\n")
  for(i in 1:ndrug){
    if(ndrug==1){
      cat(paste("Number of doses per subject (input ",i,"): ",sprintf("%.3f",mean(ndoseXid,na.rm=T))," (",sprintf("%.3f",sd(ndoseXid,na.rm=T)),"), ",sprintf("%.3f",min(ndoseXid,na.rm=T))," to ",sprintf("%.3f",max(ndoseXid,na.rm=T)),"\n",sep=""))
      cat(paste("Dose per subject (input ",i,"): ",sprintf("%.3f",mean(unlist(doseXid),na.rm=T))," (",sprintf("%.3f",sd(unlist(doseXid),na.rm=T)),"), ",sprintf("%.3f",min(unlist(doseXid),na.rm=T))," to ",sprintf("%.3f",max(unlist(doseXid),na.rm=T)),"\n",sep=""))
    } else {
      cat(paste("Number of doses per subject (input ",i,"): ",sprintf("%.3f",mean(ndoseXid[,i],na.rm=T))," (",sprintf("%.3f",sd(ndoseXid[,i],na.rm=T)),"), ",sprintf("%.3f",min(ndoseXid[,i],na.rm=T))," to ",sprintf("%.3f",max(ndoseXid[,i],na.rm=T)),"\n",sep=""))
      cat(paste("Dose (input ",i,"): ",sprintf("%.3f",mean(unlist(doseXid[,i]),na.rm=T))," (",sprintf("%.3f",sd(unlist(doseXid[,i]),na.rm=T)),"), ",sprintf("%.3f",min(unlist(doseXid[,i]),na.rm=T))," to ",sprintf("%.3f",max(unlist(doseXid[,i]),na.rm=T)),"\n",sep=""))
    }
  }
  cat("\nOUTPUTS\n")
  for(i in 1:numeqt){
    if(numeqt==1){
      nobs <- unlist(nobsXid)
      obs <- unlist(obsXid)
    } else {
      nobs <- unlist(nobsXid[,i])
      obs <- unlist(obsXid[,i])
    }
    obs <- obs[obs!=-99]
    cat(paste("Number of obs per subject (outeq ",i,"): ",sprintf("%.3f",mean(nobs,na.rm=T))," (",sprintf("%.3f",sd(nobs,na.rm=T)),"), ",sprintf("%.3f",min(nobs,na.rm=T))," to ",sprintf("%.3f",max(nobs,na.rm=T)),"\n",sep=""))
    cat(paste("Observation per subject (outeq ",i,"): ",sprintf("%.3f",mean(obs,na.rm=T))," (",sprintf("%.3f",sd(obs,na.rm=T)),"), ",sprintf("%.3f",min(obs,na.rm=T))," to ",sprintf("%.3f",max(obs,na.rm=T)),"\n",sep=""))
  }
  if(ncov>0){
    cat("\nCOVARIATES\n")
    for(i in 1:ncov){
      cat(paste(covnames[i],": ",sprintf("%.3f",mean(unlist(cov[[i]]),na.rm=T))," (",sprintf("%.3f",sd(unlist(cov[[i]]),na.rm=T)),"), ",sprintf("%.3f",min(unlist(cov[[i]]),na.rm=T))," to ",sprintf("%.3f",max(unlist(cov[[i]]),na.rm=T)),"\n",sep=""))
    }
  }
  
  if(!is.null(x$formula)){
    cat(paste("\nFormula\n",paste(rep("-",75),collapse=""),"\n",sep=""))
    print(x$formula)
  }
  cat(paste(paste(rep("-",75),collapse=""),"\nNote: See help(summary.PMmatrix) for accessing specific items by name.\n",sep=""))
  
  detach(x)
  
} #end function
