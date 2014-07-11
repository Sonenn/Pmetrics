#' Calculates the Percent Target Attainment (PTA)
#'
#' \code{makePTA} will calculate the PTA for any number of simulations, targets and definitions of success.
#' Simulations typically differ by dose, but may differ by other features such as children vs. adults.
#'
#' @title Calculation of PTAs
#' @param simdata A vector of simulator output filenames, e.g. c(\dQuote{simout1.txt},\dQuote{simout2.txt}),
#' with wildcard support, e.g. \dQuote{simout*} or \dQuote{simout?}, or
#' a list of PMsim objects made by \code{\link{SIMparse}} with suitable simulated doses and observations.  The number and times of simulated
#' observations does not have to be the same in all objects.
#' @param targets A vector of pharmacodynamic targets, such as Minimum Inhibitory Concentrations (MICs), e.g. c(0.25, 0.5,1,2,4,8,16,32)
#' @param target.type A numeric or character vector, length 1.  If numeric, must correspond to an observation time common to all PMsim objects in
#' \code{simdata}, rounded to the nearest hour.  In this case, the target statistic will be the ratio of observation at time \code{target.type} to target.  This enables 
#' testing of a specific timed concentration (e.g. one hour after a dose or C1) which may be called a peak, but is not actually the maximum drug
#' concentration.  Be sure that the time in the simulated data is used, e.g. 122 after a dose given at 120.  Character values may be one of 
#' \dQuote{time}, \dQuote{auc}, \dQuote{peak}, or \dQuote{min}, for, respectively, percent time above target within the time range
#' specified by \code{start} and \code{end}, ratio of area under the curve within the time range to target, ratio of peak concentration within the time range 
#' to target, or ratio of minimum concentration within the time range to target.  
#' @param success A single value specifying the success statistic, e.g. 0.4 for proportion time (end-start) above target, or 100 for peak:target.
#' @param outeq An integer specifying the number of the simulated output equation to use. Default is 1.
#' @param free.fraction Proportion of free, active drug.  Default is 1, i.e. 100\% free drug or 0\% protein binding.
#' @param start Specify the time to begin PTA calculations. Default is a vector with the first observation time for subjects
#' in each element of \code{simdata}, e.g. dose regimen. If specified as a vector, values will be recycled as necessary.
#' @param end Specify the time to end PTA calculations so that PTA is calculated
#' from \code{start} to \code{end}.  Default for end is the maximum observation
#' time for subjects in each element of \code{simdata}, e.g. dose regimen.  If specified as a vector, values will be recycled
#' as necessary. Subjects with insufficient data (fewer than 5 simulated observations) for a specified interval will trigger a warning.
#' Ideally then, the simulated datset should contain sufficient observations within the interval specified by \code{start} and \code{end}.
#' @return The output of \code{makePTA} is a list of class \emph{PMpta},
#' which has 2 objects:
#' \item{results }{A data frame with the following columns: simnum,id,target,ratio.  
#' \emph{simnum} is the number of the simulation; \emph{id} is the simulated profile number
#' within each simulation; \emph{target} is the specified target; and \emph{ratio} is
#' the target ratio, e.g. time > target, auc:target, etc.}
#' \item{outcome }{A data frame summarizing the results with the following columns: simnum, target, success, meanratio, and sdratio.
#' \emph{simnum} and \emph{target} are as for \code{outcome}.  The \emph{prop.success} column has the proportion with a ratio > \code{success},
#' as specified in the function call.  The \emph{mean.stat} and \emph{sd.stat} columns have the 
#' mean and standard deviation of the target statistic (e.g. proportion end-start above target, ratio of Cmax to target) for each simulation and target.}  
#' @author Michael Neely
#' @seealso \code{\link{plot.PMpta}}, \code{\link{SIMparse}}

makePTA <- function(simdata,targets,target.type,success,outeq=1,free.fraction=1,start,end){
  if(missing(simdata) | missing(target.type)) stop("Simulation output and target.type must be specified.\n")
  if(is.character(target.type) & !target.type %in% c("time","auc","peak","min")) stop("Please specify target.type as a numerical value corresponding to a common\ntime in all simulated datasets, or a character value of 'time', 'auc', 'peak' or 'min'.\n")
  if(!inherits(simdata,"list")){ #so we are dealing with names of files
    simfiles <- Sys.glob(simdata)
    if(length(simfiles)==0) stop("There are no files matching \"",simdata,"\".\n",sep="")
    simdata <- list()
    for(i in 1:length(simfiles)){
      simdata[[i]] <- tryCatch(SIMparse(simfiles[i]),
                               error=function(e) stop(paste(simfiles[i],"is not a PMsim object.\n")))
    }
  }
  #check for one PMsim object only, and if so, make it a one-item list
  if(!inherits(simdata[[1]],"PMsim")) {simdata <- list(simdata);class(simdata) <- c("PMsim","list")}
  #number of sims
  nsim <- length(simdata)
  #replicate start and end times if supplied for each simulation
  if(!missing(start)) {start <- rep(start,nsim)}
  if(!missing(end)) {end <- rep(end,nsim)}
  #number of targets
  ntarg <- length(targets)
  #the list to hold the PTA results
  results <- list()
  
  cat("\nCalculating PTA for each simulation and target...\n")
  flush.console()
  if(target.type=="time"){
    maxpb <- nsim*ntarg
  } else {maxpb <- nsim}
  pb <- txtProgressBar(min = 0, max = maxpb, style = 3)
  
  #loop through each simulation, calculating PTA
  for(simnum in 1:nsim){
    
    #get the simulated data for sim
    wrk.sim <- simdata[[simnum]]$obs
    #get the correct outeq
    wrk.sim <- wrk.sim[wrk.sim$outeq==outeq,]
    #take out missing observations
    wrk.sim <- wrk.sim[!is.na(wrk.sim$out),]
    #multiply by free fraction
    wrk.sim$out <- wrk.sim$out*free.fraction
    #simulated times
    wrk.times <- unique(wrk.sim$time)
    #if start and end times missing, set them to min/max, else use those supplied
    if(missing(start)) {wrk.start <- min(wrk.times)} else {wrk.start <- start[simnum]}
    if(missing(end)) {wrk.end <- max(wrk.times)} else {wrk.end <- end[simnum]} 
    if(wrk.start>=wrk.end) {stop(paste("For simulation ",simnum,", start is not less than end/n",sep=""))}
    #filter simulated data by start/end times
    wrk.sim <- wrk.sim[wrk.sim$time>=wrk.start & wrk.sim$time<=wrk.end,]
    if(length(wrk.sim)==0){
      cat(paste("Note: Simulation ",simnum," omitted because no simulated observations fell within the time window defined by start and end.\n",sep=""))
      next
    } 
    #recheck times after filtering
    wrk.times <- unique(wrk.sim$time)
    #number of observations
    wrk.nobs <- length(wrk.times)
    if(wrk.nobs<5) warning(paste("Only ",wrk.nobs," simulated observations available for simulation ",simnum,".\nThis can compromise estimates of target attainment.\nStrongly consider increasing the number of simulated observations.\n",sep=""))
    
    #time above target
    if(target.type=="time"){
      #function to calculate time above target for pair of times/outs
      timeabove <- function(times,outs,targ){
        #both outs are below target
        if(outs[1]<targ & outs[2]<targ) interval <- 0
        #both outs are at or above target
        if(outs[1]>=targ & outs[2]>=targ) interval <- times[2]-times[1]
        #first is below, second is at or above
        if(outs[1]<targ & outs[2]>=targ){
          lm.1 <- lm(times~outs)
          cross1 <- predict(lm.1,data.frame(outs=targ))
          interval <- times[2]-cross1
        }
        #first is at or above, second is below
        if(outs[1]>=targ & outs[2]<targ){
          lm.1 <- lm(times~outs)
          cross1 <- predict(lm.1,data.frame(outs=targ))
          interval <- cross1-times[1]
        }
        return(interval) #the time above target
      }
      
      #function to split data into blocks of 2 rows
      pairUp <- function(sim){
        outs <- lapply(1:(nrow(sim)-1),function(x) c(sim$out[x],sim$out[x+1]))
        times <- lapply(1:(nrow(sim)-1),function(x) c(sim$time[x],sim$time[x+1]))
        return(list(times,outs))
      }
      
      #function to calculate cumulative time above target
      cumTime <- function(sim,targ){
        pairs <- pairUp(sim)
        npairs <- length(pairs[[1]])
        interval <- sum(unlist(lapply(1:npairs,function(x) timeabove(times=pairs[[1]][[x]],outs=pairs[[2]][[x]],targ=targ))))
        #divide total time in the interval by the end-start interval
        return(interval/(wrk.end-wrk.start))
      }
      
      #get the results, which is initially a list [[ntarg]][nsim]
      pta <- list()
      for(t in 1:ntarg){
        targ <- targets[t]
        pta[[t]] <- by(wrk.sim,wrk.sim$id,function(x) cumTime(x,targ=targ))
        setTxtProgressBar(pb, (simnum-1)*ntarg + t)    
        
      }
      #get results into a format consistent with the others, i.e. matrix [ntarg,nsim]
      results[[simnum]] <- do.call(rbind,pta)
      if(ntarg==1) results[[simnum]] <- matrix(results[[simnum]],nrow=1)
    }
    #auc above target
    if(target.type=="auc"){      
      auc <- by(wrk.sim,wrk.sim$id,function(x) makeAUC(x,out~time)[,2])
      results[[simnum]] <- sapply(auc,function(x) x/targets) #matrix [ntarg,nsim]
      if(ntarg==1) results[[simnum]] <- matrix(results[[simnum]],nrow=1)
      setTxtProgressBar(pb, simnum)    
    }
    #peak above target
    if(target.type=="peak"){      
      peak <- tapply(wrk.sim$out,wrk.sim$id,max)
      results[[simnum]] <- sapply(peak,function(x) x/targets) #matrix [ntarg,nsim]
      if(ntarg==1) results[[simnum]] <- matrix(results[[simnum]],nrow=1)
      setTxtProgressBar(pb, simnum)    
    }
    #min above target
    if(target.type=="min"){      
      minobs <- tapply(wrk.sim$out,wrk.sim$id,min)
      results[[simnum]] <- sapply(minobs,function(x) x/targets) #matrix [ntarg,nsim]
      if(ntarg==1) results[[simnum]] <- matrix(results[[simnum]],nrow=1)
      setTxtProgressBar(pb, simnum)    
    }
    #specific obs above target
    if(is.numeric(target.type)){  #specific timed sample    
      timed <- by(wrk.sim,wrk.sim$id,function(x) x$out[round(x$time,2)==target.type])
      results[[simnum]] <- sapply(timed,function(x) x/targets) #matrix [ntarg,nsim]
      if(ntarg==1) results[[simnum]] <- matrix(results[[simnum]],nrow=1)
      setTxtProgressBar(pb, simnum)    
    }
    
  } #close simnum for loop  
  close(pb)
  require(reshape2,warn.conflicts=F,quietly=T)
  resultDF <- melt(results)
  names(resultDF) <- c("target","id","ratio","simnum")
  resultDF$target <- targets[resultDF$target]
  resultDF <- resultDF[,c("simnum","id","target","ratio")]
  
  succSimXtarg <- tapply(resultDF$ratio,list(resultDF$target,resultDF$simnum),
                         function(x) sum(x>=success)/sum(!is.na(x)))
  
  meanratio <- tapply(resultDF$ratio,list(resultDF$target,resultDF$simnum),mean,na.rm=T)
  sdratio <- tapply(resultDF$ratio,list(resultDF$target,resultDF$simnum),sd,na.rm=T)
  
  pta.outcome <- data.frame(simnum=rep(1:nsim,each=ntarg),
                            target=rep(targets,nsim),
                            prop.success=c(succSimXtarg),
                            mean.stat=c(meanratio),
                            sd.stat=c(sdratio))

  rval <- list(results=resultDF,outcome=pta.outcome)
  class(rval) <- c("PMpta","list")
  return(rval)
}

