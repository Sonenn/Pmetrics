#' Combines output(s) of makePTA to create new combination targets
#'
#' 
#'
#' @title Combine \emph{PMpta} objects
#' @param pta1 A \emph{PMpta} object created by \code{\link{makePTA}}.
#' @param pta2 An optional second \emph{PMpta} object to be combined. If \code{pta2} is provided, the number of simulations and profiles must match. If \code{success} is numerical, \code{pta2} will be ignored.
#' @param targets1 A single value or vector of values matching target(s) from \code{pta1}. If missing, all targets will be used. 
#' @param targets2 A single value or vector of values specifying target(s) to be combined with \code{targets1}. If a vector is 
#' supplied, its length must match length of \code{targets1}. If missing, all targets of \code{pta2} (if specified) will be used. If \code{success} is numerical, \code{targets2} will be ignored. 
#' @param success Logical operator to be applied or new success threshold for comparison. If a new success threshold is specified, a new \emph{PMpta} object is created using \link{\code{makePTA}}. 
#' Numerical value means new success threshold. Operators can be "AND", "OR", "NOT", "XOR", or "DIFF". \code{XOR} performs exclusive OR, returning \code{True} if one or the other but not both are true. 
#' Each operator can be preceded or followed by \dQuote{!} to invert the success value of the first or second target. See details and examples for more on how to use these to combine different types of \emph{PMpta} objects.
#' @param silent Suppress messages
#' @return The output of \code{makePTA} is a list of class \emph{PMpta}. See \code{\link{makePTA}} for details.
#' 
#' @author Jan Strojil & Michael Neely
#' @seealso \code{\link{makePTA}}


combinePTA <- function(pta1, pta2, targets1, targets2, success){
  if(missing(pta1)) stop("You need to specify at least one PMpta object.", call. = F)
  if(!inherits(pta1,"PMpta")) {
    stop("Object pta1 is not an PMpta object.", call. = F)
  }
  if(!missing(pta2)) if(!inherits(pta2,"PMpta")) {
    stop("Object pta2 is not an PMpta object.", call. = F)
  }
  #browser()
  pta1_targets <- unique(pta1$outcome$target)
  if (!all(targets1 %in% pta1_targets)) {
    cat("\nThe following targets were not found in pta1: ", paste(targets1[!(targets1 %in% pta1_targets)], collapse = ", "),sep = "")
    cat("\nThis is the list of targets in pta1: ", paste(pta1_targets, collapse = ", "))
    stop("Targets1 contains values not found in object pta1.", call. = F)
  }

  if (missing(targets2) & !missing(pta2)){
    targets2 <- unique(pta2$outcome$target) 
  } 
  else{
    if (!missing(pta2)) {
      pta2_targets <- unique(pta2$outcome$target)
      if (!all(targets2 %in% pta2_targets)) {
        cat("\nThe following targets were not found in pta1: ", paste(targets2[!(targets2 %in% pta2_targets)], collapse = ", "),sep = "")
        cat("\nThis is the list of targets in pta2: ", paste(pta2_targets, collapse = ", "))
        stop("'Targets2' contains values not found in object pta2.", call. = F)
      }
    }
    else if (!all(targets2 %in% pta1_targets)) {
      cat("\nThe following targets were not found in pta1: ", paste(targets2[!(targets2 %in% pta1_targets)], collapse = ", "),sep = "")
      cat("\nThis is the list of targets in pta1: ", paste(pta1_targets, collapse = ", "))
      stop("'Targets2' contains values not found in object pta1.", call. = F)
    }
  }
  
  if (length(targets2)>1 & length(targets2)!=length(targets1)) {
    cat("\n", length(targets2), " values were specified for 'targets2' but there are ", length(targets1), " values in 'targets1'.", sep = "")
    cat("\nIf 'targets2' contains more than 1 value, the lengths of 'targets1' and 'targets2' must match.")
    stop("Non-matching numbers of targets.", call. = F)
  }
  
  if (!missing(pta2)) sametype <-  (attr(pta1, "type")==attr(pta2, "type")) else sametype <- T

  if (length(targets2)==1) targets2 <- rep(targets2[[1]], length(targets1))
  
  if (sametype & (any(targets1==targets2))) {
    cat("Targets (", paste(targets1[targets1==targets2], collapse = ", "), ") cannot be combined with themselves.", sep = "")
    stop("'targets2' conflicts with 'targets1'.", call. = F)
  }
      
  simnum1 <- max(pta1$results$simnum)
  if (!missing(pta2)) {
    simnum2 <- max(pta2$results$simnum)
    if (simnum2 != simnum1) {
      cat("PMpta objects contain different numbers of simulations, cannot combine.", sep = "")
      stop("Numbers of simulations must match.", call. = F)
    }
  }
  
  browser() 
  if(missing(pta2))
  
  x$results <- x$results[x$results$target %in% c(targ1,targ2),]   
  diffPTA <- -1*tapply(x$results$pdi,list(x$results$simnum,x$results$id),diff)
  nsim <- dim(diffPTA)[1]
  nsub <- dim(diffPTA)[2]
  targetlab <- paste("(",targ1," to ",targ2,")",sep="")
  resultDF <- data.frame(simnum=rep(1:nsim,each=nsub),id=rep(1:nsub,nsim),target=rep(targetlab,nsub*nsim),pdi=c(t(diffPTA)))
  newsuccSimXtarg <- tapply(resultDF$pdi, list(resultDF$target, resultDF$simnum), function(x) sum(x >= success)/sum(!is.na(x)))
  newmeanpdi <- tapply(resultDF$pdi, list(resultDF$target, resultDF$simnum), mean, na.rm = T)
  newsdpdi <- tapply(resultDF$pdi, list(resultDF$target, resultDF$simnum), 
                     sd, na.rm = T)
  
  nsim <- max(resultDF$simnum)
  ntarg <- length(unique(resultDF$target))
  targets <- c(unique(resultDF$target))
  newpta.outcome <- data.frame(simnum = 1:nsim, 
                               target = rep(targetlab, nsim), prop.success = c(newsuccSimXtarg), 
                               pdi.mean = c(newmeanpdi), pdi.sd = c(newsdpdi))
  
  rval <- list(results = resultDF, outcome = newpta.outcome)
  attr(rval, "simlabels") <- attr(x, "simlabels")
  attr(rval, "simTarg") <- attr(x, "simTarg")
  class(rval) <- c("PMpta", "list")
  return(rval)
  
}
