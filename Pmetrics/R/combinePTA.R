#' Combines output(s) of makePTA to create new combination targets
#'
#' \code{combinePTA} will take existing PMpta object(s) and combine them to form a new, more complex PTA object. It can be used to compare PTAs for different targets and/or success thresholds or to combine PTAs of different type. Any target can be also inverted (e.g. time spent below a threshold).
#'
#' @title Combine \emph{PMpta} objects
#' @param pta1 A \emph{PMpta} object created by \code{\link{makePTA}}.
#' @param pta2 An optional second \emph{PMpta} object to be combined. Please note that the two PMpta objects \emph{must} come from the same simulation for the results to have any meaning. A new value for \code{success} cannot be used if \code{pta2} is specified.
#' @param targets1 A single value or vector of values matching target(s) from \code{pta1}. If missing, all targets will be used. 
#' @param targets2 A single value or vector of values specifying target(s) to be combined with \code{targets1}. If a vector is 
#' supplied, its length must match length of \code{targets1}. If missing, all targets of \code{pta2} (if specified) will be used. If \code{success} is numerical, \code{targets2} will be ignored. 
#' @param combine A string specifying the operation to apply. perators can be "AND", "OR", "XOR", "NAND", or "DIFF". \code{XOR} performs exclusive OR, returning \code{True} if one or the other but not both are true. 
#' Each operator can be preceded or followed by \dQuote{!} to invert the success value of the first or second target. See details and examples for more on how to use these to combine different types of \emph{PMpta} objects.
#' @param success A numerical value between 0 and 1. When used with \code{combine="XOR"}, \code{success} specifies a new success threshold for comparison with an existing \emph{PMpta}. When used with \code{combine="DIFF"}, it specifies the proportion of time spent between target levels.
#' @param silent Suppress messages
#' @return The output of \code{makePTA} is a list of class \emph{PMpta}. See \code{\link{makePTA}} for details.
#' @examples 
#' # building a more complex PTA for gentamicin using objects genta.peak (with success=10) and genta.min built with makePTA
#' # the following will define targets with peaks 10 times MIC and troughs below 2
#' genta.comb <- combinePTA(genta.peak, genta.min, targets1=c(0.125, 0.25, 0.5, 1, 2), targets2=c(2), combine="AND!")
#' 
#' # if we wanted to limit peaks to less than 20, we could run the function on the new object again
#' # note: since the peak PMpta in this example was defined with success=10, we use targets2=2 to limit peaks to 20 mg/l
#' # targets1 was omitted, using all targets from genta.comb
#' genta.comb <- combinePTA(genta.comb, genta.peak, targets2=2, combine="AND!")
#' 
#' # compare two ranges for digoxin, makePTA object of type "min" with targets 0.5, 0.8, 1.2, and 2
#' # following line will return proportions of patients with troughs between 0.5 and 1.2, and 0.8 and 2
#' dig.compare <- combinePTA(digPTA, targets1=c(0.5,0.8), targets2=c(1.2,2), combine="XOR")
#' 
#' # to calculate time spent in range, use the combine="DIFF" option
#' # success is e.g. 80% of time spent between 2 and 10 mg/l 
#' time_in_range <- combinePTA(PMPta, targets1=2, targets2=10, combine="DIFF", success=0.8)
#' 
#' # to compare target attainment for different success thresholds, lets assume pta.time with success=1
#' # to see how many more patients satisfy 40% above target instead of 100%, we can use:
#' pta.difference <- combinePTA(pta.time, combine="XOR", success=0.4)
#' 
#' # for more examples see the Pmetrics User Manual
#' 
#' @author Jan Strojil & Michael Neely
#' @seealso \code{\link{makePTA}}


#NOTE - success can be either new success to compare using DIFF (e.g. difference in success for 100 vs 40 above MIC)
#OR it can be success of TIME INSIDE INTERVAL!!!

combinePTA <- function(pta1, pta2=NULL, targets1, targets2=NULL, combine="AND", success, silent=F){
  
  myoperator <- function(pdi1, pdi2, combop, first_i, second_i) {
    if (first_i) {
      pdi1 <- !pdi1
    }
    if (second_i) {
      pdi2 <- !pdi2
    }
    if (combop=="xor") {
      rval <- xor(pdi1, pdi2)
    }
    if (combop=="and"){
      rval <- (pdi1 & pdi2)
    }
    if (combop=="or"){
      rval <- (pdi1|pdi2)
    }
    return(as.integer(rval))
  }
  
    library(plyr)
  # first, sanity checks if passed parameters and combinations make sense

    # check that any PTAs are PMpta class
  if(missing(pta1)) stop("You need to specify at least one PMpta object.", call. = F)
  
  if(!inherits(pta1,"PMpta")) {
    stop("Object pta1 is not a PMpta object.", call. = F)
  }
  if(!missing(pta2)) if(!inherits(pta2,"PMpta")) {
    stop("Object pta2 is not a PMpta object.", call. = F)
  }

  if (attr(pta1, "simTarg") | (!missing(pta2) & attr(pta2,"simTarg"))) {
    stop("Sorry, combinePTA does not currently work with simulated targets (coming soon).")
  } 
  
  combine <- tolower(combine)
  first_inverted <- (substr(combine,1,1)=="!") # is the first target to be inverted?
  second_inverted <- (substr(combine,nchar(combine),nchar(combine))=="!") # is the second to be inverted?
  combop <- substr(combine,as.integer(first_inverted)+1,nchar(combine)-as.integer(second_inverted)) # get rid of leading and trailing "!", if present
  if (!(combop %in% c("and", "or", "xor", "diff"))){
    stop(cat("Invalid combine instruction: ", combop, " (must be AND, OR, XOR, or DIFF)", sep = ""), call. = F)
  }
  
  # for same type value-in-range targets, "XOR" can be used instead of "AND!" since if only one target succeeds, it must be the lower one
  # DIFF is *only* used for TIME IN RANGE! (when comparing targets)
  # - the reason DIFF does not make sense for others is that they are linear and 2xtarget is the same as 2xsuccess, 
  # - and difference between targets 3 and 6 must be 2-fold by definition
  # TIME can be also be compared using XOR (when comparing success rates for different cut-offs (e.g. 40% vs 100%) or success rates for different targets (e.g. MIC and 4xMIC))
  # note - DIFF with success rate of 100% if the same as MIN above X and PEAK below Y.
  # workaround to invert just one target: make the PTA object including a target of 0 (always succeeds) and then !AND - not(T) and 1 -> not(T)
  
  # check that all targets in targets1 are found in PTA1
  pta1_targets <- unique(pta1$results$target)
  if (missing(targets1)) {
    targets1 <- pta1_targets
    if (!silent) {cat("\nNo 'targets1' specified, using all targets from 'pta1'.")}
  } else {
    if (!all(targets1 %in% pta1_targets)) {
      cat("\nThe following targets were not found in pta1: ", paste(targets1[!(targets1 %in% pta1_targets)], collapse = ", "),sep = "")
      cat("\nThis is the list of targets in pta1: ", paste(pta1_targets, collapse = ", "))
      stop("'Targets1' contains values not found in 'pta1'.", call. = F)
    }
  }

  # if PTA2 is specified, see if PTA1 and PTA2 are of the same type
  if (!missing(pta2)) {
    sametype <- (attr(pta1, "type")==attr(pta2, "type")) 
  }
  else {
    sametype <- T # if only one PTA, its the same type by default
  }

  if (!missing(success) & (!missing(pta2))){
    cat("\nIvalid set of parameters. A new success threshold can only be used with a single PMpta object.", sep = "")
    cat("\n- Use combine=\"DIFF\" and two sets of targets for time-in-range calculations.", sep = "")
    cat("\n- Use one set of targets and combine=\"XOR\" for comparing success rates with different success thresholds (e.g. 40% and 100% of time above MIC).", sep = "")
    stop("Too many PMpta objects. Function aborted.", call. = F)
  }

  if (combop=='diff') { # this means TIME-IN-RANGE
    if (attr(pta1, "simTarg")) {
      stop("Cannot use time-between for PMpta objects with simulated targets.", call. = F)
    }
    if (!attr(pta1, "type")=="time"){
      cat(paste("\nA PMpta object of type 'time' is required for time-in-range calculation.", sep = ""))
      cat(paste("\nSupplied PMpta object is of type '",attr(pta1, "type"), "'.", sep = ""))
      stop("Mismatched PMpta type for time-in-range calculation.", call. = F)
    }
    if (missing(targets2)){
      cat(paste("\nNo 'targets2' was specified.", sep = ""))
      cat(paste("\nCalculation requires the upper limits of range.", sep = ""))
      stop("Targets2 not specified.", call. = F)
    }
    if(missing(success)){
      cat("\nTime in range calculation needs a success threshold to determine PTA.", sep = "")
      cat("\nThe PMpta object contains success threshold of ", attr(pta1, "success"), sep="")
      ans <- readline(cat("\nDo you want to use this as the success threshold for time in range?\n1) yes\n2) no (abort)", sep=""))
      if (ans==1) {
        success=attr(pta1, "success")
      }else {
        stop("Function aborted.", call. = F)
      }
    }
  } # end of TIME-IN-RANGE checks
  
  # block of checks for NEW SUCCESS comparison
  if (!missing(success) & combop!="diff"){ 
    if(!missing(targets2) & !all(targets1!=targets2)){
      cat("\nTargets2 specified for new success threshold comparison.", sep = "")
      cat("\nOnly one set of targets is needed, ignoring.", sep="")
      targets2 <- targets1
      warning("Targets2 was specified but ignored.", call. = F)
    }
    if (!(combine %in% c("xor", "!and", "and!"))){
      cat("\nNew success threshold comparison is performed using option combine=\"XOR\"")
      cat("\nThis returns the proportion of patients who succeeded at one but not the other threshold.")
      cat("\nIf you use \"AND\", you will get results equal to that of the lower threshold; if you use \"OR\", you will get results equal to that of the higher threshold.")
      stop("Invalid combine operator.", call. = F)
    }
  }

  if (!missing(targets2)) { # targets2 specified
    if (!missing(pta2)) { # pta2 specified
      pta2_targets <- unique(pta2$results$target)
      if (!all(targets2 %in% pta2_targets)) {
        cat("\nThe following targets were not found in pta2: ", paste(targets2[!(targets2 %in% pta2_targets)], collapse = ", "),sep = "")
        cat("\nThis is the list of targets in pta2: ", paste(pta2_targets, collapse = ", "))
        stop("'Targets2' contains values not found in 'pta2'.", call. = F)
      }
    }
    else {   # PTA2 is missing but targets2 are specified
      if (!all(targets2 %in% pta1_targets)) { #  look for targets2 in PTA1
        cat("\nThe following targets were not found in pta1: ", paste(targets2[!(targets2 %in% pta1_targets)], collapse = ", "),sep = "")
        cat("\nThis is the list of targets in pta1: ", paste(pta1_targets, collapse = ", "))
        stop("'Targets2' contains values not found in 'pta1'.", call. = F)
      } # end of not found in PTA1
    } # end of T2 present but P2 not
  } # end of T2 not missing
  else { #T2 missing
    if (!missing(pta2)) {
      targets2 <- unique(pta2$outcome$target) # if T2 missing but P2 present, get all
      if (!silent) {cat("\nNo 'targets2' specified, using all targets from 'pta2'.")}
    }
    else {
      if (missing(success)) {
        cat("\nNo 'targets2' specified.", sep = "")
        cat("\nOne range of targets with one PMpta object is only possible with new success threshold.", sep="")
        stop("Nothing to combine. 'Targets2' not specified.", call. = F)
      }
      else {
        targets2 <- targets1 # pta2 specified, no targets2, and no success threshold
      }
    }
  }
    
 # if there is more than 1 value in targets2, the number of targets1 and 2 must be the same for 1:1 matching
  if (length(targets2)>1 & length(targets2)!=length(targets1)) {
    cat("\n", length(targets2), " values were specified for 'targets2' but there are ", length(targets1), " values in 'targets1'.", sep = "")
    cat("\nIf 'targets2' contains more than 1 value, the lengths of 'targets1' and 'targets2' must match.")
    stop("Mismatched number of targets.", call. = F)
  }
  
  # if only 1 target in targets2, expand it to the same number as targets1 for 1:1 matching
  if (length(targets2)==1) targets2 <- rep(targets2[[1]], length(targets1))
  
  # if sametype, no new success threshold and a target matches with itself, abort
  if (sametype & missing(success) & (any(targets1==targets2))) {
    cat("The following targets cannot be combined with themselves: ", paste(targets1[targets1==targets2], collapse = ", "), sep = "")
    stop("Aborted, 'targets2' conflicts with 'targets1'.", call. = F)
  }
  
  nsim1 <- max(pta1$results$simnum)
  nprofiles1 <- ddply(pta1$results,.(simnum),numcolwise(max))$id
  if (!missing(pta2)) {
    # get number of simulations and profiles and see if they are the same, abort if not
    nsim2 <- max(pta2$results$simnum)
    nprofiles2 <- ddply(pta2$results,.(simnum),numcolwise(max))$id
    if (nsim1 != nsim2) {
      cat("PMpta objects contain different numbers of simulations, cannot combine.", sep = "")
      cat("\nPlease note: To give meaningful results the PTA objects must come from the same simulation.", call. = F)
      stop("Numbers of simulations must match.", call. = F)
    }
    if (any(nprofiles1 != nprofiles2)) {
      cat("\nPMpta objects contain different numbers of profiles, cannot combine.", sep = "")
      cat("\nPlease note: To give meaningful results the PTA objects must come from the same simulation.", call. = F)
      stop("PMpta objects must come from the same simulation.", call. = F)
    }
    if (sametype) {
      if (any(pta1$results$pdi != pta2$results$pdi)) {
      cat("\nPMpta objects are of same type but contain different PDI results.", sep = "")
      stop("To give meaningful results the PTA objects must come from the same simulation.", call. = F)
      # if ptas are of the same type, we can check if they are the same (which they should be)
      # there is really no reason to have 2 PMpta objects of the same type since the targets
      # could and should be contained in a single object in this case
    }}
  } # end checks that pta2 is the same (if present)
  
  
  # for two TIME objects, there are two options:
  # do a DIFF, disregarding the SUCCESS threshold of the originals and using a new one for time BETWEEN
  #   - e.g. success <- 80% time between 0.5 and 1 mg/l
  # do an AND (OR, XOR...) using the original SUCCESS/FAILURE values
  #   - e.g. successs <-  at least 20% time below 0.5 mg/l and at least 40% time above 5 mg/l
  # for HELP file: specifying "!" for TIME pta means 100-X% below Y
  #   - e.g. if original was made using 60% time above 2 mg/l, "!" would mean AT LEAST 40% of time below 2 mg/l
  

  nsim <- max(pta1$results$simnum)              # get number of simulations
  ntarg <- length(targets1)  # get number of targets pairs
  nsub <- max(nprofiles1)  
  

  ############ TIME IN RANGE ############
  # "diff" means it MUST be time in range with just one pta
  if (combop=="diff"){
    
    # sort pairs of targets to make sure smaller one is first to avoid negative results and get nicer labels
    t2_temp <- pmax(targets1, targets2)
    targets1 <- pmin(targets1, targets2)
    targets2 <- t2_temp
    
    diffPDI <- list()
    for (x in 1:ntarg){
      diffPDI[x] <- list(pta1$results[pta1$results$target %in% targets1[x],]$pdi - pta1$results[pta1$results$target %in% targets2[x],]$pdi)
    }
    
    targetlabels <- paste("(",targets1," to ",targets2,")",sep="")
    resultsDF <- data.frame(simnum=rep(1:nsim,each=nsub),id=rep(1:nsub,nsim),target=rep(targetlabels,each=nsub*nsim),pdi=c(unlist(diffPDI)))
  } else {
    
    # ALL OTHER combinations besides TIME IN RANGE
    
    ########### NEW SUCCESS THRESHOLD ##############
    if (!missing(success)) { # this means a new success rate but not DIFF
      success1 <- attr(pta1, "success")
      success2 <- success
      resultsDF <- pta1$results[pta1$results$target %in% targets1,]  # get results from the PMpta object
      
      
      # the new success must be done before the original PDI is "destroyed" by converting it to T/F
      resultsDF <- cbind(resultsDF, pdi2=unlist(lapply(resultsDF$pdi, function(x) (x >= success2))))
      resultsDF$pdi <- unlist(lapply(resultsDF$pdi, function(x) (x >= success1)))
      resultsDF$comb <- myoperator(resultsDF$pdi, resultsDF$pdi2, combop, first_inverted, second_inverted)
      
      resultsDF$pdi <- NULL
      resultsDF$pdi2 <- NULL
      names(resultsDF)[names(resultsDF)=="comb"] <- "pdi"
      
      ### RESULTSDF created for XOR with new SUCCESS RATE"
    } else {
      
      ########## ALL REMAINING CASES #############

      success1 <- attr(pta1, "success")
      if (!missing(pta2)) {
        success2 <- attr(pta2, "success")
      } else {
        success2 <- success1
      }
      
      #let's built the dataframe by target pairs
      # NOTE: this is fairly slow and could be most likely improved by 
      # a) preallocating the DF
      # b) getting rid of the cycle (by? each? lapply?? melt??? no idea. This works and takes couple seconds.)
      resultsDF <- data.frame()
      for (x in 1:ntarg){
        if (missing(pta2)){
          resultsDF <- rbind(resultsDF, 
                             (cbind(pta1$results[pta1$results$target == targets1[x],],
                                    pdi2=pta1$results$pdi[pta1$results$target == targets2[x]])))
        }
        else{
          resultsDF <- rbind(resultsDF, 
                             (cbind(pta1$results[pta1$results$target == targets1[x],],
                                    pdi2=pta2$results$pdi[pta2$results$target == targets2[x]])))
        }
      }
      names(resultsDF)[names(resultsDF)=="pdi"] <- "pdi1"
      names(resultsDF)[names(resultsDF)=="pta2$results$pdi"] <- "pdi2"
      
      # convert PDIs to success/failure
      resultsDF$pdi1 <- unlist(lapply(resultsDF$pdi1, function(x) (x >= success1)))
      resultsDF$pdi2 <- unlist(lapply(resultsDF$pdi2, function(x) (x >= success2)))
      
      #combine PDIs using the specified operator, applying inversions, if specified
      resultsDF$pdi <- myoperator(resultsDF$pdi1, resultsDF$pdi2, combop, first_inverted, second_inverted)
      
      #delete intermediate PDI1 and PDI2
      resultsDF$pdi1 <- NULL
      resultsDF$pdi2 <- NULL
      
      ### ResultsDF CREATED for all other cases #############
      
    }
  }
  
  # if not doing time-in-range, success is 1 (=TRUE)
  if (combop!="diff"){
    success <- 1
  }
  
  newsuccSimXtarg <- tapply(resultsDF$pdi, list(resultsDF$target, resultsDF$simnum), function(x) sum(x >= success)/sum(!is.na(x)))
  newmeanpdi <- tapply(resultsDF$pdi, list(resultsDF$target, resultsDF$simnum), mean, na.rm = T)
  newsdpdi <- tapply(resultsDF$pdi, list(resultsDF$target, resultsDF$simnum), sd, na.rm = T)
  
  
  # right now the target labels are missing from outcome DF, must figure out a clean way to label them
  newpta.outcome <- data.frame(simnum = 1:nsim, 
                               target = rep(targets1, nsim), prop.success = c(newsuccSimXtarg), 
                               pdi.mean = c(newmeanpdi), pdi.sd = c(newsdpdi))
  newpta.outcome <- newpta.outcome[order(newpta.outcome$target, newpta.outcome$simnum),]

  rval <- list(results = resultsDF, outcome = newpta.outcome)
  attr(rval, "simlabels") <- attr(pta1, "simlabels")
  attr(rval, "simTarg") <- F
  attr(rval, "type") <- "combined"
  class(rval) <- c("PMpta", "list")
  
  
  ##### this block prints out the summary of what was done ######
  type1 <- attr(pta1,"type")
  success1 <- attr(pta1, "success")
  #browser()
  if (!missing(pta2)) {
    type2 <- attr(pta2,"type")
    success2 <- attr(pta2, "success")
  } else { 
    type2 <- attr(pta1,"type")
    if (!missing(success) & combop!="diff"){
      success2 <- success
    }else{
      success2 <- attr(pta1, "success")
    }
  }
  
  if (!silent) {
    cat("\nTargets1: ", paste(targets1, collapse = ", "), sep="")
    cat("\nTargets2: ", paste(targets2, collapse = ", "), sep="")
  }
  
  AB1 <- "above"
  HL1 <- "higher"
  if (first_inverted) {
    AB1 <- "below"
    HL1 <- "lower"
    if (type1 == "time") success1 <- 1-success1
  }
  if (type1=="time") {descr1 <- paste("more than ", success1*100, "% of time ", AB1, " ", sep ="" )} 
  if (type1=="peak") {descr1 <- paste("peaks ",HL1," than ", success1, " times: ", sep ="" )} 
  if (type1=="min") {descr1 <- paste("troughs ",HL1," than ", success1, " times: ", sep ="" )}
  if (is.numeric(type1)) {descr1 <- paste("concentrations at ", type1, " hours ",HL1," than ", success1, " times: ", sep ="" )}
  if (type1=="AUC") {descr1 <- paste("AUCs ",HL1," than ", success1, " times: ", sep ="" )}
  if (type1=="combined") {
    descr1 <- attr(pta1, "success")
    if (first_inverted) paste("NOT (", descr1,")", sep="")
  }
  else {descr1 <- paste(descr1, "[target1] ")}
  if (combop=="diff") {descr1 <- ""}
  

  AB2 <- "above"
  HL2 <- "higher"
  if (second_inverted) {
    AB2 <- "below"
    HL2 <- "lower"
    if (type2 == "time") success2 <- 1-success2
  }
  
  if (type2=="time") {descr2 <- paste("more than ", success2*100, "% of time ",AB2, " ", sep ="" )} 
  if (type2=="peak") {descr2 <- paste("peaks ",HL2," than ", success2, " times: ", sep ="" )} 
  if (type2=="min") {descr2 <- paste("troughs ",HL2," than ", success2, " times: ", sep ="" )}
  if (is.numeric(type2)) {descr2 <- paste("concentrations at ", type2, " hours ",HL2," than ", success2, " times: ", sep ="" )}
  if (type2=="AUC") {descr2 <- paste("AUCs ",HL2," than ", success2, " times: ", sep ="" )}
  if (combop=="diff") {descr2 <- ""}
  if (type2=="combined") {
    descr2 <- paste(attr(pta2, "success"), " ", sep="")
    if (second_inverted) paste("NOT (", descr2,")", sep="")
  }
  else {descr2 <- paste(descr2, "[target2]")}
  
  
  
  if (combine!="diff"){
    final_success <- paste(descr1, toupper(combop), " ", descr2, sep="")}
  else
  {
    final_success <- paste("at least ", success*100, "% of time between level(s) ", paste(targets1, collapse = ", "), " and ", paste(targets2, collapse = ", "), sep = "")
  }

  if (!silent) {
    cat("\nFinal success will be: ", final_success, sep="")
  }
  
  attr(rval, "success") <- final_success

  return(rval)
  
}
