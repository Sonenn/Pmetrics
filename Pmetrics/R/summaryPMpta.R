#' Summarize a Pmetrics Percent Target Attainment Object
#'
#' Summarize target statistics and success proportions in a PMpta object made by \code{\link{makePTA}}.
#'
#' @title Summarize Percent Target Attainment
#' @method summary PMpta
#' @param object A PMpta object made by \code{\link{makePTA}}.
#' @param \dots Other parameters which can be passed to \code{summary}.
#' @return A data frame with the following columns: simnum, target, success, meanratio, and sdratio.  
#' \emph{simnum} is the number of the simulation; \emph{target} is the specified target; 
#' \emph{success}  has the proportion with a ratio > \code{success}; \emph{meanratio} and \emph{sdratio} 
#' are the mean and standard deviation of the target ratios for each simulation and target.  

#' @author Michael Neely
#' @seealso \code{\link{makePTA}}
#' @export

summary.PMpta <- function(object,...){
  if(inherits(object,"PMpta")){
    print(object$outcome)
  } else {print(summary(object))}
}
