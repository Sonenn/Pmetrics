#' Write session defaults to disk so that they are persistent.
#'
#' The default behavior of any Pmetrics function can be set by the user for a session with \code{\link{setDefaults}}.  However, these
#' custom defaults will not be saved.  Use this function to write them to the Pmetrics library for session-to-session permanence.
#'
#' @title Save Pmetrics Session Defaults
#' @return None
#' @author Michael Neely
#' @seealso \code{\link{setDefaults}}, \code{\link{getDefaults}}, \code{\link{unsetDefaults}}

PMwriteDefaults <- function(){
  
  PMopt <- list()
  nopt <- length(getDefaults())
  if(nopt>0){
    for (i in 1:nopt){
      PMopt[i] <- list(getDefaults(getDefaults()[i]))
      names(PMopt)[i] <- paste(getDefaults()[i],"Default",sep=".")
    }
  }else { PMopt <- NULL}
  optFile <- paste(get("PmetricsPath",envir=PMenv),"/Pmetrics/config/PMopt.Rdata",sep="")
  save(PMopt,file=optFile)
  
}

