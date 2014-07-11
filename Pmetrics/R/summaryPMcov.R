#' Summarize a Pmetrics Covariate object
#'
#' Summarize covariates and Bayesian posterior parameter values for each subject.
#'
#' @title Summarize Covariates and Bayesian Posterior Parameter Values
#' @method summary PMcov
#' @param x A PMop object made by \code{\link{makeOP}}.
#' @param icen Summary function for covariates and posterior parameters. Default is \dQuote{median}, but can specify \dQuote{mean}.
#' @return A data frame with the summary of the PMcov object for each subject's covariates and 
#' Bayesian posterior parameter values.
#' @author Michael Neely
#' @seealso \code{\link{makeCov}}
#' @export

summary.PMcov <- function(x,icen="median"){
  if("icen" %in% names(x)){
    data <- x[x$icen==icen,]
    data <- subset(data,select=-icen) 
  } else {data <- x}
  sumCov <- aggregate(data,list(data$id),match.fun(icen),na.rm=T)
  #remove the first grouping column
  sumCov <- sumCov[,-1]
  attr(sumCov,"icen") <- icen
  return(sumCov)
}

