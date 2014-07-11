#' Print a Pmetrics Error Object
#'
#' Print the errors in a Pmetrics data file or PMmatrix object.
#'
#' @title Print Data Errors
#' @method print PMerr
#' @param x A PMerr object made by \code{\link{PMcheckMatrix}}.
#' @param \dots Other parameters which are not necessary.
#' @return A printed object.
#' @author Michael Neely
#' @seealso \code{\link{PMcheckMatrix}}
#' @export

print.PMerr <- function(x,...){
  cat("\n")
  for(i in 1:length(x)){
    if(is.na(x[[i]]$results[1])) {
      cat(paste("(",i,") ",x[[i]]$msg,"\n",sep=""))
    } else {
      if(i==16){cat(paste("\n(",i,") ",x[[i]]$msg," ",paste(x[[i]]$results,collapse=", "),"\n\n",sep=""))
      } else {cat(paste("\n(",i,") ",x[[i]]$msg," Access the rows in your data with 'data[err$",names(x)[i],"$results,]', replacing 'err' with the name of this error object.\n",paste(x[[i]]$results,collapse=", "),"\n\n",sep=""))}
    }
  }
  
}
