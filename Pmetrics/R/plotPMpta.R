#' Plots PMpta objects
#'
#' This function will plot the percent target attainment for objects made with the \code{\link{makePTA}} function.
#' For the legend, defaults that are different that the standard are:
#' \itemize{
#'   \item x Default \dQuote{topright}
#'   \item legend Default will be \dQuote{Dose 1, Dose 2,...Dose n}, where \emph{n} is the number of Dose levels in the PMpta object.  This default
#'   can be overridden by a supplied character vector of dose names.
#'   \item col The color of each Dose plot as specified by the default color scheme or \code{col}
#'   \item pch The plotting character for each Dose plot as specified by the default plotting characters or \code{pch}
#'   \item lty The line type of each Dose plot as specified by the default line types or \code{lty}
#'   \item bg Default \dQuote{white}
#' }
#'
#' @title Plot PMpta Percent Target Attainment objects
#' @method plot PMpta
#' @param x The name of an \emph{PMpta} data object read by \code{\link{makePTA}}
#' @param log Boolean operator to plot in log-log space; the default is \code{True}
#' @param pch Vector of integers which control the plotting symbol for each dose curve; the default is 1:nsim.  NA results in no symbol.
#' Use 0 for open square, 1 for open circle, 2 for open triangle, 3 for cross, 4 for X, or 5 for a diamond.
#' Other alternatives are \dQuote{*} for asterisks, \dQuote{.} for tiny dots, or \dQuote{+} for a smaller,
#' bolder cross.  These plotting symbols are standard for R (see \code{\link{par}}).
#' @param grid Either a boolean operator to plot a reference grid, or a list with elements x and y,
#' each of which is a vector specifying the native coordinates to plot grid lines; the default is \code{False}.
#' For example, grid=list(x=seq(0,24,2),y=1:10).  Defaults for missing x or y will be calculated by \code{\link{axTicks}}.
#' @param xlab Label for the x axis.  Default is \dQuote{MIC}
#' @param ylab Label for the y axis.  Default is \dQuote{Proportion with success}
#' @param col A vector of color names to be used for each dose plotted.  If the
#' length of \code{col} is too short, values will be recycled.
#' @param lty A vector of line types to be used for each dose plotted.  If the
#' length of \code{lty} is too short, values will be recycled.
#' @param lwd Line width, with default of 4.
#' @param legend Either a boolean operator or a list of parameters to be supplied to the \code{\link{legend}}
#' function (see its documentation).  If \code{False}, a legend will not be plotted.
#' If \code{True} (the default), the default legend parameters will be used, as documented in that function, with exceptions
#' as noted in \emph{Details}.
#' @param out Direct output to a PDF, EPS or image file.  Format is a named list whose first argument, 
#' \code{type} is one of the following character vectors: \dQuote{pdf}, \dQuote{eps} (maps to \code{postscript}),
#' \dQuote{\code{png}}, \dQuote{\code{tiff}}, \dQuote{\code{jpeg}}, or \dQuote{\code{bmp}}.  Other named items in the list
#' are the arguments to each graphic device. PDF and EPS are vector images acceptable to most journals
#' in a very small file size, with scalable (i.e. infinite) resolution.  The others are raster images which may be very
#' large files at publication quality dots per inch (DPI), e.g. 800 or 1200. Default value is \code{NA} which means the 
#' output will go to the current graphic device (usually the monitor). For example, to output an eps file,
#' out=list(\dQuote{eps}) will generate a 7x7 inch (default) graphic.
#' @param \dots Other parameters as found in \code{\link{plot.default}}.
#' @return Plots the object.
#' @author Michael Neely
#' @seealso \code{\link{makePTA}}, \code{\link{plot}}, \code{\link{par}}, \code{\link{axis}}
#' @export


plot.PMpta <- function(x,log=T,pch,grid,xlab,ylab,col,lty,lwd=4,legend=T,out=NA,...){

  #choose output
  if(inherits(out,"list")){
    if(out$type=="eps") {setEPS();out$type <- "postscript"}
    if(length(out)>1) {do.call(out$type,args=out[-1])} else {do.call(out$type,list())}
  }
  
  if(!inherits(x,"PMpta")) stop("Please supply a PMpta object made by makePTA().\n")
  if(missing(xlab)) xlab <- "MIC"
  if(missing(ylab)) ylab <- "Proportion with success"
  logscale <- c("","x")[1+as.numeric(log)]
  nsim <- max(x$outcome$simnum)
  if(missing(pch)) {pch <- 1:nsim} else {pch <- rep(pch,nsim)}
  if(missing(col)) {col <- rep(c("black","red","blue","green","purple","orange"),nsim)} else {col <- rep(col,nsim)}
  if(missing(lty)) {lty <- 1:nsim} else {lty <- rep(lty,nsim)}     
  if(class(legend)=="list"){
    legend$plot<- T
    if(is.null(legend$x)) legend$x <- "topright"
    if(is.null(legend$bg)) legend$bg <- "white"
    if(is.null(legend$col)) legend$col <- col
    if(is.null(legend$pch)) legend$pch <- pch
    if(is.null(legend$lty)) legend$lty <- lty
    if(is.null(legend$legend)) legend$legend <- paste("Simulation",1:nsim)
    
  } else {
    if(legend) legend <- {list(plot=T,x="topright",bg="white",col=col,lty=lty,pch=pch,legend=paste("Dose",1:nsim))} else {legend <- list(plot=F)}
  }
  
  plot(prop.success~target,x$outcome,type="n",xlab=xlab,ylab=ylab,log=logscale,xaxt="n",...)
  axis(side=1,at=x$outcome$target,labels=x$outcome$target,lwd=1,...)
  #make grid if necessary
  if(missing(grid)){
    grid <- list(x=NA,y=NA)
  } else {
    if(inherits(grid,"logical")){
      if(grid){
        grid <- list(x=x$outcome$target,y=axTicks(2))
      } else {
        grid <- list(x=NA,y=NA)
      }
    }
    if(inherits(grid,"list")){
      if(is.null(grid$x) | all(!is.na(grid$x))) grid$x <- x$outcome$target
      if(is.null(grid$y) | all(!is.na(grid$y))) grid$y <- axTicks(2)
    }
  }
  abline(v=grid$x,lty=1,col="lightgray")
  abline(h=grid$y,lty=1,col="lightgray")
  
  for(i in 1:nsim){
    lines(prop.success~target,x$outcome[x$outcome$simnum==i,],type="o",xlab=xlab,ylab=ylab,lty=lty[i],lwd=lwd,col=col[i],pch=pch[i],...)   
  }
  if(legend$plot) do.call("legend",legend)
  
  #close device if necessary
  if(inherits(out,"list")) dev.off()

}
