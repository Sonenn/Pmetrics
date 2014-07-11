#' See changelog for Pmetrics
#'
#' @title Pmetrics changelog
#' @param version Default is the current version, otherwise a character string with the starting version you wish to see up to 
#' the current, e.g. \dQuote{0.21}.
#' @return The changelog for the requested version.
#' @author Michael Neely
#' @examples
#' PMnews()

PMnews <- function(version=packageVersion("Pmetrics")){
  news(Version>=version,package="Pmetrics")
}

