#' @title
#' greater circle dist
#' 
#' @description
#' Takes in any numeric value and squares itFunction for calculating greater circle dist
#' 
#'@param x1 first longitude
#'@param x2 second longitude
#'@param y1 first latitude
#'@param y2 second latitude
#'@return greater circle distance. 
#'@export

gcdist <- function(x1, y1, x2, y2) {
  ##############################################################################################
  #a function called by several functions in the ncf library
  #Function for great circle distance -- due to T. Keitt.
  #See http://www.census.gov/cgi-bin/geo/gisfaq?Q5.1
  ##from package ncf
  ##############################################################################################
  r <- 360/(2 * pi)
  lon1 <- x1 / r
  lat1 <- y1 / r
  lon2 <- x2 / r
  lat2 <- y2 / r
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  a <- (sin(dlat/2))^2 + cos(lat1) * cos(lat2) * (sin(dlon/2))^2
  c <- 2 * atan2( sqrt(a), sqrt(1-a) )
  return(6370 * c)
}