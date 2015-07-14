#' Function that calculates least cost pathway.
#' 
#' @usage least.cost.dist(trans.mat,sites)
#' @param trans.rast  transition object as computed by trans.mat
#' @param sites A two-columns matrix or data.frame containing latitude and longitude for 2 or more locations.
#' @details least.cost.dist computes least cost distances between 2 or more locations. This function relies on the package gdistance (van Etten, 2011. http://CRAN.R-project.org/package=gdistance) and on the trans.mat function to define a range of depths where the paths are possible.
#' @return distance matrix in meters between all possible pairs of locations
#' @export

least.cost.dist <- function (trans.rast, sites) 
{
  cost <- gdistance::costDistance(trans.rast, as.matrix(sites))
  return(round(cost))
}
