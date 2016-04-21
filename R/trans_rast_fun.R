#' Create transition raster
#'
#'@description Creates a transition raster from matrix or raster to be used to compute least cost distances between sites.
#'@usage trans_rast_fun(data,min.range=0,max.range=NULL)
#'@param data A matrix of lons as rows, lats as cols and value as gradient (eg depth or elevation), or a raster
#'@param min.range  Numeric. The range of depth between which the path will be possible. The default (min.depth=0 and max.depth=NULL) indicates that the transition between cells of the grid is possible between 0 meters depth and the maximum range.
#'@param max.range	Numeric. The range of depth between which the path will be possible. The default (min.depth=0 and max.depth=NULL) indicates that the transition between cells of the grid is possible between 0 meters depth and the maximum range.
#'@return A transition raster for least_cost_dist function.
#'@note This function is slow, use courser resolution lons x lat or raster grids to speed up.
#'@export

trans_rast_fun <- function (data, min.range = 0, max.range = NULL)
{ if(!is.raster(data)){
  ras <- data
  ras[data > min.range] <- 1e-08
  ras[data <= min.range] <- 1
  if (!is.null(max.range))
    ras[data <= max.range] <- 1e-08
  lat <- as.numeric(colnames(data))
  lon <- as.numeric(rownames(data))
  r <- raster::raster(ncol = nrow(data), nrow = ncol(data),
                      xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  raster::values(r) <- as.vector(ras[, rev(1:ncol(ras))])
  trans <- gdistance::transition(r, transitionFunction = mean,
                                 directions = 16)
  transC <- gdistance::geoCorrection(trans, type = "c", multpl = FALSE)
  return(transC)
} else {
  r <- data
  trans <- gdistance::transition(r, transitionFunction = mean,
                                 directions = 16)
  transC <- gdistance::geoCorrection(trans, type = "c", multpl = FALSE)
  return(transC)
  }
}
