#' Spatiall predict a bayesian bootstrap generalised dissimilarity model.
#' @param object As derived from bbgdm function
#' @param data raster stack of the same covariates used to fit model object for the region you wish to predict too.
#' @param neigbourhood int default is three, number of neighbouring cells to estimate mean dissimilarity.
#' @param outer logical default is FALSE, if TRUE only calculates the outer edge of neighbourhoods area.
#' @param uncertainty logical if TRUE predict will return a list with two rasters the mean estimate and the uncertainty (defined as the coefficent of variation)
#' @return raster of mean turnover estimated based on neighbourhood distance.
#' @export

predict.bbgdm <- function (object, data, neighbourhood=NULL, outer=FALSE, uncertainty=TRUE,...)
{
  options(warn.FPU = FALSE)
    if (class(data) != "RasterStack" & class(data) != "RasterLayer" &
        class(data) != "RasterBrick") {
      stop("Prediction data need to be a raster object")
    }
    if (nlayers(data) != ncol(object$env.dat)) {
      stop("Number of raster layers does not equal the number used to fit the model")
    }
    for (i in 1:nlayers(data)) {
      if (names(object$env.dat)[i] != names(data)[i]) {
        stop("Raster layers don't match variables used to fit the model - check they are in the correct order")
      }
    }
  if(is.null(neighbourhood)){
    cat('using default three cell neighbourhood to estimate dissimilarity')
    neighbourhood <- 3
  }
  #create neighbour matrix
  rel<-expand.grid(seq(-neighbourhood,neighbourhood,1),seq(-neighbourhood,neighbourhood,1))
  if (outer) rel <-rel[(rel[,1]==-neighbourhood | rel[,1]==neighbourhood | rel[,2]==-neighbourhood | rel[,2]==neighbourhood),]
  rel<-as.matrix(rel[,c(2,1)])
  max.dist<-max(rel)
  rel<-cbind(rel,(1-((sqrt((rel[,1]^2)+(rel[,2]^2))/max.dist))))

  #create updated i-spline rasters (transform layers to i-spline values)
  XYdata <- as.data.frame(na.omit(rasterToPoints(data,progress = "text")))
  cells <- cellFromXY(data[[1]], cbind(XYdata$x, XYdata$y))
  data_spline <- spline.trans(XYdata[,-1:-2],spline_type = "ispline",spline_df=2)
  st_data <- as.data.frame(cbind(XYdata[,1:2], data_spline$spline))
  coordinates(st_data) <- ~x+y
  suppressWarnings(gridded(st_data) <- TRUE)
  data_stack <- stack(st_data)

  beta.r <- data[[1]]
  beta.r[cells]<-1
  intercept <- beta.r
  rs <- stack(intercept,data_stack)
  bbgdm_coef <- as.vector(object$median.coefs.se)
  raster_data <- raster::as.array(rs)
  beta_mat <- raster::as.matrix(beta.r)
  beta.r[] <- bbgdm::pred_bbgdm_cpp(raster_data,beta_mat,rel,bbgdm_coef)
  if(uncertainty){
  beta.r.lw <- data[[1]]
  beta.r.up <- data[[1]]
  beta.r.lw[cells]<-1
  beta.r.up[cells]<-1
  bbgdm_coef.lw <- as.vector(object$quantiles.coefs.se[1,])
  bbgdm_coef.up <- as.vector(object$quantiles.coefs.se[2,])
  beta.r.lw[] <- bbgdm::pred_bbgdm_cpp(raster_data,beta_mat,rel,bbgdm_coef.lw)
  beta.r.up[] <- bbgdm::pred_bbgdm_cpp(raster_data,beta_mat,rel,bbgdm_coef.up)
  pred.se <- abs(beta.r.up-beta.r.lw)/beta.r
  }
  if(uncertainty) return(list(mean.beta=beta.r,se.beta=pred.se))
  else return(beta.r)
}
