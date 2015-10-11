#' Transformation of Covariates to Splines. 
#' 
#' Splines to be used in bbGDM; bsplines and isplines can be calculated.
#' @param x data.frame of covariates to use as predictors in bbGDM.
#' @param spline_type If "bspline" calculates bs spline from spline package. If "ispline" calculates ispline.
#' @param spline_df degrees of freedom; one can specify df rather than knots. Default = 3.
#' @param spline_knots The internal breakpoints that define the spline. Default = NULL for bs spline and 1 for ispline.
#' @return A matrix of dimension c(length(x), df), where df was supplied.
#' @export


spline.trans <- function(x,spline_type ="bspline",spline_df=2,spline_knots=1){
  if(!is.data.frame(x)) x <- as.data.frame(x)
#   stop("Covariates aren't a data.frame object")
  require(splines)
  index <- 1
  spline.attr <- list()
#   n_spl <- ncol(spline_knots)+(spline_df-2)
#   spline <- matrix(NA,nrow(x),ncol(x)*n_spl)
  if(spline_type=="bspline"){
    if(is.null(spline_knots)) spline_knots <- t(apply(x,2,function(x)quantile(x,c(.5))))
    tmp <- bs(x[,1],degree=spline_df,knots=spline_knots[1])
    n_spl <- dim(tmp)[2]
    spline <- matrix(NA,nrow(x),ncol(x)*n_spl)
      for(ii in 1:ncol(x)){  
      tmp.st <- bs(x[,ii],degree=spline_df,knots=spline_knots[ii])
      spline.attr[[ii]] <- attributes(tmp.st)
      spline[,index:(index+(spline.attr[[ii]]$dim[2]-1))] <- tmp.st
      index <- index + spline.attr[[ii]]$dim[2]
    }
  } 
  if(spline_type=="ispline") {
    if(is.null(spline_knots)) spline_knots <- 1#t(apply(x,2,function(x)quantile(x,c(0,0.25,.5,.75,1))))
    tmp <- ispline(x[,1],spline.knots=spline_knots,spline.degree = spline_df)
    n_spl <- dim(tmp)[2]
    spline <- matrix(NA,nrow(x),ncol(x)*n_spl)
    for(ii in 1:ncol(x)){
      tmp.st <- ispline(x[,ii],spline.knots=spline_knots,spline.degree = spline_df)
      spline.attr[[ii]] <- attributes(tmp.st)
      spline[,index:(index+(spline.attr[[ii]]$dim[2]-1))] <- tmp.st
      index <- index + spline.attr[[ii]]$dim[2]
    }
  }
    colnames(spline) <- do.call(paste, c(as.list(rev(expand.grid(1:n_spl, colnames(x)))), sep='_'))
  return(list(spline=spline,spline.attr=spline.attr))
}