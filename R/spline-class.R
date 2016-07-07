#' @title spline objects
#' @title spline objects
#' @rdname spline
#' @description create and transform covariates to isplines or bsplines from use in \code{bbgdm}.
#' @name ispline
#' @param x data.frame of covariates to use as predictors in bbgdm.
#' @param spline.knots number knots.
#' @param knots NULL degrees of freedom; one can specify df rather than knots. Default = 3.
#' @param spline.degree degrees of freedom for ispline.
#' @return A matrix of dimension c(length(x), df), where df was supplied.
#' @export

ispline <- function (x, spline.knots = 2, knots = NULL, spline.degree = 3)
{
  if (is.null(knots)) {
    quantiles <- seq(0, 1, length = spline.knots + 2)
    interval <- quantile(x, probs = quantiles, names = FALSE)
    knots <- unique(interval)
  }
  if (length(knots) < 2) {
    I <- structure(normalise(x),
                   dim = c(length(x), 1))
    return(I)
  }
  else {
    spline.degree <- min(length(knots) - 1, spline.degree)
    isplinebase <- function (x, knots, d)
    {
      if (is.null(knots) || any(is.na(knots)) || any(diff(knots) == 0) || length(knots) <= 2)
        return(x)
      m <- length(knots)
      n <- length(x)
      interval <- findInterval(x, knots, all.inside = TRUE)
      M <- sapply(sequence(m - 1), `==`, interval)
      for (i in 2:(d + 1)) {
        tik <- c(knots[-1], rep(knots[m], i - 2))
        ti <- c(rep(knots[1], i - 2), knots[-m])
        M <- M %*% diag(1/(tik - ti))
        Dx <- Dt <- array(0, dim = c(m + i - 3, m + i - 2))
        Dx[1L + 0L:(m + i - 4L) * (m + i - 2L)] <- -1
        Dx[1L:(m + i - 3L) * (m + i - 2L)] <- 1
        Dt[1L + 0L:(m + i - 4L) * (m + i - 2L)] <- tik
        Dt[1L:(m + i - 3L) * (m + i - 2L)] <- -ti
        M <- (M * x) %*% Dx + M %*% Dt
      }
      M <- M[, -1]
      S <- array(1, dim = rep(NCOL(M), 2))
      S[upper.tri(S)] <- 0
      I <- M %*% S
      return(I)
    }
    I <- structure(isplinebase(x, knots, spline.degree),
                   splineInterval = knots, splineDegree = spline.degree)
    return(I)
  }
}

#' @title spline objects
#' @rdname spline
#' @name spline.trans
#' @param spline_type If "ispline" calculates bs spline from spline package. If "bspline" calculates ispline.
#' @param spline_df degrees of freedom; one can specify df rather than knots. Default = 3.
#' @param spline_knots The internal breakpoints that define the spline. Default = NULL for bs spline and 1 for ispline.
#' @return A matrix of dimension c(length(x), df), where df was supplied.
#' @export

spline.trans <- function(x,spline_type ="ispline",spline_df=2,spline_knots=1){
  if(!is.data.frame(x)) x <- as.data.frame(x)
  index <- 1
  spline.attr <- list()
  if(spline_type=="bspline"){
    spline_knots <- t(apply(x,2,function(x)quantile(x,c(seq(0,1,1/(1+spline_knots))[-c(1,length(seq(0,1,1/(1+spline_knots))))]))))
    tmp <- splines::bs(x[,1],degree=spline_df,knots=spline_knots[1])
    n_spl <- dim(tmp)[2]
    spline <- matrix(NA,nrow(x),ncol(x)*n_spl)
    for(ii in 1:ncol(x)){
      tmp.st <- splines::bs(x[,ii],degree=spline_df,knots=spline_knots[ii])
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
      spline.attr[[ii]]$class <- 'is'
    }
  }
  colnames(spline) <- do.call(paste, c(as.list(rev(expand.grid(1:n_spl, colnames(x)))), sep='_'))
  structure(list(spline=spline,spline.attr=spline.attr),class='spline')
}

#' @rdname dissimilarity
#' @param object spline class object
#' @export
is.spline <- function (object) {
  # test whether object is a dissim_table object
  ans <- inherits(object, "spline")
  # return the answer
  return (ans)

}

#' @rdname spline
#' @name spline_trans_for_pred
#' @param attrib Spline parameters from model.
#' @param splineInterval NULL
#' @param splineDegree NULL
#' @return I-Spline predictions
#' @export

spline_trans_for_pred <- function(x, attrib = NULL, splineInterval = NULL, splineDegree = NULL)
{
  if((attrib$class)[1]=='bs'){
    if (is.list(attrib)) {
      #       values = attrib$values
      #       standardization = attrib$standardization
      splineInterval = attrib$knots
      splineDegree = attrib$degree
      Bound.knots = attrib$Boundary.knots
    }
    if (!is.null(splineInterval))
      return(splines::bs(x, knots = splineInterval, degree = splineDegree,Boundary.knots=Bound.knots))
    # else if (!is.null(standardization)) {
    #   I <- normalise(x, standardize = standardization)
    #   attr(I, "dim") <- c(length(x), 1)
    #   return(I)
    # }
    # else if (!is.null(values)) {
    #   x <- factor(x, levels = values)
    #   I <- model.matrix(~x - 1)
    #   attr(I, "values") <- values
    #   return(I)
    # }
    # else {
    #   I <- as.numeric(x)
    #   attr(I, "dim") <- c(length(x), 1)
    #   return(I)
    # }

  } else{
    if (is.list(attrib)) {
      values = attrib$values
      standardization = attrib$standardization
      splineInterval = attrib$splineInterval
      splineDegree = attrib$splineDegree
    }
    if (!is.null(splineInterval))
      return(ispline(x, knots = splineInterval, spline.degree = splineDegree))
    # else if (!is.null(standardization)) {
    #   I <- normalise(x, standardize = standardization)
    #   attr(I, "dim") <- c(length(x), 1)
    #   return(I)
    # }
    # else if (!is.null(values)) {
    #   x <- factor(x, levels = values)
    #   I <- model.matrix(~x - 1)
    #   attr(I, "values") <- values
    #   return(I)
    # }
    # else {
    #   I <- as.numeric(x)
    #   attr(I, "dim") <- c(length(x), 1)
    #   return(I)
    # }
  }
}


normalise <- function (x)
{
  stand.a <- min(x)
  stand.b <- max(x) - stand.a
  X <- (x - stand.a)/stand.b
  attr(X, "standardization") <- list(a = stand.a, b = stand.b)
  return(X)
}
