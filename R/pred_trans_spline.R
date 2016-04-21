#' prediction spline transform
#' Transform spline data based on parameters
#' @param x list of predictor values
#' @param attrib Spline parameters from model.
#' @param values NULL
#' @param standardization NULL 
#' @param splineInterval NULL
#' @param splineDegree NULL
#' @return I-Spline predictions
#' @export


spline_trans_for_pred <- function(x, attrib = NULL, values = NULL, standardization = NULL, 
          splineInterval = NULL, splineDegree = NULL) 
{
  if(!is.null(attrib$class)){#[1]=='bs'){
    if (is.list(attrib)) {
#       values = attrib$values
#       standardization = attrib$standardization
      splineInterval = attrib$knots
      splineDegree = attrib$degree
      Bound.knots = attrib$Boundary.knots
    }
    if (!is.null(splineInterval)) 
      return(bs(x, knots = splineInterval, degree = splineDegree,Boundary.knots=Bound.knots))
    else if (!is.null(standardization)) {
      I <- normalize(x, standardize = standardization)
      attr(I, "dim") <- c(length(x), 1)
      return(I)
    }
    else if (!is.null(values)) {
      x <- factor(x, levels = values)
      I <- model.matrix(~x - 1)
      attr(I, "values") <- values
      return(I)
    }
    else {
      I <- as.numeric(x)
      attr(I, "dim") <- c(length(x), 1)
      return(I)
    }
    
  } else{
  if (is.list(attrib)) {
    values = attrib$values
    standardization = attrib$standardization
    splineInterval = attrib$splineInterval
    splineDegree = attrib$splineDegree
  }
  if (!is.null(splineInterval)) 
    return(ispline(x, knots = splineInterval, spline.degree = splineDegree))
  else if (!is.null(standardization)) {
    I <- normalize(x, standardize = standardization)
    attr(I, "dim") <- c(length(x), 1)
    return(I)
  }
  else if (!is.null(values)) {
    x <- factor(x, levels = values)
    I <- model.matrix(~x - 1)
    attr(I, "values") <- values
    return(I)
  }
  else {
    I <- as.numeric(x)
    attr(I, "dim") <- c(length(x), 1)
    return(I)
  }
  }
}
