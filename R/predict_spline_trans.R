#' Predict splines transformations
#' @param x list of predictor values
#' @param attrib Spline parameters from model.
#' @param values NULL
#' @param standardization NULL 
#' @param splineInterval NULL
#' @param splineDegree NULL
#' @return I Spline predictions
#' @export

predict.spline.trans <- function (x, attrib = NULL, values = NULL, standardization = NULL, 
          splineInterval = NULL, splineDegree = NULL) 
{
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
