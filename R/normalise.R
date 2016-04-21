#' Function to normalise isplines. 
#' prediction spline transform
#' Transform spline data based on parameters
#' @param x spline data
#' @param standardize 'zscore' 

normalise <- function (x, standardize = "zscore") 
{
  if (is.list(standardize)) {
    stand.a <- standardize$a
    stand.b <- standardize$b
  }
  else if (standardize == "zscore") {
    stand.a <- mean(x)
    stand.b <- sd(x)
    stand.b[stand.b == 0] <- 1
  }
  else if (standardize == "interval") {
    stand.a <- min(x)
    stand.b <- max(x) - stand.a
  }
  else {
    stand.a <- 0
    stand.b <- 1
  }
  if (is.na(stand.b) || stand.b == 0) {
    stand.a <- mean(x)
    stand.b <- 1
  }
  X <- (x - stand.a)/stand.b
  attr(X, "standardization") <- list(a = stand.a, b = stand.b)
  return(X)
}