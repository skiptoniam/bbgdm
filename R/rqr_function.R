#'Randomised quantile residual function
#'Computes randomised quantile resiudeals for bbgdm
#'@param X model matrix
#'@param y observed values
#'@param coefs parameter estimates
#'@return res randomised quantile residuals


rqr_function <- function(X,y,coefs,offset,link){ 
  lp <- (X %*% coefs + offset)
  if(link=='negexp')link.fun <- bbgdm::negexp()
  else link.fun <- make.link(link=link)
  if(!is.null(dim(y))){
    n <- rep(1, length(y[,1]))
    y[,1] <- n * y[,1]
    a <- pbinom(y[,1]-1, y[,2], p$linkinv(lp))#-1
    b <- pbinom(y[,1], y[,2], p$linkinv(lp))
    u <- runif(n = length(y[,1]), min = a, max = b)
    res <- qnorm(u)
#     res[res==Inf]<-max(res[res!=Inf])
#     res[res==-Inf]<-min(res[res!=-Inf])
    } else {
  n <- rep(1, length(y))
  y <- n * y
  a <- pbinom(y-1, n, pi)
  b <- pbinom(y, n, pi)
  u <- runif(n = length(y), min = a, max = b)
  res <- qnorm(u)
#   res[res==Inf]<-max(res[res!=Inf])
#   res[res==-Inf]<-min(res[res!=-Inf])
  }
  return(res)
}
