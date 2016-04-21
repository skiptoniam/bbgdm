#' log-likelihood function
#' 
#' @param params initial values see start \link[stats]{optim}
#' @param X model matrix
#' @param y response variable
#' @param wt weights
#' @param offset offset vector (currently not used)
#' @param link link function

LogLikFun <- function(params,X,y,wt,offset,link){
    lp <- X %*% c(params[1],exp(params[-1])) + offset
    if(link=='negexp') p <- bbgdm::negexp()
    else  p <- make.link(link = link) 
    lb <- dbinom(y[,1], y[,2], p$linkinv(lp),log = TRUE)*wt
    ll.contr<-sum(lb)
  return(ll.contr)
}  
