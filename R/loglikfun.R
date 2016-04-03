LogLikFun <- function(params,X,y,wt,offset,link){
#   if(!is.null(dim(y))){
    lp <- X %*% c(params[1],exp(params[-1])) + offset
    if(link=='negexp') p <- bbgdm::negexp()
    else  p <- make.link(link = link) 
    lb <- dbinom(y[,1], y[,2], p$linkinv(lp),log = TRUE)*wt
    ll.contr<-sum(lb)
  return(ll.contr)
}  
