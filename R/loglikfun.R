LogLikFun <- function(params,X,y,wt,offset){
  if(!is.null(dim(y))){
    lp <- X %*% c(params[1],exp(params[-1])) + offset
    p <- make.link(link = "logit")
    lb <- dbinom(y[,1], y[,2], p$linkinv(lp),log = TRUE)*wt
    ll.contr<-sum(lb)
    } else {
    paramsTrans <- c(params[1], exp( params[-1]))
    lp <- X %*% paramsTrans + offset
    mu <- plogis(lp)
    lprob <- dbinom(0, 1, mu, log=TRUE)
    lprob <- wt *ifelse(y,log(mu),log(1-mu))
    ll.contr <- sum(lprob)
    }
  return(ll.contr)
}  
