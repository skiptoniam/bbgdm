LogLikFun <- function(params,X,y,wt,offset){
  if(!is.null(dim(y))){
    paramsTrans <- c(params[1],exp(params[-1]))
    lp <- X %*% paramsTrans + offset
    mu <- plogis(lp)
    lprob <- dbinom( y[,1], y[,2], mu, log=TRUE)
    ll.contr <- sum(wt * lprob)
    } else {
    paramsTrans <- c(params[1], exp( params[-1]))
    lp <- X %*% paramsTrans + offset
    mu <- plogis(lp)
    lprob <- dbinom(0, 1, mu, log=TRUE)
    lprob <- wt *ifelse(y,log(mu),log(1-mu))
    ll.contr <- sum(lprob)
#     ll.contr <- sum(wt*ifelse(y,log(mu),log(1-mu)))
    }
  return(ll.contr)
}  
