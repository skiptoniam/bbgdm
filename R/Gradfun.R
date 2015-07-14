GradFun <- function(params,X,y,wt,offset){
  if(!is.null(dim(y))){
    paramsTrans <- c(params[1], exp( params[-1]))
    lp <- X %*% paramsTrans + offset
    np<-length(params)
    ns<-nrow(X)
    deri<- matrix(NA,ns,np)
    for(i in 1:ns){
      mu <- y[i,2]*exp(lp[i])/(1+exp(lp[i]))
      dldm <- y[i,1]/mu - (1-y[i,1]) / (1-mu)
      dmde <- mu*(1-mu)
      dedb <- X[i,]
      dbdg <-  c(params[1], exp( params[-1]))
      deri[i,] <- wt[i] * dldm * dmde * dedb * dbdg
    }
    sum_deri <- apply(deri,2,sum)
  } else {
    paramsTrans <- c(params[1], exp( params[-1]))
    lp <- X %*% paramsTrans + offset
    mu <- plogis(lp)
    sum_deri <- drop(wt*dlogis(lp)*ifelse(y,1/mu,-1/(1-mu)))%*%X
#   np<-length(params)
#   ns<-nrow(X)
#   deri<- matrix(NA,ns,np)
#     for(i in 1:ns){
#       mu <- exp(lp[i])/(1+exp(lp[i]))
#       dldm <- y[i]/mu - (1-y[i]) / (1-mu)
#       dmde <- mu*(1-mu)
#       dedb <- X[i,]
#       dbdg <- c(params[1], exp( params[-1]))
#       deri[i,] <- wt[i] * dldm * dmde * dedb * dbdg
#     }
#     sum_deri <- apply(deri,2,sum)
  }
  return(sum_deri)
}