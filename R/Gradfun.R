GradFun <- function(params,X,y,wt,offset,link){
#   if(!is.null(dim(y))){
    lp <- X %*% c(params[1],exp(params[-1])) + offset
    np<-length(params)
    ns<-nrow(X)
    deri<- matrix(NA,ns,np)
    if(link!='negexp')p <- negexp()
    else  p <- make.link(link = link) 
    for(i in 1:ns){
      mu <- y[i,2]*p$linkinv(lp[i])
      dldm <- (y[i,1]/mu) - ((1-y[i,1])/(1-mu))
      dmde <- mu*(1-mu)
      dedb <- X[i,]
      dbdg <- c(1,exp(params[-1])) 
      deri[i,] <- wt[i] * dldm * dmde * dedb * dbdg
    }
    sum_deri <- apply(deri,2,sum)
    return(sum_deri)
}