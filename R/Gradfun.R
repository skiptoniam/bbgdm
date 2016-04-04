GradFun <- function(params,X,y,wt,offset,link){
  #   if(!is.null(dim(y))){
  eta <- X %*% c(params[1],exp(params[-1])) + offset
  np<-length(params)
  ns<-nrow(X)
  deri<- matrix(NA,ns,np)
  if(link=='negexp')  p<- bbgdm::negexp()
  else p <- make.link(link = link)
  for(i in 1:ns){
    mu <- p$linkinv(eta[i])
    dldm <- ((y[i,1]/y[i,2])/mu) - ((1-(y[i,1]/y[i,2]))/(1-mu))
    dmde <- p$mu.eta(eta[i])*y[i,2]
    dedb <- X[i,]
    dbdg <- c(1,exp(params[-1])) 
    deri[i,] <- wt[i] * dldm * dmde * dedb * dbdg
  }
  sum_deri <- apply(deri,2,sum)
  return(sum_deri)
}