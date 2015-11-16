#' Simulate number of non-shared species across a gradient. 
#' 
#' @param form formula for simulation.
#' @param dat data to simulate dissimilarites.
#' @param theta response parameters.
#' @param nsp Number of species underlying simulation.
#' @return list containing, form,dat,theta,nsp,variance and max_sp_site.
#' @export
#' @examples
#' library(bbgdm)
#' theta <- matrix(c(.5,2),1,2)
#' form <- y~1+x
#' nsp <- 250
#' Nsite <- 50
#' ndissim <- (Nsite*(Nsite-1))/2
#' dat <- data.frame(y=rep(1,ndissim),x=seq(0,4,length.out = ndissim))#imaginary difference in a covariate
#' sim1 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),.1),1,2),nsp)
#' plot(sim1$X[,2],sim1$dism,ylim=c(0,1),pch=16,cex=.75)

simulate_dism_data <- function(form,dat,theta,nsp){
  max_sp_site <- round(runif(dim(dat)[1],50,nsp))
  X <- model.matrix(form, dat)
  lgtp <- X %*% drop(theta)
  p <- exp(lgtp)/(1 + exp(lgtp))
  sim_nnss <- rbinom(dim(X)[1], max_sp_site, p)
  sim_dism <- sim_nnss/max_sp_site
  vars <- (sim_nnss*p)*(1-p)
  return(list('X'=X,'nnss'=sim_nnss,'dism'=sim_dism,'max_sp_site'= max_sp_site,'theta'=theta,'vars'=vars,'y'=data.frame(SharedSpp_ij=sim_nnss, MaxSpp_ij=max_sp_site)))
}