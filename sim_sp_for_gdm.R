# generate thetas
require(mvtnorm)
library(MASS)
# generate "species" with responses to covariates... Create the underlying hypotheses..
# and to mix it up (skew the distributions).. Mixed responses to models..
# mean1 <- c(1.2,.5,1) #for intercepts and slopes # a hypothesis on the response of each parameter...
mean1 <- c(.5,.5,-.25) #for intercepts and slopes # a hypothesis on the response of each parameter...
mean2 <- c(.5,-1.2,-1.44)
means <- rbind(mean1,mean2)
variances <- matrix(c(rep(1.5,dim(means)[2]),rep(2,dim(means)[2])),dim(means)[1],dim(means)[2],byrow=T)
covariances <- apply(variances,1,function(x)1/(sqrt(mean(x))*sqrt(mean(x))))
mix.prop <- 0.35
nSp <- 250   #lots of species (to see pattern in distribution)
nSites <- 50
dat <- data.frame(y=rep(1,nSp),x=seq(-4,4,length.out = nSites),z=runif(nSites,2,6))
# dat$x2 <- dat$x^2
species_theta_generation_mixed_proportions <- function(means,variances,covariances,
                                                       nSp,dat,mix.prop,dist='negbin',
                                                       phi=NULL,plot=TRUE){
  nSigma <- dim(means)[1]
  dimSigma <- dim(means)[2]
  sigmas <- array(0L,c(dimSigma,dimSigma,nSigma))  
  species_parameters <- array(0L,c(nSp,dimSigma,nSigma))
  for(i in 1:nSigma){
    tmp <-sigmas[,,i] 
    tmp <- diag(variances[i,])
    tmp[row(tmp)!=col(tmp)] <- covariances[i]
    sigmas[,,i]<-tmp
    species_parameters[,,i] <- mvtnorm::rmvnorm(nSp, means[i,], sigmas[,,i])
  }
  grp <- rbinom(nSp, c(nSigma)-1, mix.prop)
  df <- data.frame(grp=grp,do.call(cbind,plyr::alply(species_parameters,3)))
  theta <- matrix(0L,nSp,dimSigma)
  for(i in 1:nSigma){
    row_tmp <- which(grp == c(i-1))
    theta[row_tmp,] <- as.matrix(subset(df, grp == c(i-1), select=c(c(i+1):c(i+1+nSigma))))
  }
  out <- matrix(0, dim(dat)[1], nSp)
  X <- as.matrix(dat)
  for (s in 1:nSp){
    if (dist == "bernoulli") {
      lgtp <- X %*% theta[s, ]
      p <- exp(lgtp)/(1 + exp(lgtp))
      out[, s] <- rbinom(dim(X)[1], 1, p)
    }
    if (dist == "negbin") {
      if(is.null(phi))phi<-rep(1,nSp)
      tmp <- rep(1e+05, dim(X)[1])
      while (max(tmp, na.rm = T) > 5000 | sum(tmp) < 2) {
        lpd_sp <- X %*% theta[s,]  
        p <- exp(lpd_sp)
        tmp<-rnbinom(dim(X)[1],mu=p,size=1/phi[s])
      }
      out[,s] <- tmp
    }
    if (dist == "gaussian") {
      tmp <- rep(1e+05, dim(X)[1])
      while (max(tmp, na.rm = T) > 50000 | sum(tmp) < 100) {
        lgtp <- X %*% theta[s, ]
        p <- (lgtp)
        tmp <- rnorm(dim(X)[1], mean = p, sd = 1)
      }
      out[, s] <- tmp
    }
  }
  if(plot)for(i in 1:dim(theta)[2]-1)filled.contour(MASS::kde2d(theta[,1], theta[,i+1]), xlab="Intercepts", ylab=paste0("Slopes",i))
  return(list(sp_data=out, thetas=theta, mu=means,sigma=sigmas))
}

sim_data <- species_theta_generation_mixed_proportions(means,variances,covariances,nSp,dat,mix.prop,dist='bernoulli',plot=TRUE)
sim_data$sp_data
head(dat)
form <- ~ 1 + x + z
fm1 <- bbgdm::gdm.bb(form,sim_data$sp_data,dat,nboot = 100,geo = FALSE)





