# generate thetas
library(mvtnorm)
library(MASS)
# generate "species" with responses to covariates... Create the underlying hypotheses..
# and to mix it up (skew the distributions).. Mixed responses to models..
# mean1 <- c(1.2,.5,1) #for intercepts and slopes # a hypothesis on the response of each parameter...
mean1 <- c(1.5,1.5) #for intercepts and slopes # a hypothesis on the response of each parameter...
mean2 <- c(-1.5,-1.5)
means <- rbind(mean1,mean2)
variances <- matrix(c(rep(.56,dim(means)[2]),rep(.56,dim(means)[2])),dim(means)[1],dim(means)[2],byrow=T)
variances <- t(variances)
covariances <- c(-.50,-.50) #apply(variances,1,function(x)1/(sqrt(mean(x))*sqrt(mean(x))))
mix.prop <- 0.5
nSp <- 250   #lots of species (to see pattern in distribution)
nSites <- 30
dat <- data.frame(y=rep(1,nSites),x=seq(-2,2,length.out = nSites))
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
    theta[row_tmp,] <- as.matrix(subset(df, grp == c(i-1), select=c(c(i+1):c(i+nSigma))))
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
      while (max(tmp, na.rm = T) > 50000 | sum(tmp) < 2) {
        lpd_sp <- X %*% theta[s,]  
        p <- exp(lpd_sp)
        tmp<-rnbinom(dim(X)[1],mu=p,size=1/phi[s])
      }
      out[,s] <- tmp
    }
    if (dist == "poisson") {
      tmp <- rep(1e+05, dim(X)[1])
      while (max(tmp, na.rm = T) > 50000 | sum(tmp) < 2) {
        lpd_sp <- X %*% theta[s,]  
        p <- exp(lpd_sp)
        tmp<-rpois(dim(X)[1],lambda =p)
      }
      out[,s] <- tmp
    }
  }
  if(plot)for(i in 1:dim(theta)[2]-1)filled.contour(MASS::kde2d(theta[,1], theta[,i+1]), xlab="Intercepts", ylab=paste0("Slopes",i))
  return(list(sp_data=out, thetas=theta, mu=means,sigma=sigmas))
}

set.seed(42)
sim_data <- species_theta_generation_mixed_proportions(means,variances,covariances,nSp,dat,mix.prop,dist='poisson',plot=TRUE)
sim_data$sp_data

head(dat)
form <- ~ 1 + x
library(bbgdm)
fm1 <- bbgdm::gdm.bb(form,sim_data$sp_data,dat,nboot = 100,geo = FALSE,optim.meth = 'admb',control=logit_glm_control(start=rbeta(4,2,2)))
bbgdm::bb.gdm.check(fm1)
bbgdm::plotResponse(fm1)
bbgdm::bbgdm.wald.test(fm1)
bbgdm::plot.gdm.bb(fm1)

library(plyr)
bbgdm_coefs_sim <- fm1$all.coefs.se
link.fun <- make.link("logit")
par(mfrow=c(1,2))
x <- seq(0,4,length.out = 105)                                                
xs <- bbgdm::spline.trans(x,spline_type = 'ispline')

x1<-diff_table_cpp(as.matrix(dat[,2]))


plot(x1,fm1$starting_gdm$y[,1]/fm1$starting_gdm$y[,2],ylim=c(0.2,.8),xlim=c(0,4),col = "grey60",cex=.6,pch=16,ylab='dissimilarity')
link.fun <- make.link("logit")
all_results <- matrix(0,100,105)
for(i in 1:nrow(bbgdm_coefs_sim)) {
  lines(x, link.fun$linkinv(cbind(1,xs$spline)%*%as.vector(as.matrix(bbgdm_coefs_sim[i,]))), col = "#00000030")
  # lines(x, link.fun$linkinv(cbind(1,xs$spline)%*%as.vector(as.matrix(bbgdm_coefs_var_sim[i,1:4]))), col = "#00000030")
  # lines(x, link.fun$linkinv(cbind(1,xs$spline)%*%as.vector(as.matrix(bbgdm_coefs_var_sim[i,5:8]))), col = "#00000030")
  all_results[i,] <- link.fun$linkinv(cbind(1,xs$spline)%*%as.vector(as.matrix(bbgdm_coefs_sim[i,])))
}
lines(x, apply(all_results,2,mean), col = "red",lwd=2)
