#' Function to perform a non-parametric wald-test on bootstrap parameter estimates
#' 
#' @param object Returned mode from \link[bbgdm]{gdm.bb}. 
#' @return A table of Wald-Test statistics
#' @export

# 1) Do BBootstrap to get B samples of the estimates.  Store the sample in a B by p matrix A.
# 2) Calculate variance of estimates, using esti.var <- var( A).
# 2a) Calculate the estimates, either from the mean of the bootstraps or the point estimate from the working model.  Call them beta
# 3) Calculate the Wald statistic as t( beta) %*% solve( esti.var) %*% beta.
# 4) Compare the Wald stat against a Chi-Squared with Q degrees of freedom.
# 5) If the Wald stat is bigger than the 0.95th percentile of the Chi-sq(Q) distribution then term is important.
# 
# Things to note:
# *) This still may not work as 1) we don't know Q, 2) Wald tests are not that great anyway.
# *) R is the identity matrix here.  The notation in the wikipedia article is technically correct, but not illuminating for most readers (it even took me a while).
# *) You cannot trust the hessian as a good measure of variance -- it is based on the working (dissimilarity) model, which we know is wrong for measuring the amount of information.

bbgdm.wald.test <- function(object,H0=0){
  # H0: Hypothesis test = 0
  # IM: Identiy Matrix
  # beta: parameter estimates from model
  # vcov: Variance-covariance matrix estimated from BB
  
  A <- object$all.coefs.se #matrix of B bootstrap coeficient estimates.
  esti.var <- var(A) #make sure this is a matrix
  beta <- object$median.coefs.se #medians of coef estimates
  
  #Intercept
  intercept_IM <- matrix(c(1,rep(0,length(beta)-1)),nrow=1)
  wd_inter <- t(intercept_IM%*%beta-H0) %*% solve(intercept_IM%*%esti.var%*%t(intercept_IM))%*%(intercept_IM%*%beta-H0)
  pval_i = 1-pchisq(wd_inter,1)
  
  #Splines
  splineLength <- sapply(object$dissim_dat_params, `[[`, "dim")[2,]
  val1 <- seq(2,length(beta),splineLength[1])
  val2 <- seq(1+splineLength[1],length(beta),splineLength[1])
  wd_vals <- matrix(NA,length(val1)+1,3)
  wd_vals[1,]<- c(wd_inter,1,pval_i)
  for(i in 1:length(splineLength)){
    w <- splineLength[1]
    L <- matrix(rep(0, length(beta) * w), ncol = length(beta))
    Terms <- seq(val1[i],val2[i],1)
      for (ii in 1:w) L[ii, Terms[ii]] <- 1
      wd <- t(L %*% beta - H0) %*% solve(L %*% esti.var %*% t(L)) %*% (L %*% beta - H0)
      pv <- 1 - pchisq(wd, df = w)
      wd_vals[1+i,]<- c(wd,w,pv)
    }
  colnames(wd_vals) <- c("W","df","p-value")
  if(object$geo){ rownames(wd_vals)<-c('intercept','geo',names(object$env.dat)[-c(1:2)])
  } else { rownames(wd_vals)<-c('intercept',names(object$env.dat))
  }
  return(wd_vals)
}
  
