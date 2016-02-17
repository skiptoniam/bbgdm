#' Function to perform a non-parametric wald-test on bootstrap parameter estimates
#' 
#' @param object Returned model from \link[bbgdm]{bbgdm}.
#' @param H0 A numeric value giving the null hypothesis for the test. Generally zero.
#' @param gdm Logic if true calculates the Wald-test using variance-covariance matrix derived from the hessian matrix. Note: These estimates are probably wrong due to hessian matrix being calulated with respect to the likelihoods.
#' @return A table of Wald-Test statistics
#' @export

bbgdm.wald.test <- function(object,H0=0,gdm=FALSE){
  # H0: Hypothesis test = 0
  # IM: Identiy Matrix
  # beta: parameter estimates from model
  # vcov: Variance-covariance matrix estimated from BB
  
  #for gdm
  if(gdm){
    vcov <-object$starting_gdm$var
    beta <- object$starting_gdm$coef
    #Intercept
    intercept_IM <- matrix(c(1,rep(0,length(beta)-1)),nrow=1)
    wd_inter <- t(intercept_IM%*%beta-H0) %*% solve(intercept_IM%*%vcov%*%t(intercept_IM))%*%(intercept_IM%*%beta-H0)
    pval_i = 1-pchisq(wd_inter,1)
    
    #Splines
    splineLength <- sapply(object$dissim_dat_params, `[[`, "dim")[2,]
    val1 <- seq(2,length(beta),splineLength[1])
    val2 <- seq(1+splineLength[1],length(beta),splineLength[1])
    wd_vals_gdm <- matrix(NA,length(val1)+1,3)
    wd_vals_gdm[1,]<- c(wd_inter,1,pval_i)
    for(i in 1:length(splineLength)){
      w <- splineLength[1]
      L <- matrix(rep(0, length(beta) * w), ncol = length(beta))
      Terms <- seq(val1[i],val2[i],1)
      for (ii in 1:w) L[ii, Terms[ii]] <- 1
      wd <- t(L %*% beta - H0) %*% solve(L %*% vcov %*% t(L)) %*% (L %*% beta - H0)
      pv <- 1 - pchisq(wd, df = w)
      wd_vals_gdm[1+i,]<- c(wd,w,pv)
    }
    colnames(wd_vals_gdm) <- c("gdm_W","gdm_df","gdm_p-value")
    if(object$geo){ rownames(wd_vals_gdm)<-c('intercept','geo',names(object$env.dat)[-c(1:2)])
    } else { rownames(wd_vals_gdm)<-c('intercept',names(object$env.dat))
    }
    
  }
    
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
  colnames(wd_vals) <- c("bbgdm_W","bbgdm_df","bbgdm_p-value")
  if(object$geo){ rownames(wd_vals)<-c('intercept','geo',names(object$env.dat)[-c(1:2)])
  } else { rownames(wd_vals)<-c('intercept',names(object$env.dat))
  }
  if(gdm) return(cbind(wd_vals_gdm,wd_vals))
  else return(wd_vals)
}
  
