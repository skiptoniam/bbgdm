#' Function that plots diagnostic plots from bb.gdm object
#' 
#' @param object bb.gdm model output
#' @param rl.col colour of lines
#' @param family logical "binomial" or "beta_dist"
#' @return distances Five column matrix of sites pairs (cols 1:4) and distance as col 5.
#' @export
#' @examples
#' x <- matrix(rbinom(20*10,1,.6),20,10)# presence absence matrix
#' y <- simulate_covariates(x,2)
#' form <- ~ 1 + covar_1 + covar_2
#' test.gdm.bb <- gdm.bb(form,sp.dat=x, env.dat=y,family="binomial", dism_metric="number_shared",nboot=10, scale_covar=T)
#' bb.gdm.check(test.gdm.bb)

bb.gdm.check <- function (object, plots.mfrow = c(2, 2),...) 
{
  family <-object$family
  X <- object$starting_gdm$X
  y <- object$starting_gdm$y
  offset <- object$offset
  coefs<-object$starting_gdm$coef
  par(mfrow=plots.mfrow)
  res <- rqr_function(X,y,coefs,offset)
  qqnorm(res)
  qqline(res, col = 'red')
  hist(res, xlab = "Residuals", main = "Histogram of residuals")
  if(is.null(dim(y))) {
    lp<- X%*%coefs + offset
    pi <- plogis(lp)
    plot(pi, res, main = "Fitted vs. Linear pred.",xlab = "fitted", ylab = "residuals",...)
    plot(y,pi, xlab = "Observed Values", ylab = "Response", 
         main = "Observed vs. Fitted Values",...)
  } else {
    lp<- X%*%coefs + offset
    pi <- plogis(lp)
    plot(pi*y[,2], res, xlab = "Predicted Values", ylab = "Residuals", 
         main = "Fitted Values Vs. Residuals",...)
    plot(y[,1],pi*y[,2], xlab = "Observed Values", ylab = "Predicted Values", 
         main = "Observed Vs. Fitted Values",...)
  }
  mtext("Original GDM diagnostics",outer = TRUE, cex=1.1,col="black",font=2,line=-1)
  ## Now for bayesian boot strap estimates
  par(mfrow=plots.mfrow)
  if(is.null(dim(y))){
    coefs_bb <- object$mean.coefs.se[1:ncol(X)]
    lp_bb <- X%*%coefs_bb + offset
    pi_bb <- plogis(lp_bb)
    res <- rqr_function(X,y,coefs_bb,offset)
    qqnorm(res)
    qqline(res, col = 2)
    hist(res, xlab = "Residuals", main = "Histogram of residuals")
    plot(pi_bb, res, main = "Fitted Values Vs. Residuals", 
         xlab = "Predicted Values", ylab = "Residuals",...)
    plot(y,pi_bb, xlab = "Observed Values", ylab = "Predicted Values", 
         main = "Observed Vs. Fitted Values",...)
    mtext("Bayesian Bootstrap GDM diagnostics",outer = TRUE, cex=1.1,col="black",font=2,line=-1)
  } else {
  coefs_bb <- object$mean.coefs.se[1:ncol(X)]
  lp_bb <- X%*%coefs_bb + offset
  pi_bb <- plogis(lp_bb)
  res <- rqr_function(X,y,coefs_bb,offset)
  qqnorm(res)
  qqline(res, col = 2)
  hist(res, xlab = "Residuals", main = "Histogram of residuals")
  plot(pi_bb*y[,2], res, main = "Fitted Values Vs. Residuals", 
                                xlab = "Predicted Values", ylab = "Residuals",...)
  plot(y[,1],pi_bb*y[,2], xlab = "Observed Values", ylab = "Predicted Values", 
     main = "Observed Vs. Fitted Values",...)
  }
mtext("Bayesian Bootstrap GDM diagnostics",outer = TRUE, cex=1.1,col="black",font=2,line=-1)
}
 
