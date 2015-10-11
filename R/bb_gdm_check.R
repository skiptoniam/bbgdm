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
  p <- make.link(link = "logit")
  family <-object$family
  X <- object$starting_gdm$X
  y <- object$starting_gdm$y
  offset <- object$offset
  coefs<-object$starting_gdm$coef
  par(mfrow=plots.mfrow)
  res <- rqr_function(X,y,coefs,offset)
  qqnorm(res)
  qqline(res, col = 'red')
  hist(res, xlab = "Residuals", main = "Histogram of residuals",...)
  lp<- X%*%coefs + offset
  pi <- p$linkinv(lp)
  plot(pi, res, xlab="Predicted Dissimilarity",ylab="Random Quantile Residuals",cex=.5,pch=16,...)
  plot(pi,y[,1]/y[,2], xlab="Predicted Dissimilarity",ylab="Observed Dissimilarity",cex=.5,pch=16,...)
  mtext("Original GDM diagnostics",outer = TRUE, cex=1.1,col="black",font=2,line=-1)
  ## Now for bayesian boot strap estimates
  par(mfrow=plots.mfrow)
  npreds <- object$nboots
   preds<- matrix(0L,npreds,nrow(y))
    for(i in 1:npreds)preds[i,] <- p$linkinv(X%*%t(as.matrix(object$all.coefs.se[i,])))
    pi<-apply(preds,2,mean)
    a <- pbinom(y[,1]-1, y[,2], pi)#-1
    b <- pbinom(y[,1], y[,2], pi)
    u <- runif(n = length(y[,1]), min = a, max = b)
    res <- qnorm(u)
    qqnorm(res,ylab='Random Quantile Residuals',main = "")
    qqline(res, col = 'red')
    hist(res, xlab = "Random Quantile Residuals",main = "",...)
    plot(pi,res,cex=.5,pch=16, xlab="Predicted Dissimilarity",ylab="Random Quantile Residuals",...)
    plot(pi,y[,1]/y[,2],cex=.5,pch=16, xlab="Predicted Dissimilarity",ylab="Observed Dissimilarity",...)
    mtext("Bayesian Bootstrap GDM diagnostics",outer = TRUE, cex=1.1,col="black",font=2,line=-1)

}
 
