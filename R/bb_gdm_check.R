#' Function that plots diagnostic plots from bbgdm object
#' 
#' @param object bbgdm model output
#' @param rl.col colour of lines
#' @param family logical "binomial" or "beta_dist"
#' @return distances Five column matrix of sites pairs (cols 1:4) and distance as col 5.
#' @export
#' @examples
#' x <- matrix(rbinom(20*10,1,.6),20,10)# presence absence matrix
#' y <- simulate_covariates(x,2)
#' form <- ~ 1 + covar_1 + covar_2
#' test.bbgdm <- bbgdm(form,sp.dat=x, env.dat=y,family="binomial", 
#'                     dism_metric="number_non_shared",nboot=10,geo=FALSE, scale_covar=FALSE)
#' bbgdm.check(test.bbgdm)

bbgdm.check <- function (object, plots.mfrow = c(2, 2),...) 
{
  link <-object$link
  if(link=='negexp'){link.fun <- bbgdm::negexp()
  }else{ link.fun <- make.link(link=link)
  }
  family <-object$family
  X <- object$starting_gdm$X
  y <- object$starting_gdm$y
  offset <- object$offset
  par(mfrow=plots.mfrow)
  npreds <- object$nboots
   preds<- matrix(0L,npreds,nrow(y))
    for(i in 1:npreds)preds[i,] <- link.fun$linkinv(X%*%t(as.matrix(object$all.coefs.se[i,])))
    pi<-apply(preds,2,mean)
    a <- pbinom(y[,1]-1, y[,2], pi)#-1
    b <- pbinom(y[,1], y[,2], pi)
    u <- runif(n = length(y[,1]), min = a, max = b)
    u[u==1] <- u[u==1]-1e-5
    u[u==0] <- u[u==0]+1e-5
    res <- qnorm(u)
    qqnorm(res,ylab='Random Quantile Residuals',main = "")
    qqline(res, col = 'red')
    hist(res, xlab = "Random Quantile Residuals",main = "",...)
    plot(pi,res,cex=.5,pch=16, xlab="Predicted Dissimilarity",ylab="Random Quantile Residuals",...)
    plot(pi,y[,1]/y[,2],cex=.5,pch=16, xlab="Predicted Dissimilarity",ylab="Observed Dissimilarity",...)
    mtext("Bayesian Bootstrap GDM diagnostics",outer = TRUE, cex=1.1,col="black",font=2,line=-1)
}
 
