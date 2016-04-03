#' Function for plotting GDM w/ Bayesian Bootstrap.
#'
#' Plots a Generalised dissimilarity model with bayesian bootstrap.
#' @param bbgdm_object As derived from bbgdm function
#' @param plot.layout par settings for plot. default is c(2,2) (four plots).
#' @return a plots from bbgdm
#' @export
#' @examples
#' x <- matrix(rbinom(20*10,1,.6),20,10)# presence absence matrix
#' y <- simulate_covariates(x,2)
#' form <- ~ 1 + covar_1 + covar_2
#' test.bbgdm <- bbgdm(form,sp.dat=x, env.dat=y,family="binomial", dism_metric="bray_curtis",geo=FALSE,nboot=10, scale_covar=T)
#' plot(test.bbgdm)
#' 
## plotting
plot.bbgdm <- function(object,plot.layout = c(1,2),plot.colour='black',plot.linewidth=2,line.col='red'){
  if(nrow(object$X)>200) pred_sample <- 200
  else pred_sample <- nrow(object$X)
    par(mfrow=plot.layout)
    link.fun <- make.link("logit")
    X <- object$X
    Y <- object$starting_gdm$y
    offset <- object$offset
    bb.lp <- X%*%object$median.coefs.se + offset
    bb.lp.05 <- X%*%object$quantiles.coefs.se[1,] + offset
    bb.lp.95 <- X%*%object$quantiles.coefs.se[2,] + offset
    bb.pred <- link.fun$linkinv(bb.lp) 
    
    plot(bb.lp,Y[,1]/Y[,2],main='', xlab = "Linear Predictor", 
         ylab = "Observed Compositional Dissimilarity", type = "n",ylim=c(0,1))
    points(bb.lp,Y[,1]/Y[,2], pch = 20, cex = 0.25, 
           col = plot.colour)
    y.pred <- Y[sample(pred_sample),]
    y.pred <- y.pred[order(y.pred[,2], y.pred[,1]),]
    overlayX.bb <- seq(from = min(bb.lp), to = max(bb.lp), 
                    length = pred_sample)
    overlayY.bb <- link.fun$linkinv(overlayX.bb)
    lines(overlayX.bb, overlayY.bb,lwd = plot.linewidth,col=line.col)

    plot(bb.lp,Y[,1],xlab = "Linear Predictor",ylim=c(0,max(Y[,1])),#xlim=c(-5,5), 
           ylab = "Predicted non_shared species", type = "n")
    points(bb.lp,Y[,1], pch = 20, cex = 0.25, 
             col = plot.colour)
    y.pred <- Y[sample(pred_sample),]
    y.pred <- y.pred[order(y.pred[,2], y.pred[,1]),]
    overlayX.bb <- seq(from = min(bb.lp), to = max(bb.lp), 
                       length = pred_sample)
    overlayY.bb <- link.fun$linkinv(overlayX.bb)
    lines(overlayX.bb, overlayY.bb*y.pred[,2], lwd = plot.linewidth,col=line.col)
    lines(overlayX.org, overlayY.org*y.pred[,2], lwd = plot.linewidth,col="blue")
}
