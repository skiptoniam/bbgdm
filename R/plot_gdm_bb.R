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
#' test.bbgdm <- bbgdm(form,sp.dat=x, env.dat=y,family="binomial", dism_metric="number_shared",nboot=10, scale_covar=T)
#' plot(test.bbgdm)
#' 
## plotting
plot.bbgdm <- function(object,plot.layout = c(2,2),plot.colour='black',plot.linewidth=2,line.col='red'){
  if(nrow(object$X)>200) pred_sample <- 200
  else pred_sample <- nrow(object$X)
    par(mfrow=plot.layout)
    link.fun <- make.link("logit")
    X <- object$X
    Y <- object$starting_gdm$y
    offset <- object$offset
    gdm.lp <-  X%*%object$starting_gdm$coef + offset
    gdm.pred <- link.fun$linkinv(gdm.lp)
    plot(gdm.lp,Y[,1] ,main='starting GDM', xlab = "Predicted Ecological Distance", 
         ylab = "Observed Non-Shared Species", type = "n")
    points(gdm.lp,Y[,1] , pch = 20, cex = 0.25,col = plot.colour)
    y.pred <- Y[sample(pred_sample),]
    y.pred <- y.pred[order(y.pred[,2], y.pred[,1]),]
    overlayX.org <- seq(from = min(gdm.lp), to = max(gdm.lp), 
                    length = pred_sample)
    overlayY.org <- link.fun$linkinv(overlayX.org)
    lines(overlayX.org, overlayY.org*y.pred[,2], lwd = plot.linewidth,col=line.col)
    
    plot(gdm.lp,Y[,1]/Y[,2] ,main='starting GDM', xlab = "Predicted Ecological Distance", 
       ylab = "Observed Compositional Dissimilarity", type = "n", ylim=c(0,1))
    points(gdm.lp,Y[,1]/Y[,2] , pch = 20, cex = 0.25,col = plot.colour)
    overlayX.org <- seq(from = min(gdm.lp), to = max(gdm.lp), 
                        length = pred_sample)
    overlayY.org <- link.fun$linkinv(overlayX.org)
    lines(overlayX.org, overlayY.org, lwd = plot.linewidth,col=line.col)    

    bb.lp <- X%*%object$median.coefs.se + offset
    bb.lp.05 <- X%*%object$quantiles.coefs.se[1,] + offset
    bb.lp.95 <- X%*%object$quantiles.coefs.se[2,] + offset
    bb.pred <- link.fun$linkinv(bb.lp) #inverse_logit 1/(1+exp(-1*X%*%object$stats.median[1:length(coef(object$starting_gdm))]))
    
    plot(bb.lp,Y[,1]/Y[,2],main='bb.gdm estimate', xlab = "Linear Predictor", 
         ylab = "Observed Compositional Dissimilarity", type = "n",ylim=c(0,1))
    points(bb.lp,Y[,1]/Y[,2], pch = 20, cex = 0.25, 
           col = plot.colour)
    y.pred <- Y[sample(pred_sample),]
    y.pred <- y.pred[order(y.pred[,2], y.pred[,1]),]
    overlayX.bb <- seq(from = min(bb.lp), to = max(bb.lp), 
                    length = pred_sample)
    overlayY.bb <- link.fun$linkinv(overlayX.bb)
#     lines(overlayX, overlayY*y.pred[,2], lwd = plot.linewidth)
    lines(overlayX.bb, overlayY.bb,lwd = plot.linewidth,col=line.col)
#     lines(overlayX.org, overlayY.org,lwd = plot.linewidth,col="blue")

    plot(bb.lp,Y[,1],xlab = "Linear Predictor",ylim=c(0,max(Y[,1])),#xlim=c(-5,5), 
           ylab = "Predicted non-shared species", type = "n")
    points(bb.lp,Y[,1], pch = 20, cex = 0.25, 
             col = plot.colour)
    y.pred <- Y[sample(pred_sample),]
    y.pred <- y.pred[order(y.pred[,2], y.pred[,1]),]
    overlayX.bb <- seq(from = min(bb.lp), to = max(bb.lp), 
                       length = pred_sample)
    overlayY.bb <- link.fun$linkinv(overlayX.bb)
    lines(overlayX.bb, overlayY.bb*y.pred[,2], lwd = plot.linewidth,col=line.col)
    lines(overlayX.org, overlayY.org*y.pred[,2], lwd = plot.linewidth,col="blue")
    legend("bottomright", c("BBGDM mean","Single GDM mean"), pch = 16,col = c("red", "blue"))    
    
    y.est.runs <- matrix(NA,nrow(X),object$nboots)
    lp.runs <- matrix(NA,nrow(X),object$nboots)
    for(ii in 1:object$nboots){
      lp <- X%*%object$bb_gdms[[ii]]$coef
      link.fun$linkinv(lp)
      y.est.runs[,ii] <- link.fun$linkinv(lp)
      lp.runs[,ii]<-lp
    }
    bb.lp <- X%*%object$median.coefs.se + offset
    quantiles.bb.lp.05 <- X%*%object$quantiles.coefs.se[1,] + offset
    quantiles.bb.lp.95 <- X%*%object$quantiles.coefs.se[2,] + offset
    
    yhat_mean <- link.fun$linkinv(bb.lp) #inverse_logit 1/(1+exp(-1*X%*%object$stats.median[1:length(coef(object$starting_gdm))]))
    yhat_q05 <- link.fun$linkinv(quantiles.bb.lp.05)
    yhat_q95 <- link.fun$linkinv(quantiles.bb.lp.95)
    
#     plot(-1,ylab='Predicted Compositional Dissimilarity',xlab='Observed Compositional Dissimilarity',ylim=range(Y[,1]),xlim=range(c(yhat_mean*Y[,2],yhat_q05*Y[,2],yhat_q95*Y[,2])),type = "n")
#     points(yhat_q05*Y[,2], Y[,1], pch=16, cex=.5,col='blue')
#     points(yhat_q95*Y[,2], Y[,1], pch=16, cex=.5,,col='red')
#     points(yhat_mean*Y[,2], Y[,1], pch=16, cex=.5,col = "orange")
#     legend("bottomright", c("5% quantile","mean","95% quantile"), pch = 16,col = c("blue", "orange", "red"))    

}
