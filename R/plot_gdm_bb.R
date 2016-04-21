#' Function for plotting GDM w/ Bayesian Bootstrap.
#'
#' Plots a Generalised dissimilarity model with bayesian bootstrap.
#' @param object As derived from bbgdm function
#' @param plot.layout par settings for plot. default is c(2,2) (four plots).
#' @param plot.colour colour of plot points
#' @param plot.linewidth line weight
#' @param line.col line colour
#' @param ... other plot arguments
#' @return a plots from bbgdm
#' @export
#' @examples
#' x <- matrix(rbinom(20*10,1,.6),20,10)# presence absence matrix
#' y <- simulate_covariates(x,2)
#' form <- ~ 1 + covar_1 + covar_2
#' test.bbgdm <- bbgdm(form,sp.dat=x, env.dat=y,family="binomial", dism_metric="number_non_shared",
#'                    nboot=10,geo=FALSE)
#' plot(test.bbgdm)
#'
## plotting
plot.bbgdm <- function(object,plot.layout = c(1,1),plot.colour='black',plot.linewidth=2,line.col='red',...){
  if(nrow(object$X)>200) pred_sample <- 200
  else pred_sample <- nrow(object$X)
    par(mfrow=plot.layout,oma = c(0, 0, 2, 0))
    link <-object$link
    if(link=='negexp') link.fun <- bbgdm::negexp()
    else link.fun <- make.link(link=link)
    X <- object$X
    Y <- object$starting_gdm$y
    offset <- object$offset
    bb.eta <- X%*%object$median.coefs.se + offset
    bb.eta.05 <- X%*%object$quantiles.coefs.se[1,] + offset
    bb.eta.95 <- X%*%object$quantiles.coefs.se[2,] + offset
    bb.pred <- link.fun$linkinv(bb.eta)
    plot(bb.eta,Y[,1]/Y[,2],main="Observed compositional dissimilarity vs. the linear predictor", xlab = "Linear Predictor",
         ylab = paste0("Observed ",object$dism_metric," dissimilarity"), type = "n",ylim=c(0,1),...)
    points(bb.eta,Y[,1]/Y[,2], pch = 20, cex = 0.25,
           col = plot.colour)
    y.pred <- Y[sample(pred_sample),]
    y.pred <- y.pred[order(y.pred[,2], y.pred[,1]),]
    overlayX.bb <- seq(from = min(bb.eta), to = max(bb.eta),
                    length = pred_sample)
    overlayY.bb <- link.fun$linkinv(overlayX.bb)
    lines(overlayX.bb, overlayY.bb,lwd = plot.linewidth,col=line.col)
}
