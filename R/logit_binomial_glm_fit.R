#' Function to perform postive binomial logistic regression.
#' 
#' @param X Model matirx of predictors, see \link[stats]{model.matrix}. 
#' @param y Model response, see \link[stats]{model.response}.
#' @param wt weights for model.
#' @param scale.covar logical If TRUE centre predictor variables.
#' @param control control option from optim see \link[BayesbootGDM]{logit_glm_control} or \link[stats]{optim}
#' @return fit fitted logistic binomial model as per optim methods
#' @export
 


logit_glm_fit <- function(X, y, wt=NULL,offset, control=logit_glm_control(...),...){
  my.fun <- function(x) {
    -LogLikFun(x, X, y, wt, offset)
  }
  my.grad <- function(x) {
    -GradFun(x, X, y, wt, offset)
  }
  if(is.null(wt)){
    if(!is.null(dim(y))) wt <- rep(1,nrow(y))
    else wt <- rep(1,length(y))
  }
  nbeta <- ncol(X)
  method <- control$method
  hessian <- control$hessian
  start <- control$start
  fsmaxit <- control$fsmaxit
  fstol <- control$fstol
  control$method <- control$hessian <- control$start <- control$fsmaxit <- control$fstol <- NULL
  if(is.null(start)){ 
    fm.b <- glm.fit(X,y,family=binomial("logit"))
    start <- c(fm.b$coefficients)
  }
  fit <- optim(par = start, fn = my.fun, gr =my.grad,  
                 method = method, hessian = hessian, control = control)
#   fit <- nlminb(start=start,objective=my.fun,gradient=my.grad)
  invisible(fit)
  fit$par <- c(fit$par[1],exp(fit$par[-1]))
  lp <- (X %*% fit$par + offset)
  p<-make.link(link = "logit")
  pi <- p$linkinv(lp)
  dev.weights <- rep(1,length(wt))
  var <- NULL
  if (hessian) {
    print("Calculating the variance of the estiamtes")
    var <- solve(fit$hessian)#(nH2(pt = esti$par, fun = my.fun))
    colnames(var) <- rownames(var) <- names(fit$par)
  }
  out <- list(coef = fit$par, logl = fit$value, 
              counts=fit$counts, fitted = pi, var=var,#residuals = res, 
              X=X, y=y,method = method, control=control)
  class(out) <- "logit_bin_GLM"
  return(out)
}
