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
#   var <- NULL
#   if (est.var) {
#     print("Calculating the variance of the estiamtes")
#     var <- solve(fit$hessian)#(nH2(pt = esti$par, fun = my.fun))
#     colnames(var) <- rownames(var) <- names(esti$par)
#   }
  if(!is.null(dim(y))){
    a <- pbinom((y[,1]-1), y[,2], pi)#-1
    b <- pbinom(y[,1], y[,2], pi)
    u <- runif(n = length(y[,1]), min = a, max = b)
    res <- qnorm(u)
    res[res==Inf]<-max(res[res!=Inf])
    res[res==-Inf]<-min(res[res!=-Inf])
    wtdmu <- sum(dev.weights * (y[,1]/y[,2]))/sum(dev.weights)
    d <- sum(binomial()$dev.resids(y[,1]/y[,2],pi, dev.weights),na.rm=T)
    dn <- sum(binomial()$dev.resids(y[,1]/y[,2],wtdmu, dev.weights),na.rm=T)
    dev_exp <- (dn - d) /dn
  } else {
    n <- rep(1, length(y))
    y <- n * y
    a <- pbinom(y-1, n, pi)
    b <- pbinom(y, n, pi)
    u <- runif(n = length(y), min = a, max = b)
    res <- qnorm(u)
    res[res==Inf]<-max(res[res!=Inf])
    res[res==-Inf]<-min(res[res!=-Inf])
    wtdmu <- sum(dev.weights * y)/sum(dev.weights)
    d <- sum(binomial()$dev.resids(y,pi, dev.weights),na.rm=T)
    dn <- sum(binomial()$dev.resids(y,wtdmu, dev.weights),na.rm=T)
    dev_exp <- (dn - d) /dn
  }
  AIC <- 2 * fit$value + 2 * length(fit$par)
  BIC <- 2 * fit$value + log(nrow(X)) * length(fit$par)
  out <- list(coef = fit$par, logl = fit$value, AIC = AIC, BIC=BIC,
              gdm.deviance=d,null.deviance=dn,deviance.explained=dev_exp,
              counts=fit$counts, fitted = pi, residuals = res, 
              X=X, y=y,method = method, control=control)
  class(out) <- "logit_bin_GLM"
  return(out)
}
