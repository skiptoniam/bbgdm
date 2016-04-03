#' Function to perform postive binomial logistic regression.
#' 
#' @param X Model matirx of predictors, see \link[stats]{model.matrix}. 
#' @param y Model response, see \link[stats]{model.response}.
#' @param wt weights for model.
#' @param scale.covar logical If TRUE centre predictor variables.
#' @param control control option from optim see \link[BayesbootGDM]{gdm_fit_control} or \link[stats]{optim}
#' @return fit fitted logistic binomial model as per optim methods
#' @export
 


gdm_fit <- function(X, y, wt=NULL,offset,optim.meth="optim", est.var=TRUE, trace=FALSE,prior=FALSE,
                          link, control=gdm_fit_control(...),...){
  
  my.fun <- function(x) {
    -LogLikFun(x, X, y, wt, offset,link=link)
  }
  my.grad <- function(x) {
    -GradFun(x, X, y, wt, offset,link=link)
  }
  if(is.null(wt)){
    if(!is.null(dim(y))) wt <- rep(1,nrow(y))
    else wt <- rep(1,length(y))
  }
  if(optim.meth=="optim"){
    if(trace)control$trace <- 1
  method <- control$method
  hessian <- control$hessian
  init.par <- control$start
  fsmaxit <- control$fsmaxit
  fstol <- control$fstol
  control$method <- control$hessian <- control$start <- control$fsmaxit <- control$fstol <- NULL
  if(is.null(init.par)){ 
    fm.b <- glm.fit(X,y,family=binomial("logit"))
    init.par <- c(fm.b$coefficients)
    if(prior)init.par<-rbeta(length(init.par),2,2)
  }
  if(est.var)hessian <- TRUE
  fit <- optim(par = init.par, fn = my.fun, gr =my.grad,  
                 method = method, hessian = hessian, control = control)
  }
  if (optim.meth=="nlmnib"){
    init.par <- control$start
    if(is.null(init.par)){
      fm.b <- glm.fit(X,y,family=binomial(link))
      init.par <- c(fm.b$coefficients)
      if(prior)init.par<-rbeta(length(init.par),2,2)
    }
    fit <- nlminb(start=init.par,objective=my.fun,gradient=my.grad,control = list(trace=trace))
    fit$value <- fit$objective
    fit$counts <- fit$evaluations
  }
  if (optim.meth=='admb'){
    init.par <- control$start
    if(is.null(init.par)){
      init.par <- rep(0,ncol(X))
      if(prior)init.par<-rbeta(length(init.par),2,2)
    }
    dyn.load(dynlib("logit_reg"))
    obj <- MakeADFun(
      data = list(x = X[,-1], y = y, w = wt,offset=offset), 
      parameters = list(a =init.par[1], b = init.par[-1]),
      DLL = "logit_reg",hessian=TRUE,silent=TRUE)
      fit <- suppressWarnings(nlminb(obj$par,obj$fn,obj$gr,control =list(trace=trace)))
      fit$value <- fit$objective
      fit$counts <- fit$evaluations
  }
  invisible(fit)
  fit$par <- c(fit$par[1],exp(fit$par[-1]))
  names(fit$par) <- colnames(X)
  var <- NULL
  if (est.var) {
    # cat("Calculating the variance of the estimates.\n")
    if(optim.meth=='optim') var <- solve(fit$hessian)
    if(optim.meth=='nlmnib') var <- solve(numDeriv::hessian(my.fun,init.par))
    if(optim.meth=='admb') var <- solve(numDeriv::hessian(obj$fn,obj$par))
    colnames(var) <- rownames(var) <- names(fit$par)
  }
  out <- list(coef = fit$par, logl = fit$value, 
              counts=fit$counts, fitted = pi, var=var,
              X=X, y=y,control=control)
  class(out) <- "gdm"
  return(out)
}
