#' @title bbgdm.fit objects
#' @rdname bbgdm.fit
#' @name bbgdm.fit
#' @description function for fitting a \code{bbgdm} model. This is essential a single gdm
#' and builds the core of \code{bbgdm}. Like \code{\link[stats]{glm.fit}} a number of options
#' can be defined. Like: optimisation method \code{\link[stats]{optim}} or \code{\link[stats]{nlminb}},
#' link function \code{\link[stats]{make.link}}.
#' @param \dots for \code{bbgdm.fit()} or more \code{bbgdm.control()}.
#' @param X Model matirx of predictors, see \code{\link[stats]{model.matrix}}.
#' @param y Model response, see \code{\link[stats]{model.response}}.
#' @param wt weights for model.
#' @param link character link functions. default is 'logit', can call other
#' link functions \code{\link[stats]{make.link}} or a 'negexp' custom link function (e.g., Ferrier etal., 2007).
#' @param optim.meth optimisation method options avaliable are 'optim' and 'nlmnib',
#' calls either method \code{\link[stats]{optim}} or \code{\link[stats]{nlminb}}.
#' @param est.var logical if true estimated parameter variance using optimiser.
#' @param trace options looks at \code{\link[stats]{optim}} for details.
#' @param prior numeric vector of starting values for intercept and splines.
#' @param control control option from optim see \code{\link[bbgdm]{bbgdm.control}} or \code{\link[stats]{optim}}.
#' @return fit single gdm
#' @export
#' @author Skipton Woolley

bbgdm.fit <- function(X, y, wt=NULL, link, optim.meth="optim", est.var=TRUE, trace=FALSE,prior=FALSE,
                    control=bbgdm.control()){

  loglike <- function(x) {
    -loglikelihood(x, X, y, wt, link)
  }
  grad <- function(x) {
    -gradient(x, X, y, wt, link)
  }
  if(is.null(wt)){
    if(!is.null(dim(y))) wt <- rep(1,nrow(y))
    else wt <- rep(1,length(y))
  }
  if(length(wt)!=nrow(X))stop('Weights are not the same size as model matrix')
  if(optim.meth=="optim"){
    if(trace)control$trace <- 1
    method <- control$method
    hessian <- control$hessian
    init.par <- control$start
    fsmaxit <- control$fsmaxit
    fstol <- control$fstol
    control$method <- control$hessian <- control$start <- control$fsmaxit <- control$fstol <- NULL
    if(is.null(init.par)){
      if(link=='negexp')fm.b <- glm.fit(X,y,family=binomial(link=negexp()))
      else fm.b <- glm.fit(X,y,family=binomial(link=link))
      init.par <- c(fm.b$coefficients)
      if(prior)init.par<-rbeta(length(init.par),2,2)
    }
    if(est.var)hessian <- TRUE
    fit <- optim(par = init.par, fn = loglike, gr =grad,
                 method = method, hessian = hessian, control = control)
  }
  if (optim.meth=="nlmnib"){
    init.par <- control$start
    if(is.null(init.par)){
      if(link=='negexp')fm.b <- glm.fit(X,y,family=binomial(link=negexp()))
      else fm.b <- glm.fit(X,y,family=binomial(link=link))
      init.par <- c(fm.b$coefficients)
      if(prior)init.par<-rbeta(length(init.par),2,2)
    }
    fit <- nlminb(start=init.par,objective=loglike,gradient=grad,control = list(trace=trace))
    fit$value <- fit$objective
    fit$counts <- fit$evaluations
  }
  invisible(fit)
  fit$par <- c(fit$par[1],exp(fit$par[-1]))
  names(fit$par) <- colnames(X)
  var <- NULL
  if (est.var) {
    if(optim.meth=='optim') var <- solve(fit$hessian)
    if(optim.meth=='nlmnib') var <- solve(numDeriv::hessian(loglike,init.par))
    colnames(var) <- rownames(var) <- names(fit$par)
  }
  out <- list(coef = fit$par, logl = fit$value,
              counts=fit$counts, fitted = pi, var=var,
              X=X, y=y,control=control)
  class(out) <- "gdm"
  return(out)
}

#'@rdname bbgdm.fit
#'@name bbgdm.control
#'@param method characters string specifying the method argument passed to optim.
#'@param maxit  integer specifying the maxit argument (maximal number of iterations) passed to optim.
#'@param hessian	logical. Should the numerical Hessian matrix from the optim output be used for estimation of the covariance matrix? By default the analytical solution is employed. For details see below.
#'@param start	an optional vector with starting values for all parameters.
#'@param fsmaxit	integer specifying maximal number of additional (quasi) Fisher scoring iterations. For details see below.
#'@param fstol	numeric tolerance for convergence in (quasi) Fisher scoring. For details see \code{\link[stats]{optim}}.
#'@param cores the number of cores to use in fitting.
#'@export

bbgdm.control <- function (method = "BFGS", maxit = 1000, hessian = FALSE,
                         trace = FALSE, start = NULL, fsmaxit = 20, fstol = 1e-05,cores=1,
                         ...)
{
  rval <- list(method = method, maxit = maxit, hessian = hessian, trace = trace, start = start, fsmaxit = fsmaxit, fstol = fstol,cores=cores)
  rval <- c(rval, list(...))
  if (!is.null(rval$fnscale))
    warning("fnscale must not be modified")
  rval$fnscale <- 1
  if (is.null(rval$reltol))
    rval$reltol <- 1e-05
  rval
}

loglikelihood <- function(params,X,y,wt,link){
  lp <- X %*% c(params[1],exp(params[-1]))
  if(link=='negexp') p <- bbgdm::negexp()
  else  p <- make.link(link = link)
  lb <- dbinom(y[,1], y[,2], p$linkinv(lp),log = TRUE)*wt
  ll.contr<-sum(lb)
  return(ll.contr)
}

gradient <- function(params,X,y,wt,link){
  eta <- X %*% c(params[1],exp(params[-1]))
  np<-length(params)
  ns<-nrow(X)
  deri<- matrix(NA,ns,np)
  if(link=='negexp')  p<- bbgdm::negexp()
  else p <- make.link(link = link)
  for(i in 1:ns){
    mu <- p$linkinv(eta[i])
    dldm <- ((y[i,1]/y[i,2])/mu) - ((1-(y[i,1]/y[i,2]))/(1-mu))
    dmde <- p$mu.eta(eta[i])*y[i,2]
    dedb <- X[i,]
    dbdg <- c(1,exp(params[-1]))
    deri[i,] <- wt[i] * dldm * dmde * dedb * dbdg #chain-rule baby!
  }
  sum_deri <- apply(deri,2,sum)
  return(sum_deri)
}

#' @rdname bbgdm.fit
#' @name negexp
#' @export

negexp<- function()
{
  linkfun <- function(mu) -log(1-mu)
  linkinv <- function(eta) 1-exp(-eta)
  mu.eta <- function(eta) exp(-eta)
  valideta <- function(eta) all(is.finite(eta))
  link <- paste0("negexp")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link),
            class = "link-glm")
}
