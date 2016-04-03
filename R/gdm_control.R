#' Control Parameters for Binomial Logistic Regression

#' Various parameters that control fitting of binomial regression models using BBgdm.

#'@param method characters string specifying the method argument passed to optim.
#'@param maxit  integer specifying the maxit argument (maximal number of iterations) passed to optim.
#'@param trace	logical or integer controlling whether tracing information on the progress of the optimization should be produced (passed to optim).
#'@param hessian	logical. Should the numerical Hessian matrix from the optim output be used for estimation of the covariance matrix? By default the analytical solution is employed. For details see below.
#'@param start	an optional vector with starting values for all parameters.
#'@param fsmaxit	integer specifying maximal number of additional (quasi) Fisher scoring iterations. For details see below.
#'@param fstol	numeric tolerance for convergence in (quasi) Fisher scoring. For details see \link[stats]{optim}.
#'@param ...	arguments passed to optim.
#'@export

gdm_control <- function (method = "BFGS", maxit = 1000, hessian = FALSE, 
                               trace = FALSE, start = NULL, fsmaxit = 20, fstol = 1e-05, 
                               ...) 
{
  rval <- list(method = method, maxit = maxit, hessian = hessian, trace = trace, start = start, fsmaxit = fsmaxit, 
               fstol = fstol)
  rval <- c(rval, list(...))
  if (!is.null(rval$fnscale)) 
    warning("fnscale must not be modified")
  rval$fnscale <- 1
  if (is.null(rval$reltol)) 
    rval$reltol <- 1e-05
  rval
}
