#' Negative exponential link function.
#'
#' A object of class \code{\link[stats]{glm}}, a list with components
#' @return negexp() link function
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
