#' Negative exponentional link function.
#' A object of class \link[stats]{link-glm}, a list with components
#' @param linkfun	Link function function(mu)
#' @param linkinv	Inverse link function function(eta)
#' @param mu.eta	Derivative function(eta) dmu/deta
#' @param valideta	function(eta){ TRUE if eta is in the domain of linkinv }.
#' @param name a name to be used for the link
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