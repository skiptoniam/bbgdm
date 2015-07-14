#' Function to create a invexp link function for model.
#' For use in bbGDM
#' @export
#' @examples
#' inex <- invexp()
#' inex$linkfun(inex$linkinv(.5))  ## check invertibility
#' library("numDeriv")
#' all.equal(grad(inex$linkinv,.2),inex$mu.eta(.2))  ## check derivative
# 
#' set.seed(101)
#' n <- 100                       
#' y <- sort(runif(n,0.1,1))
#' x <- sort(runif(n,1,10))
#' plot(-1,ylab='dissimilarity',xlab='distance',ylim=c(0,1),xlim=c(0,100),xaxt = "n")
#' lines(sort(predict(glm(y~x,family=binomial(link=inex)),type="response")))

invexp <- function() {
  ## link
  linkfun <- function(y) 1-exp(-1*(y))
  ## inverse link
  linkinv <- function(eta) -log(1-eta)
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) { 1/(1-eta) }
  valideta <- function(eta) TRUE
  link <- "invexp"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}

inex <- invexp()
inex$linkfun(inex$linkinv(.5))  ## check invertibility
library("numDeriv")
all.equal(grad(inex$linkinv,.2),inex$mu.eta(.2))  ## check derivative
