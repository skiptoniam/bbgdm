#' Simulate Covariates for bbGDM
#' 
#' Creates N covarites for testing bbGDM.
#' @param x Presence Absence Matrix.
#' @param N number of covariates to simulate
#' @return beta N covariates for number of sites of x.
#' @export
#' @examples
#' x <- matrix(sample(0:1,16, replace=T),4,4) #toy presence absence matrix
#' covr <- simulate_covariates(x,6)

simulate_covariates <- function(x,N){
  Nsite <- nrow(x)
  Ncov <- N
  beta <- matrix(rnorm(Nsite*Ncov, 0, 2),nrow=Nsite, ncol=Ncov)
  colnames(beta) <- paste0("covar_",1:Ncov)
  return(as.data.frame(beta))
}