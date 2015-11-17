#' Simulate bbgdms
#' 
#' @param X model matrix with intercept = 1
#' @param y no-non-shared-sp as two column response fro binomial model
#' @param Nsite No orginal sites
#' @param nboots number of bootstraps
#' @return a gdm.bb model object
#' @examples
#' library(bbgdm)
#' theta <- matrix(c(.5,2),1,2)
#' form <- y~1+x
#' nsp <- 250
#' Nsite <- 50
#' ndissim <- (Nsite*(Nsite-1))/2
#' dat <- data.frame(y=rep(1,ndissim),x=seq(0,4,length.out = ndissim))#imaginary difference in a covariate
#' sim1 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),.1),1,2),nsp)
#' y <- data.frame(SharedSpp_ij=sim1$nnss, MaxSpp_ij=sim1$max_sp_site)
#' test1 <- sim_bbgdm(sim1$X,y,Nsite,10)

sim_bbgdm <- function(X,y,offset=NULL,Nsite,nboots=10,control=logit_glm_control(...),...){
  x <- X[,-1,drop=FALSE]
  x <- spline.trans(x,spline_type = 'ispline')
  Xs <- cbind(1,x$spline)
  mods <- list()
  nreps <- nboots
  if(is.null(offset))offset <- rep(0,dim(X)[1])
  cat('start your engines\n')
  mod <- logit_glm_fit(Xs,y,offset=offset,control=control)
  for (ii in 1:nreps){
    w <- gtools::rdirichlet(Nsite, rep(1/Nsite,Nsite))
    wij <- w%*%t(w)
    wij <- wij[upper.tri(wij)]
    mods[[ii]] <- logit_glm_fit(Xs,y,wt=wij,offset=offset,control=control)
    if(ii %% 5 ==0 ) cat("Bayesian bootstrap ", ii, " iterations\n")
  }
  
  all.stats.ll.aic.bic.deviance <- plyr::ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained))
  median.ll.aic.bic.deviance <- apply( plyr::ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained)),2,median,na.rm=T)
  quantiles.ll.aic.bic.deviance <- apply( plyr::ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained)),2,function(x)quantile(x,c(.05,.95),na.rm=T))
  all.coefs.se <-  plyr::ldply(mods, function(x) c(x$coef))
  median.coefs.se <- apply(plyr::ldply(mods, function(x) c(x$coef)),2,median,na.rm=T)
  quantiles.coefs.se <- apply(plyr::ldply(mods, function(x) c(x$coef)),2,function(x)quantile(x,c(.05,.95),na.rm=T))
  #   vcov <- apply()
  gdm.bb.results <- list()
  gdm.bb.results$starting_gdm <- mod
  gdm.bb.results$bb_gdms <- mods
  gdm.bb.results$all.stats.ll.aic.bic.deviance <- all.stats.ll.aic.bic.deviance
  gdm.bb.results$median.ll.aic.bic.deviance <- median.ll.aic.bic.deviance
  gdm.bb.results$quantiles.ll.aic.bic.deviance <- quantiles.ll.aic.bic.deviance
  gdm.bb.results$all.coefs.se <- all.coefs.se
  gdm.bb.results$median.coefs.se <- median.coefs.se
  gdm.bb.results$quantiles.coefs.se <- quantiles.coefs.se
  gdm.bb.results$nboots <-nboots
#   gdm.bb.results$formula <- form
  gdm.bb.results$dism_metric <-  "number_shared"#dism_metric
#   gdm.bb.results$sp.dat <- sp.dat
  gdm.bb.results$env.dat <- X
  gdm.bb.results$X <- Xs
  gdm.bb.results$y <- y
  gdm.bb.results$offset <- offset
  gdm.bb.results$dissim_dat <- x$spline
  gdm.bb.results$dissim_dat_params <- x$spline.attr
  gdm.bb.results$family <- 'binomial'
  gdm.bb.results$geo <- FALSE
  gdm.bb.results$scale_covar <- FALSE#scale_covar
  class(gdm.bb.results) <- "gdm.bb"
  return(gdm.bb.results)
}