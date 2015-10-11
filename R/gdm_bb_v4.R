#' Function to perform GDM w/ Bayesian Bootstrap.
#' 
#' Runs a Generalised dissimilarity model with bayesian bootstrap.
#' @param form formula for gdm.bb model 
#' @param family a description of the error distribution and link function to be used in the model. Currently "binomial" suppported. This can be a character string naming a family function, a family function or the result of a call to a family function.
#' @param dism_metric dissimilarity metric to calculate for model. "bray_curtis" or "number_shared" currently avaliable.
#' @param nboot number of Bayesian Bootstraps to run, this is used to estimate variance around GDM models. Default is 100 iterations.
#' @param sp.dat presence absence matrix, sp as columns sites as rows.
#' @param env.dat environmental or spatial covariates at each site.
#' @param spline_type type of spline to use in GDM model. Default is monotonic isplines. Options are: "ispline" or "bspline".
#' @return a gdm.bb model object
#' @export
#' @examples
#' sp.dat <- matrix(rbinom(200,1,.6),20,10)# presence absence matrix
#' env.dat <- simulate_covariates(x,2)
#' form <- ~ 1 + covar_1 + covar_2
#' test.gdm.bb <- gdm.bb(form,sp.dat, env.dat,family="binomial", dism_metric="number_shared",nboot=100, scale_covar=F)

gdm.bb <- function(form, sp.dat, env.dat, family="binomial", dism_metric="number_shared", nboot=100, 
                   spline_type="ispline",spline_df=2,spline_knots=1,scale_covar=FALSE,
                   geo=TRUE,geo.type='euclidean',coord.names=c("X","Y"),
                   lc_data=NULL,minr=0,maxr=NULL,control=logit_glm_control()){
  cat(family,"regression is on the way. \n")
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
  if(dism_metric=="number_shared") left <- "cbind(SharedSpp_ij,MaxSpp_ij)"
  if(dism_metric=="bray_curtis") left <- "dissimilarity"
  if(dism_metric=="simpsons") left <- "dissimilarity"
  if(geo) { form <- update.formula(form, ~ X + Y + .)
            if(!all(coord.names %in% colnames(env.dat)))
            {
              stop(cat("Coordinates are missing, add coordinate data as columns called 'X' & 'Y'\n"))
            }
  }          
  env.dat <- model.frame(form, as.data.frame(env.dat))
  mean.env.dat <- sd.env.dat <- NA
  if (scale_covar) {
    X.t <- env.dat
    mean.env.dat <- apply(X.t, 2, mean)
    sd.env.dat <- apply(X.t, 2, sd)
    env.dat <- scale(X.t)        
  }
  env.dat <- model.frame(as.data.frame(env.dat))
  offset <- model.offset(env.dat)
  if(!is.null(offset)) {
    offset_name <- colnames(env.dat)[length(colnames(env.dat))]
    env.dat <- env.dat[colnames(env.dat)[1:(length(colnames(env.dat))-1)]]
  }
  dissim_dat <- dissim_table(sp.dat,env.dat,dism_metric=dism_metric,spline_type=spline_type,spline_df=spline_df,spline_knots=spline_knots,
                             geo=geo,geo.type=geo.type,lc_data=lc_data,minr=minr,maxr=maxr)
  dissim_dat_table <- as.data.frame(dissim_dat$diff_table)
  dissim_dat_params <- dissim_dat$diff_table_params
  if(dism_metric=="number_shared") preds <- colnames(dissim_dat_table[,3:ncol(dissim_dat_table)])
  else preds <- colnames(dissim_dat_table[,2:ncol(dissim_dat_table)])
  if(!is.null(offset)){ form <- update.formula(form, paste0(left," ~ ",paste(preds, collapse=" + "),"+ offset(offset_ij)"))
                        offset_ij <- diff_table_cpp(as.matrix(offset))
                        dissim_dat_table <- as.data.frame(cbind(dissim_dat_table,offset_ij))
                        colnames(dissim_dat_table)[length(colnames(dissim_dat_table))] <- "offset_ij"
  } else { form <- update.formula(form, paste0(left," ~ ",paste(preds, collapse=" + ")))
  }
  temp <- model.frame(form, as.data.frame(dissim_dat_table))
  y <- model.response(temp)
  offset <- model.offset(temp)
  if (is.null(offset)) 
    offset <- rep(0, nrow(temp))
  if (!is.null(offset) && length(offset) != NROW(y)) {
    stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                  length(offset), NROW(y)), domain = NA)}
    X <- model.matrix(form, as.data.frame(dissim_dat_table))
    mod  <- logit_glm_fit(X,y,offset=offset,control=control)
    Nsite <- nrow(sp.dat)
    nreps <- nboot
    mods <- list()
    for (ii in 1:nreps){
      w <- gtools::rdirichlet(Nsite, rep(1/Nsite,Nsite))
      wij <- w%*%t(w)
      wij <- wij[upper.tri(wij)]
      mods[[ii]] <- logit_glm_fit(X,y,wt=wij,offset=offset,control=control)
      cat(ii,"\n")
      if(ii %% 20 ==0 ) cat("Bayesian bootstrap ", ii, " iterations\n")
    }
  #summary stats
  library(plyr)
  all.stats.ll.aic.bic.deviance <- ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained))
  median.ll.aic.bic.deviance <- apply( ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained)),2,median,na.rm=T)
  quantiles.ll.aic.bic.deviance <- apply( ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained)),2,function(x)quantile(x,c(.05,.95),na.rm=T))
  all.coefs.se <-  ldply(mods, function(x) c(x$coef))
  median.coefs.se <- apply(ldply(mods, function(x) c(x$coef)),2,median,na.rm=T)
  quantiles.coefs.se <- apply(ldply(mods, function(x) c(x$coef)),2,function(x)quantile(x,c(.05,.95),na.rm=T))
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
  gdm.bb.results$nboots <-nboot
  gdm.bb.results$formula <- form
  gdm.bb.results$dism_metric <- dism_metric
  gdm.bb.results$sp.dat <- sp.dat
  gdm.bb.results$env.dat <- env.dat
  gdm.bb.results$X <- X
  gdm.bb.results$y <- y
  gdm.bb.results$offset <- offset
  gdm.bb.results$dissim_dat <- dissim_dat_table
  gdm.bb.results$dissim_dat_params <- dissim_dat_params
  gdm.bb.results$family <- as.character(family)[1]
  gdm.bb.results$geo <- geo
  gdm.bb.results$scale_covar <- scale_covar
  if(scale_covar){
    gdm.bb.results$mean.env.dat <- mean.env.dat
    gdm.bb.results$sd.env.dat <- sd.env.dat
  }
  if(geo){
    gdm.bb.results$geo.type <- geo.type
    gdm.bb.results$lc_data=lc_data
    gdm.bb.results$minr=minr
    gdm.bb.results$maxr=maxr
  }
  class(gdm.bb.results) <- "gdm.bb"
  return(gdm.bb.results)
}