#' Function to perform GDM w/ Bayesian Bootstrap.
#' 
#' Runs a Generalised dissimilarity model with bayesian bootstrap.
#' @param form formula for bbgdm model 
#' @param family a description of the error distribution and link function to be used in the model. Currently "binomial" suppported. This can be a character string naming a family function, a family function or the result of a call to a family function.
#' @param dism_metric dissimilarity metric to calculate for model. "bray_curtis" or "number_non_shared" currently avaliable.
#' @param nboot number of Bayesian Bootstraps to run, this is used to estimate variance around GDM models. Default is 100 iterations.
#' @param sp.dat presence absence matrix, sp as columns sites as rows.
#' @param env.dat environmental or spatial covariates at each site.
#' @param spline_type type of spline to use in GDM model. Default is monotonic isplines. Options are: "ispline" or "bspline".
#' @param optim.meth "optim", "nlmnib" or "admb" for different optimization approaches. 
#' @return a bbgdm model object
#' @export
#' @examples
#' 
#' sp.dat <- matrix(rbinom(200,1,.6),20,10)# presence absence matrix
#' env.dat <- simulate_covariates(sp.dat,2)
#' form <- ~ 1 + covar_1 + covar_2
#' test.bbgdm <- bbgdm(form,sp.dat, env.dat,family="binomial",dism_metric="number_non_shared",nboot=10, scale_covar=F,geo=F,optim.meth='optim')

bbgdm <- function(form, sp.dat, env.dat, family="binomial", dism_metric="number_non_shared", nboot=100, 
                   spline_type="ispline",spline_df=2,spline_knots=1,scale_covar=FALSE,
                   geo=TRUE,geo.type='euclidean',coord.names=c("X","Y"),
                   lc_data=NULL,minr=0,maxr=NULL,
                   optim.meth="nlmnib", est.var=FALSE, trace=FALSE,prior=FALSE,
                   control=logit_glm_control()){
  cat(family,"regression is on the way. \n")
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
  if(dism_metric=="number_non_shared") left <- "cbind(nonsharedspp_ij,sumspp_ij)"
  if(dism_metric=="bray_curtis") left <- "cbind(dissimilarity,100)"
#   if(dism_metric=="simpsons") left <- "dissimilarity"
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
  if(dism_metric=="number_non_shared") preds <- colnames(dissim_dat_table[,3:ncol(dissim_dat_table)])
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
    mod  <- logit_glm_fit(X,y,offset=offset,optim.meth=optim.meth, est.var=TRUE, trace=trace,prior=prior,control=control)
    Nsite <- nrow(sp.dat)
    nreps <- nboot
    boot_print <- nboot/10
    mods <- list()
    for (ii in 1:nreps){
      w <- gtools::rdirichlet(Nsite, rep(1/Nsite,Nsite))
      wij <- w%*%t(w)
      wij <- wij[upper.tri(wij)]
      mods[[ii]] <- logit_glm_fit(X,y,wt=wij,offset=offset,optim.meth=optim.meth, est.var=est.var, trace=trace,prior=prior,control=control)
#       cat(ii,"\n")
      if(ii %% boot_print ==0 ) cat("Bayesian bootstrap ", ii, " iterations\n")
    }
  #summary stats
  all.stats.ll <- ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained))
  median.ll <- apply( ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained)),2,median,na.rm=T)
  quantiles.ll <- apply( ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained)),2,function(x)quantile(x,c(.05,.95),na.rm=T))
  all.coefs.se <-  ldply(mods, function(x) c(x$coef))
  median.coefs.se <- apply(ldply(mods, function(x) c(x$coef)),2,median,na.rm=T)
  quantiles.coefs.se <- apply(ldply(mods, function(x) c(x$coef)),2,function(x)quantile(x,c(.05,.95),na.rm=T))
#   vcov <- apply()
  
  bbgdm.results <- list()
  bbgdm.results$starting_gdm <- mod
  bbgdm.results$bb_gdms <- mods
  bbgdm.results$all.stats.ll <- all.stats.ll
  bbgdm.results$median.ll <- median.ll
  bbgdm.results$quantiles.ll <- quantiles.ll
  bbgdm.results$all.coefs.se <- all.coefs.se
  bbgdm.results$median.coefs.se <- median.coefs.se
  bbgdm.results$quantiles.coefs.se <- quantiles.coefs.se
  bbgdm.results$nboots <-nboot
  bbgdm.results$formula <- form
  bbgdm.results$dism_metric <- dism_metric
  bbgdm.results$sp.dat <- sp.dat
  bbgdm.results$env.dat <- env.dat
  bbgdm.results$X <- X
  bbgdm.results$y <- y
  bbgdm.results$offset <- offset
  bbgdm.results$dissim_dat <- dissim_dat_table
  bbgdm.results$dissim_dat_params <- dissim_dat_params
  bbgdm.results$family <- as.character(family)[1]
  bbgdm.results$geo <- geo
  bbgdm.results$scale_covar <- scale_covar
  if(scale_covar){
    bbgdm.results$mean.env.dat <- mean.env.dat
    bbgdm.results$sd.env.dat <- sd.env.dat
  }
  if(geo){
    bbgdm.results$geo.type <- geo.type
    bbgdm.results$lc_data=lc_data
    bbgdm.results$minr=minr
    bbgdm.results$maxr=maxr
  }
  class(bbgdm.results) <- "bbgdm"
  return(bbgdm.results)
}