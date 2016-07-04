#' @title diagnostics objects
#' @rdname diagnostics
#' @name diagnostics
#' @param object bbgdm model output
#' @export
#' @examples
#' x <- matrix(rbinom(20*10,1,.6),20,10)# presence absence matrix
#' y <- simulate_covariates(x,2)
#' form <- ~ 1 + covar_1 + covar_2
#' test.bbgdm <- bbgdm(form,sp.dat=x, env.dat=y,family="binomial",
#'                     dism_metric="number_non_shared",nboot=10,geo=FALSE)
#' resids <- diagnostics(test.bbgdm)
#'

diagnostics <- function (object)
{
  link <-object$link
  if(link=='negexp'){link.fun <- bbgdm::negexp()
  }else{ link.fun <- make.link(link=link)
  }
  family <-object$family
  X <- object$starting_gdm$X
  y <- object$starting_gdm$y
  offset <- object$offset
  npreds <- object$nboots
  preds<- matrix(0L,npreds,nrow(y))
  for(i in 1:npreds)preds[i,] <- link.fun$linkinv(X%*%t(as.matrix(object$all.coefs.se[i,])))
  pi<-apply(preds,2,mean)
  a <- pbinom(y[,1]-1, y[,2], pi)#-1
  b <- pbinom(y[,1], y[,2], pi)
  u <- runif(n = length(y[,1]), min = a, max = b)
  u[u==1] <- u[u==1]-1e-5
  u[u==0] <- u[u==0]+1e-5
  res <- qnorm(u)
  structure(list(res=res,pi=pi,y=y),class='diagnostics')
}

#' @rdname diagnostics
#' @name plot.diagnostics
#' @export
#' @examples
#' #plot residuals
#' par(mfrow=c(2,2))
#' plot(resids)
#'

plot.diagnostics<- function(object,...){
  qqnorm(object$res,ylab='Random Quantile Residuals',main = "")
  qqline(object$res, col = 'red')
  hist(object$res, xlab = "Random Quantile Residuals",main = "",...)
  plot(object$res,object$pi,xlab="Predicted Dissimilarity",ylab="Random Quantile Residuals",...)
  plot(object$pi,object$y[,1]/object$y[,2],xlab="Predicted Dissimilarity",ylab="Observed Dissimilarity",...)
}

#' @rdname diagnostics
#' @name bbgdm.wald.test
#' @param object Returned model from \link[bbgdm]{bbgdm}.
#' @param H0 A numeric value giving the null hypothesis for the test. Generally zero.
#' @param gdm Logic if true calculates the Wald-test using variance-covariance matrix derived from the hessian matrix. Note: These estimates are probably wrong due to hessian matrix being calulated with respect to the likelihoods.
#' @export
#' @examples
#' #wald test on bbdgm
#' bbgdm.wald.test(test.bbgdm)
#'

bbgdm.wald.test <- function(object,H0=0,gdm=FALSE){
  # H0: Hypothesis test = 0
  # IM: Identiy Matrix
  # beta: parameter estimates from model
  # vcov: Variance-covariance matrix estimated from BB

  #for gdm
  if(gdm){
    vcov <-object$starting_gdm$var
    beta <- object$starting_gdm$coef
    #Intercept
    intercept_IM <- matrix(c(1,rep(0,length(beta)-1)),nrow=1)
    wd_inter <- t(intercept_IM%*%beta-H0) %*% solve(intercept_IM%*%vcov%*%t(intercept_IM))%*%(intercept_IM%*%beta-H0)
    pval_i = 1-pchisq(wd_inter,1)

    #Splines
    splineLength <- sapply(object$dissim_dat_params, `[[`, "dim")[2,]
    val1 <- seq(2,length(beta),splineLength[1])
    val2 <- seq(1+splineLength[1],length(beta),splineLength[1])
    wd_vals_gdm <- matrix(NA,length(val1)+1,3)
    wd_vals_gdm[1,]<- c(wd_inter,1,pval_i)
    for(i in 1:length(splineLength)){
      w <- splineLength[1]
      L <- matrix(rep(0, length(beta) * w), ncol = length(beta))
      Terms <- seq(val1[i],val2[i],1)
      for (ii in 1:w) L[ii, Terms[ii]] <- 1
      vcov1 <- try(solve(L %*% vcov %*% t(L)),silent = TRUE)
      if ((class(vcov1)=="try-error")||(class(vcov1)=="try-error"))
      {
        wd <- 0
      }else{
        wd <- t(L %*% beta - H0) %*% vcov1 %*% (L %*% beta - H0)
      }
      pv <- 1 - pchisq(wd, df = w)
      wd_vals_gdm[1+i,]<- c(wd,w,pv)
    }
    colnames(wd_vals_gdm) <- c("gdm_W","gdm_df","gdm_p-value")
    if(object$geo){ rownames(wd_vals_gdm)<-c('intercept','geo',names(object$env.dat)[-c(1:2)])
    } else { rownames(wd_vals_gdm)<-c('intercept',names(object$env.dat))
    }

  }

  A <- object$all.coefs.se #matrix of B bootstrap coeficient estimates.
  esti.var <- var(A) #make sure this is a matrix
  beta <- object$median.coefs.se #medians of coef estimates

  #Intercept
  intercept_IM <- matrix(c(1,rep(0,length(beta)-1)),nrow=1)
  wd_inter <- t(intercept_IM%*%beta-H0) %*% solve(intercept_IM%*%esti.var%*%t(intercept_IM))%*%(intercept_IM%*%beta-H0)
  pval_i = 1-pchisq(wd_inter,1)

  #Splines
  splineLength <- sapply(object$dissim_dat_params, `[[`, "dim")[2,]
  val1 <- seq(2,length(beta),splineLength[1])
  val2 <- seq(1+splineLength[1],length(beta),splineLength[1])
  wd_vals <- matrix(NA,length(val1)+1,3)
  wd_vals[1,]<- c(wd_inter,1,pval_i)
  for(i in 1:length(splineLength)){
    w <- splineLength[1]
    L <- matrix(rep(0, length(beta) * w), ncol = length(beta))
    Terms <- seq(val1[i],val2[i],1)
    for (ii in 1:w) L[ii, Terms[ii]] <- 1
    vcov2 <- try(solve(L %*% esti.var %*% t(L)),silent = TRUE)
    if ((class(vcov2)=="try-error")||(class(vcov2)=="try-error"))
    {
      wd <- 0
    }else{
      wd <- t(L %*% beta - H0) %*% vcov2 %*% (L %*% beta - H0)
    }
    pv <- 1 - pchisq(wd, df = w)
    wd_vals[1+i,]<- c(wd,w,pv)
  }
  colnames(wd_vals) <- c("bbgdm_W","bbgdm_df","bbgdm_p-value")
  if(object$geo){ rownames(wd_vals)<-c('intercept','geo',names(object$env.dat)[-c(1:2)])
  } else { rownames(wd_vals)<-c('intercept',names(object$env.dat))
  }
  if(gdm) return(cbind(wd_vals_gdm,wd_vals))
  else return(wd_vals)
}

#' @rdname diagnostics
#' @name ggplot.bbgdm
#' @param pars name of parameters to plot.
#' @param pars_labels vector of parameter labels.
#' @param zero_line logical. Whether or not to the zero #' medians.
#' @param horizontal logical. Whether or not you would like the lines to be
#' @param xlims custom xlims
#' @importFrom ggplot2 ggplot
#' @examples
#' #plot covariate catepillar plot
#' library(ggplot2)
#' ggplot(test.bbgdm)
#'
#' @export

ggplot.bbgdm <- function(object, pars=NULL, pars_labels = NULL,
                         horizontal = TRUE,
                         zero_line=TRUE,
                         xlims = NULL,...){
  # Extract all simulations
  sims <- as.data.frame(object$all.coefs.se)

  # Extract only desired parameters
  par_names <- names(sims)
  if(is.null(pars)){
    sims_subset <- sims
    cat('plotting all spline bases \n')
  } else {
    sims_subset <- sims[, par_names %in% pars]
  }
  if (ncol(sims_subset) == 0) {
    stop("No parameters selected. \n", call. = FALSE)
  }
  # Gather for plotting
  gathered <- tidyr::gather(sims_subset, variable, value)

  # Add labels
  if (!is.null(pars_labels)) {
    message("\nEnsure that your parameter labels are in the same order as the parameters.\n")
    if (length(pars_labels) !=
        length(unique(gathered$variable))) {
      stop("pars_labels must equal the number of plotted parameters.",
           call. = FALSE)
    }
    gathered$variable <- factor(gathered$variable,
                                labels = pars_labels)
  }

  # Find central interval
  gathered <- dplyr::group_by(gathered, variable)
  lower95 <- dplyr::summarise(gathered, quantile(value, 0.025))
  lower90 <- dplyr::summarise(gathered, quantile(value, 0.05))
  upper90 <- dplyr::summarise(gathered, quantile(value, 0.95))
  upper95 <- dplyr::summarise(gathered, quantile(value, 0.975))

  # Find medians
  medians <- dplyr::summarise(gathered, median(value))

  # Merge
  comb <- suppressMessages(dplyr::inner_join(lower95, lower90))
  comb <- suppressMessages(dplyr::inner_join(comb, medians))
  comb <- suppressMessages(dplyr::inner_join(comb, upper90))
  comb <- suppressMessages(dplyr::inner_join(comb, upper95))
  names(comb) <- c('pars', 'lower95', 'lower90', 'medians', 'upper90',
                   'upper95')

  # Plot
  pp <- ggplot2::ggplot(comb, ggplot2::aes(x = medians, y = pars,
                                           xmin = lower95,
                                           xmax = upper95),...) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_segment(ggplot2::aes(x = lower95, xend = upper95, yend = pars),
                          size = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = lower90, xend = upper90,
                                       yend = pars), size = 1.5) +
    ggplot2::xlab('parameter quantiles from bayesian bootstrap samples')+
    ggplot2::ylab('') +
    if(!is.null(xlims))ggplot2::xlim(xlims)+
    ggthemes::theme_few()+
    ggplot2::theme(plot.margin=unit(c(5,5,5,5),"mm"))

  if (isTRUE(zero_line)) pp <- pp + ggplot2::geom_vline(xintercept = 0)
  if (!isTRUE(horizontal)) pp <- pp + ggplot2::coord_flip()
  return(pp)
}

utils::globalVariables(base::c("variable", "value"))

