#' @title response objects
#' @rdname response
#' @name as.response
#' @param object As derived from bbgdm function
#' @return values for plotting i-spline responses from bbgdm
#' @export
#' @examples
#' x <-matrix(rbinom(1:100,1,.6),10,10)# presence absence matrix
#' y <- simulate_covariates(x,2)
#' form <- ~ 1 + covar_1 + covar_2
#' test.bbgdm <- bbgdm(form, sp.dat=x, env.dat=y,family="binomial",
#'                    geo=FALSE,dism_metric="number_non_shared", nboot=10)
#' responses <- as.response(test.bbgdm)


as.response <- function(object, ...){
  Xold <- data.matrix(object$dissim_dat)
  splineLength <- sapply(object$dissim_dat_params, `[[`, "dim")[2,]
  betas <- object$median.coefs.se[2:length(object$median.coefs.se)]
  betas.quantiles <- object$quantiles.coefs.se[,2:length(object$starting_gdm$coef),drop=FALSE]
  if(object$geo){ k <- ncol(object$env.dat)-1
                  coords <- as.matrix(object$env.dat[,1:2])
                  lc_data<-object$lc_data
                  minr<-object$minr
                  maxr<-object$maxr
                  geos <- calc_geo_dist(coords,geo.type=object$geo.type)
                  nr_df<-((nrow(object$sp.dat)^2)-nrow(object$sp.dat))/2
                  nc_dt<-ncol(env.dat)
                  ne<-ncol(env.dat)
                  diff_table <- diff_table_cpp(as.matrix(env.dat))
                  diff_table <- cbind(geos[,5],diff_table)
                  colnames(diff_table) <-c('geo',colnames(env.dat))
  } else {
    k <- ncol(object$env.dat)
    env.dat <- object$env.dat
    nr_df<-((nrow(object$sp.dat)^2)-nrow(object$sp.dat))/2
    nc_dt<-ncol(env.dat)
    ne<-ncol(env.dat)
    diff_table <- diff_table_cpp(as.matrix(env.dat))
    colnames(diff_table) <-c(colnames(env.dat))
    }
      grid <- matrix(rep(seq(0, 1, length.out = 100), k), ncol = k)
      grid <- t(t(grid) * as.vector(diff(apply(diff_table, 2, range))) + apply(diff_table, 2, min))
      grid <- data.frame(grid)
      min_env <- apply(env.dat,2,function(x)min(abs(x)))
      max_env <- apply(env.dat,2,function(x)max(abs(x)))
      if(object$geo){
         min_env <- c(min(grid[,1]),min_env)
         max_env <- c(max(grid[,1]),max_env)
      }
      grid_real <- grid
      for(i in 1:ncol(grid)) grid_real[,i] <-  scales::rescale(grid[,i],to=c(min_env[i],max_env[i]))
      X <- mapply(spline_trans_for_pred, grid, attrib = object$dissim_dat_params,
                  SIMPLIFY = FALSE)
      l <- sapply(object$dissim_dat_params, `[[`, "dim")[2, ]
      at <- unlist(mapply(rep, 1:k, l))
      bspl <- lapply(1:k, function(i) betas[at == i])
      bspl.05 <- lapply(1:k, function(i) betas.quantiles[1,at == i])
      bspl.95 <- lapply(1:k, function(i) betas.quantiles[2,at == i])
      structure(list(bspl=bspl,bspl.05=bspl.05,bspl.95=bspl.95,X=X,
                     grid_real=grid_real,diff_table=diff_table),class='response')
      }

#' @rdname response
#' @name plot.response
#' @export
#' @examples
#' #plot responses
#' par(mfrow=c(1,2))
#' plot(responses)

plot.response <- function(object,...){
        Splinesum <- Spline.05 <- Spline.95 <- NULL
        Splinessum <- mapply(`%*%`, object$X, object$bspl)
        Splines.05 <- mapply(`%*%`, object$X, object$bspl.05)
        Splines.95 <- mapply(`%*%`, object$X, object$bspl.95)
         for (i in 1:ncol(Splinessum)) {
              plot(object$grid_real[,i], Splinessum[,i],type='l', ylab = paste0("f(",colnames(object$diff_table)[i],")"),
                   xlab = colnames(object$diff_table)[i],ylim = range(c(Splinessum,Splines.05,Splines.95)),...)
              polygon(c(object$grid_real[, i],rev(object$grid_real[,i])),c(Splines.05[,i],rev(Splines.95[,i])),col="grey70",border=NA)
              lines(object$grid_real[,i], Splinessum[,i], col = "black",type = "l", lwd=2)
              # lines(object$grid_real[, i], Splines.05[,i], col = "black",type = "l", lty = 3)
              # lines(object$grid_real[, i], Splines.95[,i], col = "black",type = "l", lty = 3)
              mtext(paste0("(",letters[i],")"),adj = 0)
         }
}
