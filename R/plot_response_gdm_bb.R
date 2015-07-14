#' Function for plotting bb.GDM response curves.
#'
#' Plots a Generalised dissimilarity model with bayesian bootstrap.
#' @param gdm.bb_object As derived from gdm.bb function
#' @param plot.layout par settings for plot. default is c(2,2) (four plots).
#' @return a plots from gdm.bb
#' @export
#' @examples
#' x<-matrix(rbinom(1:100,1,.6),10,10)# presence absence matrix
#' #' env.dat <- simulate_covariates(x,2)
#' form <- y ~ covar_1 + covar_2
#' test.gdm.bb <- gdm.bb(form,family="binomial", dism_metric="number_shared", nboot=100)
#' plotResponse(test.gdm.bb,plot.error=T)


plotResponse <- function(object, plotdim = c(2, 2),plot.error=F,monotonic=F,...){
  nplot <- prod(plotdim)
  Xold <- data.matrix(object$dissim_dat)
  splineLength <- sapply(object$dissim_dat_params, `[[`, "dim")[2,]
  if(object$family=='beta_dist'){
    betas <- object$median.coefs.se[2:(length(coef(object$starting_gdm))-1)]
    betas.quantiles <- object$quantiles.coefs.se[,2:(length(coef(object$starting_gdm))-1)]
    } else {
  betas <- object$mean.coefs.se[2:length(object$starting_gdm$coef)]
  betas.quantiles <- object$quantiles.coefs.se[,2:length(object$starting_gdm$coef)]
  }
  if(monotonic){
    betas<- ifelse(betas<0,0,betas)
    betas.quantiles <- ifelse(betas.quantiles<0,0,betas.quantiles)+1e-06
  }
  if(object$geo){ k <- ncol(object$env.dat)-1
                  coords <- as.matrix(object$env.dat[,1:2])
                  lc_data<-object$lc_data
                  minr<-object$minr
                  maxr<-object$maxr
                  geos <- calc_geo_dist(coords,geo.type=object$geo.type,lc_data=lc_data,minr=minr,maxr=maxr) 
                  env.dat <- object$env.dat[,-c(1:2)]
                  nr_df<-((nrow(object$sp.dat)^2)-nrow(object$sp.dat))/2
                  nc_dt<-ncol(env.dat)
                  ne<-ncol(env.dat)
                  diff_table<-matrix(NA,nr_df,nc_dt)  
                  pair<-1
                  for(i_site in 1:(nrow(env.dat)-1))
                  {
                    for(j_site in (i_site+1):nrow(env.dat))
                    {
                      for(var in 1:ne)
                      {
                        diff_table[pair,var]<-abs(env.dat[i_site,var]-env.dat[j_site,var])
                      } 
                      pair<-pair+1 
                    } 
                  }
                  diff_table <- cbind(geos[,5],diff_table)
                  colnames(diff_table) <-c('geo',colnames(env.dat))
  } else { 
    k <- ncol(object$env.dat)
    env.dat <- object$env.dat
    nr_df<-((nrow(object$sp.dat)^2)-nrow(object$sp.dat))/2
    nc_dt<-ncol(env.dat)
    ne<-ncol(env.dat)
    diff_table<-matrix(NA,nr_df,nc_dt)  
    pair<-1
    for(i_site in 1:(nrow(env.dat)-1))
    {
      for(j_site in (i_site+1):nrow(env.dat))
      {
        for(var in 1:ne)
          {
            diff_table[pair,var]<-abs(env.dat[i_site,var]-env.dat[j_site,var])
          } 
          pair<-pair+1 
        } 
      }
    colnames(diff_table) <-c(colnames(env.dat))
    }
      grid <- matrix(rep(seq(0, 1, length.out = 100), k), ncol = k)
      grid <- t(t(grid) * as.vector(diff(apply(diff_table, 2, range))) + apply(diff_table, 2, min))
      grid <- data.frame(grid)
      X <- mapply(spline_trans_for_pred, grid, attrib = object$dissim_dat_params, 
                  SIMPLIFY = FALSE)
      l <- sapply(object$dissim_dat_params, `[[`, "dim")[2, ]
      at <- unlist(mapply(rep, 1:k, l))
      bspl <- lapply(1:k, function(i) betas[at == i])
      bspl.05 <- lapply(1:k, function(i) betas.quantiles[1,at == i])
      bspl.95 <- lapply(1:k, function(i) betas.quantiles[2,at == i])
      Splinesum <- Spline.05 <- Spline.95 <- NULL
      for (j in 1:ceiling(k/nplot)) {
        mini <- nplot * (j - 1) + 1
        maxi <- min(k, j * nplot)
        noi <- maxi - mini + 1
        plotlayout <- rep(0, nplot)
        plotlayout[1:noi] <- 1:noi
        plotlayout <- rbind(matrix(plotlayout, byrow = TRUE, 
                                   nrow = plotdim[1], ncol = plotdim[2]), noi + 1)
        layout(plotlayout, heights = c(rep(1, plotdim[1]), 0.2))
        par(ask = j > 1, mar = c(2.1, 4.1, 3.1, 2.1), mgp = c(2,1, 0))
        Splinessum <- mapply(`%*%`, X, bspl)
        Splines.05 <- mapply(`%*%`, X, bspl.05)
        Splines.95 <- mapply(`%*%`, X, bspl.95)
         if (is.matrix(Splinessum)) {
           Splinessum <- data.frame(Splinessum)
           Splines.05 <- data.frame(Splines.05)
           Splines.95 <- data.frame(Splines.95)
         }
          if(plot.error){
            for (i in mini:maxi) {
              plot(grid[, i], Splinessum[,i], col = "black",type = "l", ylab = paste0("f(",colnames(diff_table)[i],")"),
                   ylim = range(c(Splinessum[,i],Splines.05[,i],Splines.95[,i])), main = colnames(diff_table)[i],xlab = "",...)#,Splines.pred.05[,i],Splines.pred.95[,i]))
              lines(grid[, i], Splines.05[,i], col = "black",type = "l", lty = 3,...)
              lines(grid[, i], Splines.95[,i], col = "black",type = "l", lty = 3,...)
                } 
              par(mar = c(0, 0, 0, 0),...)
              plot(0, 0, type = "n", axes = FALSE, bty = "n", pty = "m",xlab = "", ylab = "",...)
              legend("center", c("BB quantiles estimates","BB mean estimates"), horiz = TRUE, lty = c(3,1), col = c("black","black"), bg = "white", bty = "o", xpd = T,...)
            } else {
              for (i in mini:maxi) {
                  plot(grid[, i], Splinessum[,i], col = "black",type = "l", ylab = paste0("f(",colnames(diff_table)[i],")"),
                       ylim = range(c(Splinessum[,i])), main = colnames(diff_table)[i],xlab = "",...)#,Splines.pred.05[,i],Splines.pred.95[,i]))
              }
              par(mar = c(0, 0, 0, 0),...)
              plot(0, 0, type = "n", axes = FALSE, bty = "n", pty = "m",xlab = "", ylab = "",...)
              legend("center", c("BB mean spline"), horiz = TRUE, lty = c(1), col = c("black"), bg = "white", bty = "o", xpd = T,...)
            }
          }
    par(mfrow = c(1, 1), mgp = c(3, 1, 0), mar = c(4, 5, 5,2) + 0.1,...)
}
