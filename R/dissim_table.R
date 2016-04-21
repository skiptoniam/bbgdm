#' Create dissimilarity table for bbGDM.
#'
#' Calculates dissimilarity table for bbGMD from presences absence matrix (pam).
#' @param sp.dat Presence absence matrix.
#' @param env.dat Matrix of covariates for bbGDM; if using own data as covariates make sure to set sim_covar=FALSE,
#' otherwise simulated covariates will be caluclated instead.
#' @param dism_metric Dissimilarity to caluclate; If "bray_curtis" calculates Bray-Curtis dissimilarity,
#' if "number_non_shared" returns number of shared species at site_ij which can be used in a binomial model for GDM.
#' @param spline_type If "bspline" calculates bs spline from spline package. If "ispline" calculates ispline.
#' @param spline_df degrees of freedom; one can specify df rather than knots. Default = 3.
#' @param spline_knots The internal breakpoints that define the spline. Default = NULL for bs spline and 1 for ispline.
#' @param geo logical If true geographic distance is calculated if
#' @param geo.type type of geographic distance to estimate, can call 'euclidean','greater_circle' or 'least_cost'. If least_cost is called extra parameters are required (lc_data, minr and maxr).
#' @param coord.names c("X","Y") character.vector names of coordinates.
#' @param lc_data NULL lc_cost data layer, in the form of a raster.
#' @param minr NULL range of values for marine data within the scope of the lc_cost raster. eg. min depth.
#' @param maxr NULL range of values for marine data within the scope of the lc_cost raster. eg. max depth.

#' @return diff_table dissimilarity table; creates a table of dissimilarities, and difference of covariates between each sites_ij.
#'  If select dissim="number_non_shared", will return "nonsharedspp_ij","sumspp_ij",
#'  which can be used as response variables in a binomial bbGDM.
#' @export
#' @examples
#' x <- matrix(rbinom(100,1,.5),10,10) #toy presence absence matrix
#' y <- simulate_covariates(x,4)
#' diff_table <- dissim_table(x,y,dism_metric="number_non_shared")

dissim_table <- function(sp.dat,env.dat,dism_metric="number_non_shared",spline_type="bspline",spline_df=1,
                         spline_knots=2, coord.names=c("X","Y"),
                         geo=FALSE,geo.type="euclidean",lc_data=NULL,minr=NULL,maxr=NULL){
  if(!is.matrix(env.dat)) env.dat <- as.matrix(env.dat)
  if(geo){
    coords <- as.matrix(env.dat[,which(colnames(env.dat)%in%coord.names)])
    geos <- calc_geo_dist(coords,geo.type=geo.type,lc_data=lc_data,minr=minr,maxr=maxr)
    env.dat <- env.dat[,-which(colnames(env.dat)%in%coord.names),drop=FALSE]
  }
    if(dism_metric=="bray_curtis"){
    xdism <- pam2dissim(sp.dat,dism_metric)
    nr_df<-((nrow(sp.dat)^2)-nrow(sp.dat))/2
    nc_dt<-2+(ncol(env.dat))
    ne<-ncol(env.dat)
    diff_table<-matrix(NA,nr_df,nc_dt)
    colnames(diff_table)<-c("ID","dissimilarity",colnames(env.dat))
    pair<-1
    for(i_site in 1:(ncol(xdism)-1))
    {
      for(j_site in (i_site+1):nrow(xdism))
      {
        diff_table[pair,1]<-pair
        diff_table[pair,2]<-xdism[i_site,j_site]
        for(var in 1:ne)
        {
          diff_table[pair,(2+var)]<-abs(env.dat[i_site,var]-env.dat[j_site,var])
        }
        pair<-pair+1
      }
    }
  }
  if(dism_metric=="number_non_shared"){
    xdism <- pam2dissim(sp.dat,dism_metric)
    nr_df<-((nrow(sp.dat)^2)-nrow(sp.dat))/2
    nc_dt<-3+(ncol(env.dat))
    ne<-ncol(env.dat)
    diff_table<-matrix(NA,nr_df,nc_dt)
    maxsppsite <- apply(sp.dat,1,sum)
    colnames(diff_table)<-c("ID","nonsharedspp_ij","sumspp_ij",colnames(env.dat))
    pair<-1
    for(i_site in 1:(nrow(xdism$no_share)-1))
    {
      for(j_site in (i_site+1):nrow(xdism$no_share))
      {
        diff_table[pair,1]<-pair
        diff_table[pair,2]<-xdism$no_share[i_site,j_site]
        diff_table[pair,3]<-xdism$sum_share[i_site,j_site]
#         diff_table[pair,3]<-max(maxsppsite[i_site],maxsppsite[j_site])
        for(var in 1:ne)
        {
          diff_table[pair,(3+var)]<-abs(env.dat[i_site,var]-env.dat[j_site,var])
        }
        pair<-pair+1
      }
    }
  }
  cat("Transforming covariates to ",spline_type," with ",spline_df,"degrees of freedom.\n")
  if(dism_metric=='number_non_shared'){
  if(geo){
    diff_table <- cbind(diff_table[,2:3],geos=geos[,5],diff_table[,4:ncol(diff_table),drop=FALSE])
    tmp <- spline.trans(x=diff_table[,3:ncol(diff_table)],spline_type=spline_type,spline_df=spline_df,spline_knots=spline_knots)
    diff_table_final <- tmp$spline
    diff_table_final <- cbind(diff_table[,1:2],diff_table_final)
    diff_table_params <- tmp$spline.attr
    } else {
    tmp <- spline.trans(x=diff_table[,4:ncol(diff_table)],spline_type=spline_type,spline_df=spline_df,spline_knots=spline_knots)
    diff_table_final <- tmp$spline
    diff_table_final <- cbind(diff_table[,2:3],diff_table_final)
    diff_table_params <- tmp$spline.attr
    }
  } else {
    if(geo){
      diff_table <- cbind(dissimilarity=diff_table[,2],geos=geos[,5],diff_table[,3:ncol(diff_table)])
      tmp<-spline.trans(x=diff_table[,2:ncol(diff_table)],spline_type=spline_type,spline_df=spline_df,spline_knots=spline_knots)
      diff_table_final <- tmp$spline
      diff_table_final <- cbind(dissimilarity=diff_table[,1],diff_table_final)
      diff_table_params <- tmp$spline.attr
      } else {
      tmp <- spline.trans(x=diff_table[,3:ncol(diff_table)],spline_type=spline_type,spline_df=spline_df,spline_knots=spline_knots)
      diff_table_final <- tmp$spline
      diff_table_final <- cbind(dissimilarity=diff_table[,2],diff_table_final)
      diff_table_params <- tmp$spline.attr
      }
  }
 return(list(diff_table=diff_table_final,diff_table_params=diff_table_params))
}

