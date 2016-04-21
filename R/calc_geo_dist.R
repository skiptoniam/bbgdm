#' Function that calculates geographic distances.
#' 
#' Functions for estimating geographic distances.
#' @param coords A set of lats and lons for sites
#' @param geo.type geo.type of geographic distance to estimate, 'euclidean', 'greater circle' and 'least cost pathway'. 
#' @param lc_data is a matrix of lons as rows and lats as cols, with a values of gradient as z values (eg,depth,elevation) or a raster of depth, elevation ect.
#' @param minr The minimum range to set least cost pathway, defaul is 0 to max of range.
#' @param maxr The maximum range to set least cost pathway, defaul is NULL to max of range.
#' @return distances Five column matrix of sites pairs (cols 1:4) and distance as col 5.
#' @export


calc_geo_dist <- function(coords,geo.type=c('euclidean','greater_circle','least_cost'),
                          lc_data=NULL,minr=0,maxr=NULL){
        
  geo.type <- match.arg(geo.type,c('euclidean','greater_circle','least_cost'))
      if(geo.type=="euclidean"){
          cat('calculating euclidean distance\n')
          nr_df<-((nrow(coords)^2)-nrow(coords))/2
          nc_dt<-2+(ncol(coords))
          diff_table<-matrix(NA,nr_df,nc_dt+1)
          pair<-1
            for(i_site in 1:(nrow(coords)-1))
            {
              for(j_site in (i_site+1):nrow(coords))
              {
                diff_table[pair,c(1,2)]<-as.numeric(coords[i_site,])
                diff_table[pair,3:4]<-as.numeric(coords[j_site,])
                diff_table[pair,5]<-sqrt(sum((coords[i_site,]-coords[j_site,])^2))
                pair<-pair+1 
              } 
            }
          }
          if(geo.type=="greater_circle"){
            cat('calculating greater circle distance\n')
            nr_df<-((nrow(coords)^2)-nrow(coords))/2
            nc_dt<-2+(ncol(coords))
            diff_table<-matrix(NA,nr_df,nc_dt+1)
            pair<-1
            for(i_site in 1:(nrow(coords)-1))
            {
              for(j_site in (i_site+1):nrow(coords))
              {
                diff_table[pair,1:2]<-as.numeric(coords[i_site,])
                diff_table[pair,3:4]<-as.numeric(coords[j_site,])
                diff_table[pair,5]<-gcdist(coords[i_site,1],coords[i_site,2],coords[j_site,1],coords[i_site,2])
                pair<-pair+1 
              } 
            }
          }
          if(geo.type=="least_cost"){
            cat('calculating least cost distance\n')
          if(is.null(lc_data)) stop("Please include matrix or raster to build least cost pathway, see help")
          cat('This might take a while, grab a cuppa. \n')
          if(is.null(maxr)) cat('least cost estimated between', minr,' and max range of raster values\n')
          else cat('least cost estimated between',minr,' and ',maxr,'range of matrix/raster values\n')
          tr <- trans_rast_fun(lc_data, min.range = minr, max.range = maxr)   
          lc <- least.cost.dist(tr,coords)
          lcm <- as.matrix(lc)
          nr_df<-((nrow(coords)^2)-nrow(coords))/2
          nc_dt<-2+(ncol(coords))
          diff_table<-matrix(NA,nr_df,nc_dt+1)
          pair<-1
            for(i_site in 1:(nrow(coords)-1))
            {
              for(j_site in (i_site+1):nrow(coords))
              {
                diff_table[pair,1:2]<-as.numeric(coords[i_site,])
                diff_table[pair,3:4]<-as.numeric(coords[j_site,])
                diff_table[pair,5]<-lcm[i_site,j_site]
                pair<-pair+1 
              } 
            }
          }
          colnames(diff_table) <- c("X0","Y0","X1","Y1",geo.type)
          diff_table[diff_table[,5]==0,5]<-1e-06
          return(diff_table)
}
