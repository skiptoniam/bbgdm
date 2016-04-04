#' Calculate Dissimilarities from Presence Absence Matrix
#' 
#' Creates a vector of dissimilarities from presence absence matrix.
#' @param x Presence Absence Matrix.
#' @param dism_metric If "bray_curtis" calculates Bray-Curtis dissimilarity, if "number_non_shared " returns number of shared species at site_ij which can be used in a binomial model for GDM.
#' @return dism a vector of dissimilarities or count of shared species calcualted for the upper-triangle of a dissimilarity matrix.
#' @export

pam2dissim <- function(x,dism_metric="number_non_shared"){
    cat(paste0("Calculating ",dism_metric," dissimilarity\n"))  
    if(dism_metric=="bray_curtis"){
#     if(!Abund){
    x <- as.matrix(x)
    dissimilarity.data<-matrix(NA,nrow(x),nrow(x))
    for(i_site in 1:nrow(x))
    {
      richness_i_site<-0
      for(i_spp in 1:ncol(x))
      {
        richness_i_site<-richness_i_site + x[i_site,i_spp] 
      } ## end for i_spp
      for(j_site in 1:nrow(x))
      {
        richness_j_site<-0    
        sharedspecies<-0
        for(i_spp in 1:ncol(x))
        {
          richness_j_site<-richness_j_site + x[j_site,i_spp]
          if(x[i_site,i_spp]*x[j_site,i_spp] == 1)
          {      
            sharedspecies<-sharedspecies + 1
          } ## end if  
        } ## end for i_spp    
        dissimilarity.data[i_site,j_site]<-(1 - ((2 * sharedspecies)/(richness_i_site + richness_j_site)))
      } ## end for j_site
    } ## end for i_site
      dism <- round(dissimilarity.data*100)
  } 
#   if(dissim=="simpsons"){
#     cat("Calculating Simpsons Dissimilarity.\n")
#     x <- as.matrix(x)
#     shared <- x %*% t(x)
#     not.shared <- abs(sweep(shared, 2, diag(shared))) ## end for i_site
#     min.not.shared <- pmin(not.shared, t(not.shared))
#     beta.sim <- min.not.shared/(min.not.shared + shared)
#     dism <- beta.sim
#    }
  if(dism_metric=="number_non_shared") {
    # cat("Calculating Number of Shared Species for Binomial Model.\n")
    x <- as.matrix(x)
    no_share <- matrix(NA,nrow(x),nrow(x))
    sum_share <- matrix(NA,nrow(x),nrow(x))
#     max_sp <- matrix(NA,nrow(x),nrow(x))
    for(i_site in 1:(nrow(x)-1)){
      for(j_site in c(i_site+1):c(nrow(x))){
#         no_share[i_site,j_site] <- sum(apply(x[c(i_site,j_site),],2,sum)==2)
        no_share[i_site,j_site] <- sum(apply(x[c(i_site,j_site),],2,sum)==1)#orig
        sum_share[i_site,j_site] <- sum(apply(x[c(i_site,j_site),],2,sum)!=0)
#         max_sp[i_site,j_site] <- max(apply(x[c(i_site,j_site),],1,sum))
      }
    }
  }
  if(dism_metric=="number_non_shared") {
    dism <- list(no_share=no_share,sum_share=sum_share)#,max_sp=max_sp)
    return(dism)
  }else{
    return(dism)
  }
}