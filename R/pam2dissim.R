#' Calculate Dissimilarities from Presence Absence Matrix
#' 
#' Creates a vector of dissimilarities from presence absence matrix.
#' @param x Presence Absence Matrix.
#' @param dissim If "bray_curtis" calculates Bray-Curtis dissimilarity, if "number_shared " returns number of shared species at site_ij which can be used in a binomial model for GDM.
#' @return dism a vector of dissimilarities or count of shared species calcualted for the upper-triangle of a dissimilarity matrix.
#' @export

pam2dissim <- function(x,dissim="bray_curtis"){
  if(dissim=="bray_curtis"){
#     if(!Abund){
    cat("Calculating Bray-Curtis Dissimilarity.\n")
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
# #     } else {
# #       x <- as.matrix(x)
#       result <- matrix(nrow = nrow(x), ncol = nrow(x))
#       rownames(result) <- rownames(x)
#       colnames(result) <- rownames(x)
#       for (i in 1:nrow(x)) {
#         for (j in i:nrow(x)) {
#           A <- sum(pmin(x[i, ], x[j, ]))
#           B <- sum(x[i, ]) - sum(pmin(x[i, ], x[j, ]))
#           C <- sum(x[j, ]) - sum(pmin(x[i, ], x[j, ]))
#           result[i, j] <- min(B, C)/(A + min(B, C))
#           result[j, i] <- (B + C)/(2 * A + B + C)
#         }
#       }
      dism <- dissimilarity.data
  } 
  if(dissim=="simpsons"){
    cat("Calculating Simpsons Dissimilarity.\n")
    x <- as.matrix(x)
    shared <- x %*% t(x)
    not.shared <- abs(sweep(shared, 2, diag(shared))) ## end for i_site
    min.not.shared <- pmin(not.shared, t(not.shared))
    beta.sim <- min.not.shared/(min.not.shared + shared)
    dism <- beta.sim
   }
  if(dissim=="number_shared") {
    cat("Calculating Number of Shared Species for Binomial Model.\n")
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
  if(dissim=="number_shared") {
    dism <- list(no_share=no_share,sum_share=sum_share)#,max_sp=max_sp)
    return(dism)
  }else{
    return(dism)
  }
}