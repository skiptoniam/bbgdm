#' @title data.prep objects
#' @rdname data.prep
#' @name table2pam
#' @param x data to be imported into R. A table of species occurrences and site records.
#' @param site.id Name of column that contains site names /IDs.
#' @param sp.id Name of column that contains species information/ID.
#' @param abund LOGICAL if TRUE calculate abundances.
#' @param abund.col Name of column that contains abundance/counts data.
#' @param siteXsp LOGICAL if TRUE Returns sites as rows, sites as columns, if FALSE returns species as rows and sites as cols.
#' @return A species by sites matrix for use in dissimilarity calculation
#' @export


table2pam <- function (x, site.id = "site.id", sp.id = "sp.id", abund = FALSE, abund.col = " No.of.specimens",siteXsp=TRUE)
{
  a <- site.id
  nr <- length(levels(as.factor(x[, a])))
  rn <- levels(as.factor(x[, a]))
  z <- sp.id
  cn <- levels(as.factor(x[, z]))
  nc <- length(cn)
  nm <- matrix(0, nr, nc, dimnames = list(rn, cn))
  for (i in 1:length(x[, 1])) {
    m <- as.character(x[i, a])
    n <- as.character(x[i, z])
    if (is.na(m) == TRUE | is.null(m) == TRUE | is.na(n) ==
        TRUE | is.null(n) == TRUE)
      (next)(i)
    if (m == "" | m == " " | n == "" | n == " ")
      (next)(i)
    if (abund == TRUE)
      nm[m, n] <- nm[m, n] + x[i, abund.col]
    else nm[m, n] <- 1
  }
  fm <- nm[rowSums(nm) > 0, ]
  if(siteXsp){ return(as.matrix(fm))
  } else {
    return(as.matrix(t(fm)))
  }
}

#' @title data.prep objects
#' @rdname data.prep
#' @name pam2dissim
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

#' @title data.prep objects
#' @rdname data.prep
#' @name simulate_covariates
#' @export
#' @examples
#' x <- matrix(rbinom(16,1,.6),4,4) #toy presence absence matrix
#' covr <- simulate_covariates(x,6)

simulate_covariates <- function(x,N){
  Nsite <- nrow(x)
  Ncov <- N
  beta <- matrix(rnorm(Nsite*Ncov, 0, 2),nrow=Nsite, ncol=Ncov)
  colnames(beta) <- paste0("covar_",1:Ncov)
  return(as.data.frame(beta))
}
