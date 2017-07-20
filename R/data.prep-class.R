#' @title data.prep objects
#' @rdname data.prep
#' @name table2pam
#' @description prepares table of occurrence records or precence-absence matrix for use in \code{bbgdm}.
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
#' @param y presence-absence matrix of site x sp.
#' @param dism_metric If "bray_curtis" calculates Bray-Curtis dissimilarity, if "number_non_shared " returns number of shared species at site_ij which can be used in a binomial model for GDM.
#' @return dism a vector of dissimilarities or count of shared species calcualted for the upper-triangle of a dissimilarity matrix.
#' @export

pam2dissim <- function(y,dism_metric="number_non_shared"){

  if(!length(which(y>1))==0)stop('bbgdm only works with presence-absence - transform your abundances to binary data\n')
  cat(paste0("Calculating ",dism_metric," dissimilarity\n"))
  if(dism_metric=="bray_curtis"){
    y <- as.matrix(y)
    dissimilarity.data<-matrix(NA,nrow(y),nrow(y))
    for(i_site in 1:nrow(y))
    {
      richness_i_site<-0
      for(i_spp in 1:ncol(y))
      {
        richness_i_site<-richness_i_site + y[i_site,i_spp]
      } ## end for i_spp
      for(j_site in 1:nrow(y))
      {
        richness_j_site<-0
        sharedspecies<-0
        for(i_spp in 1:ncol(y))
        {
          richness_j_site<-richness_j_site + y[j_site,i_spp]
          if(y[i_site,i_spp]*y[j_site,i_spp] == 1)
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
    y <- as.matrix(y)
    no_share <- matrix(NA,nrow(y),nrow(y))
    sum_share <- matrix(NA,nrow(y),nrow(y))
     for(i_site in 1:(nrow(y)-1)){
      for(j_site in c(i_site+1):c(nrow(y))){
         no_share[i_site,j_site] <- sum(apply(y[c(i_site,j_site),],2,sum)==1)
        sum_share[i_site,j_site] <- sum(apply(y[c(i_site,j_site),],2,sum)!=0)
       }
    }
  }
  if(dism_metric=="number_non_shared") {
    dism <- list(no_share=no_share,sum_share=sum_share)
    return(dism)
  }else{
    return(dism)
  }
}

#' @title data.prep objects
#' @rdname data.prep
#' @name simulate_covariates
#' @param n number of covariates to simulate
#' @export
#' @examples
#' x <- matrix(rbinom(16,1,.6),4,4) #toy presence absence matrix
#' covr <- simulate_covariates(x,6)

simulate_covariates <- function(y,n){
  nsite <- nrow(y)
  ncov <- n
  beta <- matrix(rnorm(nsite*ncov, 0, 2),nrow=nsite, ncol=ncov)
  colnames(beta) <- paste0("covar_",1:ncov)
  return(as.data.frame(beta))
}
