#include "Rcpp.h"
using namespace Rcpp;

//' diff_table cpp function
//' @param env_dat matrix of sites x environmental covariates.
//' @return diff_table difference between covariates between site-pairs
//' @export
//[[Rcpp::export]]

NumericMatrix diff_table_cpp(NumericMatrix env_dat){

      double diff_var;
      int nr = env_dat.nrow();
      //int nc = env_dat.ncol();
      int nr_dt = (nr*(nr-1))/2;
      int nc_dt = env_dat.ncol();
      //Rcpp::Rcout <<  nr << " " << nr_dt << " " << nc_dt << std::endl;
      NumericMatrix diff_table(nr_dt,nc_dt);
      int pair = 0;

        for(int i_site = 0; i_site<(nr-1); i_site++)
          {
            for(int j_site = i_site+1; j_site<nr; j_site++)
              {
                for(int var = 0; var<nc_dt; var++)
                  {
                    diff_var = std::abs(env_dat(i_site,var))-std::abs(env_dat(j_site,var));
                    diff_table(pair,var) = std::abs(diff_var);
                   }
                    pair++;
               }
                  //printf("%d ", pair);
          }
return diff_table;
}
