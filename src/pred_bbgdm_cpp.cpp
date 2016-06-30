#include "RcppArmadillo.h"
using namespace Rcpp;

//' predict bbgdm cpp function
//' @param raster_stack raster stack of i-spline transformed environmental layers.
//' @param output_raster binary raster one for cells to predict to and NA for areas not to predict too.
//' @param neighbourhood numeric matrix which creates the neighbourhood average.
//' @param bbgdm_coef I-spline coefficents from bbgdm model fit.
//' @export
// [[Rcpp::export]]

NumericMatrix pred_bbgdm_cpp(NumericVector raster_stack, NumericMatrix output_raster, NumericMatrix neighbourhood, NumericVector bbgdm_coef){

    NumericVector vecArray(raster_stack);
    IntegerVector arrayDims = vecArray.attr("dim");
    arma::cube raster_layers(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
    arma::mat raster = as<arma::mat>(output_raster);
    int ncol = raster.n_cols;
    int nrow = raster.n_rows;
    arma::mat beta_raster(nrow,ncol);
    beta_raster.fill(NA_REAL);
    arma::mat rel = as<arma::mat>(neighbourhood);
    int n_rel = rel.n_rows;
    arma::vec coef = as<arma::vec>(bbgdm_coef);
    int ncoef = coef.size();
    arma::vec prob(ncoef);
    std::vector<double> site_one_trans(arrayDims[2]), site_two_trans(arrayDims[2]);//, site_one_two_trans(arrayDims[2]);
    double dissimilarity_one_two, mean_weight_one_two, d_dissimilarity_sum, d_weighted_pair_sum;

    for(int i_row = 0; i_row<nrow; i_row++)
    {
      for(int i_col = 0; i_col<ncol; i_col++)
      {
        if(arma::is_finite(raster(i_row,i_col)))
        {
          // Now loop through each pair of cells in the neighbourhood around the focal cell
          d_dissimilarity_sum = 0;
          d_weighted_pair_sum = 0;
          for(int i_site_one = 0; i_site_one<(n_rel-2); i_site_one++)
          {
            int i_row_site_one = i_row + rel(i_site_one,0);
            int i_col_site_one = i_col + rel(i_site_one,1);
            if(i_row_site_one > -1 && i_row_site_one < nrow && i_col_site_one > -1 && i_col_site_one < ncol)
            {
              if(arma::is_finite(raster(i_row_site_one,i_col_site_one)))
              {
                for(int i_site_two = i_site_one; i_site_two<(n_rel-1); i_site_two++)
                {
                  int i_row_site_two = i_row + rel(i_site_two,0);
                  int i_col_site_two = i_col + rel(i_site_two,1);
                  //Rcpp::Rcout <<  i_site_one << " " << i_site_two << std::endl;
                  //Rcpp::Rcout <<  i_row_site_one << " " << i_col_site_one << " " << i_row_site_two << " " << i_col_site_two << std::endl;
                  if(i_row_site_two > -1 && i_row_site_two < nrow && i_col_site_two > -1 && i_col_site_two < ncol)
                  {
                    if(arma::is_finite(raster(i_row_site_two,i_col_site_two)))
                    {
                      //calculate dissimilarity
                      arma::rowvec site_one_trans = raster_layers.subcube(arma::span(i_row_site_one),arma::span(i_col_site_one),arma::span());
                      arma::rowvec site_two_trans = raster_layers.subcube(arma::span(i_row_site_two),arma::span(i_col_site_two),arma::span());
                      arma::rowvec ecol_dist = abs(site_one_trans-site_two_trans);
                      for(int nvar=0; nvar<ncoef; nvar++)
                      {
                        prob(nvar) = ecol_dist(nvar)*coef(nvar);
                      }
                      double probs = sum(prob);
                      dissimilarity_one_two = 1/(1+exp(-probs));
                      // calculate site-pair weighting
                      mean_weight_one_two = (rel(i_site_one,2) + rel(i_site_two,2)) / 2;
                      // add weighted dissimilarity & weighting to the tally
                      d_dissimilarity_sum = d_dissimilarity_sum + (dissimilarity_one_two * mean_weight_one_two);
                      d_weighted_pair_sum = d_weighted_pair_sum + mean_weight_one_two;
                    } // end if !is.na(raster(i_row_site_two,i_col_site_two))
                  }// end if i_site_two valid
                } // end for i_site_two
              } // end if !is.na(raster(i_row_site_one,i_col_site_one))
            } // end if i_site_one valid
          } // end for i_site_one
          if(d_weighted_pair_sum < 1)
          {
            beta_raster(i_row,i_col) = 0;
          } // end if d_weighted_pair_sum <= 0
          else{
            beta_raster(i_row,i_col) = d_dissimilarity_sum/d_weighted_pair_sum;
          } // end else d_weighted_pair_sum <= 0
        } // end if raster[i_row,i_col] > -1
        //Rcpp::Rcout << i_row << " " << i_col << std::endl;
      } // end for i_col
    } // ## end for i_row
return(wrap(beta_raster));
}
