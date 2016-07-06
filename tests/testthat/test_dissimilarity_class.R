context('dissimilarity-class')

#generate some fake data
set.seed(12345)
sp.dat <- matrix(rbinom(200,1,.5),20,10)
env.dat <- simulate_covariates(sp.dat,2)

test_that('check diff_table works', {

  # test dissimilarities
  diff_table_nns <- dissim_table(sp.dat,env.dat,dism_metric="number_non_shared")
  diff_table_bc <- dissim_table(sp.dat,env.dat,dism_metric="bray_curtis")
  testthat::expect_true(is.dissim_table(diff_table_nns))
  testthat::expect_true(is.dissim_table(diff_table_bc))

  # miss type dissimilarity name - expect error.
  testthat::expect_error(dissim_table(sp.dat,env.dat,dism_metric="bray_curti"))
  testthat::expect_error(dissim_table(sp.dat,env.dat,dism_metric="number_no_shared"))

  #expect error because no coordinates
  testthat::expect_error(diff_table_nns <- dissim_table(sp.dat,env.dat,geo = TRUE, dism_metric="number_non_shared"))

  #add coords - all ok.
  env.dat$X <- runif(nrow(env.dat), min=-20, max=20)
  env.dat$Y <- runif(nrow(env.dat), min=-20, max=20)
  diff_table_nns <- dissim_table(sp.dat,env.dat,geo = TRUE, dism_metric="number_non_shared")

  #test euclidean
  diff_table_nns_e_d <- dissim_table(sp.dat,env.dat,geo = TRUE,geo.type = 'euclidean', dism_metric="number_non_shared")

  #test gc
  diff_table_nns_gc_d <- dissim_table(sp.dat,env.dat,geo = TRUE,geo.type = 'greater_circle', dism_metric="number_non_shared")

  #Miss type the distance metric
  testthat::expect_error(dissim_table(sp.dat,env.dat,geo = TRUE,geo.type = 'eucliean', dism_metric="bray_curtis"))
  testthat::expect_error(dissim_table(sp.dat,env.dat,geo = TRUE,geo.type = 'greatercircle',dism_metric="number_non_shared"))

  ## Check splines work
  #bsplines
  diff_table_nns_e_bs <- dissim_table(sp.dat,env.dat,geo = TRUE,geo.type = 'euclidean',
                                      spline_type = 'bspline', dism_metric="number_non_shared")
  diff_table_nns_gc_bs <- dissim_table(sp.dat,env.dat,geo = TRUE,geo.type = 'greater_circle',
                                       spline_type = 'bspline', dism_metric="number_non_shared")
  #isplines
  diff_table_nns_e_is <- dissim_table(sp.dat,env.dat,geo = TRUE,geo.type = 'euclidean',
                                     spline_type = 'ispline', dism_metric="number_non_shared")
  diff_table_nns_gc_is <- dissim_table(sp.dat,env.dat,geo = TRUE,geo.type = 'greater_circle',
                                      spline_type = 'ispline', dism_metric="number_non_shared")

  testthat::expect_equal(class(diff_table_nns_gc_is$diff_table), 'matrix')
  testthat::expect_equal(class(diff_table_nns_gc_is$diff_table_params), 'list')
  testthat::expect_equal(dim(diff_table_nns_gc_is$diff_table), c(190, 11))
})

