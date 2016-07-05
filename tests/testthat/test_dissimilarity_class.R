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

  # diff_table_nns <- dissim_table(sp.dat,env.dat,geo = TRUE, dism_metric="number_non_shared")

})

