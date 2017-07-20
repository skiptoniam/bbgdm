context('response-class')

#generate some fake data
set.seed(12345)
sp.dat <- matrix(rbinom(200,1,.5),20,10)
env.dat <- simulate_covariates(sp.dat,2)
form <- ~ 1 + covar_1 + covar_2
fm <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
            nboot=10,geo=FALSE,optim.meth='optim')
env.dat$X <- runif(nrow(env.dat), min=-20, max=20)
env.dat$Y <- runif(nrow(env.dat), min=-20, max=20)

fm1 <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
            nboot=10,geo=TRUE,optim.meth='optim')

test_that('check response works', {

  # test basic calls
  responses <- as.response(fm)
  responses <- as.response(fm1)
  testthat::expect_true(is.response(responses))
  testthat::expect_false(is.response(c(1,1,1,1,4)))

  # but in wrong object.
  testthat::expect_error(as.response(dissim_table(sp.dat,env.dat,dism_metric="bray_curtis")))
  testthat::expect_error(as.response(matrix(NA,2,2)))

  testthat::expect_equal(class(responses$bspl), 'list')
  testthat::expect_equal(class(responses$bspl.05), 'list')
  testthat::expect_equal(class(responses$bspl.95), 'list')
  testthat::expect_equal(class(responses$X), 'list')
  testthat::expect_equal(class(responses$grid_real), 'data.frame')
  testthat::expect_equal(class(responses$diff_table), 'data.frame')

})

test_that('response plot works', {
  set.seed(12345)
  sp.dat <- matrix(rbinom(200,1,.5),20,10)
  env.dat <- simulate_covariates(sp.dat,2)
  form <- ~ 1 + covar_1 + covar_2
  fm <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
              nboot=10,geo=FALSE,optim.meth='optim')
  responses <- as.response(fm)
  img <- function() {
    plot(responses)
  }
  expect_identical(plot(responses), img())
  expect_equal(class(plot(responses)),'NULL')

})
