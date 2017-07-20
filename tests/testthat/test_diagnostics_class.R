context('diagnostics-class')

#generate some fake data
set.seed(12345)
sp.dat <- matrix(rbinom(200,1,.5),20,10)
env.dat <- simulate_covariates(sp.dat,2)
form <- ~ 1 + covar_1 + covar_2
fm <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                    nboot=10,geo=FALSE,optim.meth='optim')
env.dat$X <- runif(nrow(env.dat), min=-20, max=20)
env.dat$Y <- runif(nrow(env.dat), min=-20, max=20)
fm1 <- bbgdm(form,sp.dat,env.dat,family="binomial",link='negexp',dism_metric="number_non_shared",
            nboot=10,geo=TRUE,optim.meth='optim')

test_that('check diagnostics works', {

  # test basic calls
  resids <- diagnostics(fm)
  testthat::expect_true(is.diagnostics(resids))
  resids1 <- diagnostics(fm1)
  testthat::expect_true(is.diagnostics(resids1))
  testthat::expect_false(is.diagnostics(c(1,1,1,1,4)))

  # but in wrong object.
  testthat::expect_error(diagnostics(dissim_table(sp.dat,env.dat,dism_metric="bray_curtis")))
  testthat::expect_error(diagnostics(matrix(NA,2,2)))

  testthat::expect_equal(class(resids$res), 'numeric')
  testthat::expect_equal(class(resids$pi), 'numeric')
  testthat::expect_equal(class(resids$y), 'matrix')

})

test_that('check bbgdm.wald.test works', {

  # test basic calls
  wt <- bbgdm.wald.test(fm)
  wt1 <- bbgdm.wald.test(fm,H0 = 2)
  wt2 <- bbgdm.wald.test(fm,H0 = 0,gdm = TRUE)

  testthat::expect_false(is.diagnostics(wt))

  wt <- bbgdm.wald.test(fm1)
  wt1 <- bbgdm.wald.test(fm1,H0 = 2)
  wt2 <- bbgdm.wald.test(fm1,H0 = 0,gdm = TRUE)

  testthat::expect_false(is.diagnostics(wt1))
  # but in wrong object.
  testthat::expect_error(bbgdm.wald.test(c(0,1,1)))
  testthat::expect_error(bbgdm.wald.test(matrix(NA,2,2)))
  testthat::expect_error(bbgdm.wald.test(fm,H0 = 'a'))
  testthat::expect_error(bbgdm.wald.test(fm,H0 = c(1,2)))
  testthat::expect_error(bbgdm.wald.test(fm,gdm = TRU))

  testthat::expect_equal(class(wt), 'matrix')

})

test_that('residual plot works', {

  resids <- diagnostics(fm)
  img <- function() {
    plot(resids)
  }
  expect_identical(plot(resids), img())
  expect_equal(class(plot(resids)),'NULL')

})
