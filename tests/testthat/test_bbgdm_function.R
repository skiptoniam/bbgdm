context("tests on bbgdm inputs")

test_that("tests for species data variable",{
  set.seed(12345)
  sp.dat <- matrix(rbinom(200,1,.5),20,10)
  env.dat <- simulate_covariates(sp.dat,2)
  form <- ~ 1 + covar_1 + covar_2
  fm1 <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
               nboot=10,geo=F,optim.meth='optim')
  testthat::expect_is(fm1,'bbgdm')
})

