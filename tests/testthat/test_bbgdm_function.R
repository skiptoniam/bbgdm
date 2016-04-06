context("tests on bbgdm inputs")

test_that("tests for species data variable",{
  set.seed(12345)
  sp.dat <- matrix(rpois(200,1),20,10)
  env.dat <- simulate_covariates(sp.dat,2)
  form <- ~ 1 + covar_1 + covar_2
  expect_that(bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",nboot=10, 
                    scale_covar=F,geo=F,optim.meth='optim'),throws_error("Species data must be presence-absence!"))
  })

