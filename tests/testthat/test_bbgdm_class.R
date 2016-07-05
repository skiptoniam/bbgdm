context('bbgdm-class')

#generate some fake data
set.seed(12345)
sp.dat <- matrix(rbinom(200,1,.5),20,10)
env.dat <- simulate_covariates(sp.dat,2)
form <- ~ 1 + covar_1 + covar_2

test_that('bbgdm optimisers work', {

  # simple models with different optimisers
  fm_optim <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
               nboot=10,geo=FALSE,optim.meth='optim')
  fm1_nlmnib <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
               nboot=10,geo=FALSE,optim.meth = 'nlmnib')

})

test_that('bbgdm model fitting options work', {

  # simple models with different dissimilarity
  fm_optim_nns <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                    nboot=10,geo=FALSE,optim.meth='optim')
  fm1_optim_bc <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="bray_curtis",
                      nboot=10,geo=FALSE,optim.meth = 'optim')
  fm_nlmnib_nns <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                        nboot=10,geo=FALSE,optim.meth='nlmnib')
  fm1_nlmnib_bc <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="bray_curtis",
                        nboot=10,geo=FALSE,optim.meth = 'nlmnib')

  # isplines options
  fm_isp_df1_k1 <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                         nboot=10,geo=FALSE,optim.meth='optim',
                         spline_type = "ispline",spline_df = 1, spline_knots = 1)
  fm_isp_df2_k2 <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                         nboot=10,geo=FALSE,optim.meth='optim',
                         spline_type = "ispline",spline_df = 2, spline_knots = 2)
  fm_isp_df2_k3 <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                         nboot=10,geo=FALSE,optim.meth='optim',
                         spline_type = "ispline",spline_df = 2, spline_knots = 3)
  fm_isp_df3_k1 <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                        nboot=10,geo=FALSE,optim.meth='optim',
                        spline_type = "ispline",spline_df = 3, spline_knots = 1)
  fm_isp_df3_k3 <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                         nboot=10,geo=FALSE,optim.meth='optim',
                         spline_type = "ispline",spline_df = 3, spline_knots = 2)

  #b-splines now working.
  fm_bsp_df3_k1 <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                         nboot=10,geo=FALSE,optim.meth='optim',
                         spline_type = "bspline",spline_df = 2, spline_knots = 1)

  fm_bsp_df3_k3 <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                         nboot=10,geo=FALSE,optim.meth='optim',
                         spline_type = "bspline",spline_df = 3, spline_knots = 2)

  #check models when including geographic distance
  testthat::expect_error(fm_e <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                        nboot=10,geo=TRUE,optim.meth='optim'))
  testthat::expect_error(fm_e <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                                       nboot=10,geo=TRUE,geo.type = 'greater_circle',optim.meth='optim'))
  # add in coordinates
  env.dat$X <- runif(nrow(env.dat), min=-20, max=20)
  env.dat$Y <- runif(nrow(env.dat), min=-20, max=20)
  fm_e <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                                       nboot=10,geo=TRUE,optim.meth='optim')
  fm_gc <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                                       nboot=10,geo=TRUE,geo.type = 'greater_circle',optim.meth='optim')


})
