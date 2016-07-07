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
  expect_error(fm1_nlmnib <- bbgdm(form,sp.dat,env.dat,family="binol",dism_metric="number_non_shared",
                                     nboot=10,geo=FALSE,optim.meth = 'nlmnib'))

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


test_that('model prediction works', {

  #generate some fake data
  library(raster)
  set.seed(123)
  xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
  d <- as.matrix(dist(xy))
  w <- exp(-1/nrow(xy) * d)
  ww <- chol(w)
  xy$z <- t(ww) %*% rnorm(nrow(xy), 0, 0.1)
  xy$z <- scales::rescale(xy$z,range(env.dat$covar_1))
  coordinates(xy) <- ~x+y
  r <- rasterize(xy, raster(points2grid(xy)), 'z')
  #give it the same name as variable in bbgdm model.
  names(r)<- 'covar_1'
  r2 <- raster(r)
  res(r2) <- 0.05
  r2 <- resample(r, r2)
  r3 <- r2
  names(r3)<- 'covar_2'

  set.seed(12345)
  sp.dat <- matrix(rbinom(200,1,.5),20,10)
  env.dat <- simulate_covariates(sp.dat,2)
  form <- ~ 1 + covar_1 + covar_2
  fm <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
                    nboot=10,geo=FALSE,optim.meth='optim')

  #expect error is name of layer doesn't match covariates in model.
  r4<- r3
  names(r4) <- 'foo'
  testthat::expect_error(pred.fm.ut <- predict(fm,stack(r2,r4),uncertainty = TRUE))

  #use these layer to predict turnover.
  pred.fm.ut <- predict(fm,stack(r2,r3),uncertainty = TRUE)
  pred.fm.uf <- predict(fm,stack(r2,r3),uncertainty = FALSE)

  #expect error if dims miss match
  expect_error(pred.fm.ut <- predict(fm,r2,uncertainty = TRUE))

  # if wrong prediction data is submitted - should be raster, stack or brick.
  testthat::expect_error(pred.fm.ut <- predict(fm,as.matrix(stack(r2,r3)),uncertainty = TRUE))
  testthat::expect_error(pred.fm.uf <- predict(fm,c(stack(r2,r3)),uncertainty = FALSE))

  # if wrong neighbourhood is submitted - should be an int.
  testthat::expect_error(pred.fm.ut <- predict(fm,stack(r2,r3),neighbourhood=as.character(1),uncertainty = TRUE))
  testthat::expect_error(pred.fm.uf <- predict(fm,stack(r2,r3),neighbourhood=c(3,2,1),uncertainty = FALSE))

  # check classes
  expect_equal(class(pred.fm.uf)[1], 'RasterLayer')
  expect_equal(class(pred.fm.ut), 'list')
  expect_equal(class(pred.fm.ut[[1]])[1], 'RasterLayer')
  expect_equal(class(pred.fm.ut[[2]])[1], 'RasterLayer')

  # check dimensions of prediction match input.
  expect_equal(dim(pred.fm.uf), dim(r2))

})


test_that('model plot works', {

   fm <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
              nboot=10,geo=FALSE,optim.meth='optim')
   img <- function() {
    plot(fm)
  }
   expect_identical(plot(fm), img())
   expect_equal(class(plot(fm)),'NULL')

})

test_that('model print works', {
  set.seed(12345)
  sp.dat <- matrix(rbinom(200,1,.5),20,10)
  env.dat <- simulate_covariates(sp.dat,2)
  form <- ~ 1 + covar_1 + covar_2
  fm <- bbgdm(form,sp.dat,env.dat,family="binomial",dism_metric="number_non_shared",
              nboot=10,geo=FALSE,optim.meth='optim')
  # check print.dynamic works
  expect_equal(class(print(fm)),'NULL')

  expected <- c(" A Bayesian Bootstrap GDM fitted against:",
              " 20 sites,"," 10 species and ",
              " 190 dissimilarities used as observations in the model.","",
              " A total of 10 Bayesian Bootstraps were run.")

  expect_equal(capture.output(print(fm))[1:6], expected)

})
