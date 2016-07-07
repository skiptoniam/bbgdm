context('spline-class')

#generate some fake data
set.seed(12345)
sp.dat <- matrix(rbinom(200,1,.5),20,10)
x <- simulate_covariates(sp.dat,2)
form <- ~ 1 + covar_1 + covar_2
fm <- bbgdm(form,sp.dat,x,family="binomial",dism_metric="number_non_shared",
            nboot=10,geo=FALSE,optim.meth='optim')
fm_bs <- bbgdm(form,sp.dat,x,family="binomial",dism_metric="number_non_shared",
               spline_type = 'bspline', nboot=10,geo=FALSE,optim.meth='optim')

test_that('check ispline works', {

  x1 <- ispline(x[,1], spline.knots = 2, knots = NULL, spline.degree = 3)
  x2 <- ispline(x[,1], spline.knots = 1, knots = NULL, spline.degree = 3)
  x3 <- ispline(x[,1], spline.knots = 1, knots = 1, spline.degree = 1)

  #worng inputs
  expect_error(x3 <- ispline(x[,1], spline.knots = 'a', knots = NULL, spline.degree = 3))
  expect_error(x3 <- ispline(x[,1], spline.knots = 1, knots = x, spline.degree = 3))
  expect_error(x3 <- ispline(x[,1], spline.knots = 1, knots = NULL, spline.degree = letters[seq(1,6)]))

  #error if matrix
  testthat::expect_error(x <- ispline(x, spline.knots = 2, knots = NULL, spline.degree = 3))
  testthat::expect_error(x <- ispline(c("a","b"), spline.knots = 2, knots = NULL, spline.degree = 3))

  #test outputs
  testthat::expect_equal(class(x1), 'matrix')
  testthat::expect_equal(class(attributes(x1)), 'list')

})


test_that( "spline.trans works", {

  x1 <- spline.trans(x, spline_type = "ispline", spline_df = 2, spline_knots = 1)
  x1 <- spline.trans(x, spline_type = "ispline", spline_df = 2, spline_knots = NULL)
  x1 <- spline.trans(x, spline_type = "bspline", spline_df = 2, spline_knots = 1)
  x1 <- spline.trans(x, spline_type = "ispline", spline_df = 4, spline_knots = 1)
  x1 <- spline.trans(x, spline_type = "bspline", spline_df = 2, spline_knots = 2)

  # #worng inputs
  expect_error(x3 <- spline.trans(x, spline_type = "ipline", spline_df = 2, spline_knots = 1))
  expect_error(x3 <- spline.trans(x, spline_type = "ispline", spline_df = 'a', spline_knots = 1))
  expect_error(x3 <- spline.trans(x, spline_type = "ispline", spline_df = 2, spline_knots = 'a'))

  #test outputs
  expect_true(is.spline(x1))
  testthat::expect_equal(class(x1), 'spline')
  testthat::expect_equal(class(attributes(x1)), 'list')
  testthat::expect_equal(class(x1$spline), 'matrix')
  testthat::expect_equal(class(x1$spline.attr), 'list')

})

test_that( "Tests for spline_trans_for_pred", {

  object <- fm
  Xold <- data.matrix(object$dissim_dat)
  splineLength <- sapply(object$dissim_dat_params, `[[`, "dim")[2,]
  betas <- object$median.coefs.se[2:length(object$median.coefs.se)]
  betas.quantiles <- object$quantiles.coefs.se[,2:length(object$starting_gdm$coef),drop=FALSE]
  k <- ncol(object$env.dat)
  nr_df<-((nrow(object$sp.dat)^2)-nrow(object$sp.dat))/2
  nc_dt<-ncol(object$env.dat)
  ne<-ncol(object$env.dat)
  diff_table <- diff_table_cpp(as.matrix(object$env.dat))
  colnames(diff_table) <-c(colnames(object$env.dat))
  grid <- matrix(rep(seq(0, 1, length.out = 100), k), ncol = k)
  grid <- t(t(grid) * as.vector(diff(apply(diff_table, 2, range))) + apply(diff_table, 2, min))
  grid <- data.frame(grid)
  min_env <- apply(object$env.dat,2,function(x)min(abs(x)))
  max_env <- apply(object$env.dat,2,function(x)max(abs(x)))
  grid_real <- grid
  for(i in 1:ncol(grid)) grid_real[,i] <-  scales::rescale(grid[,i],to=c(min_env[i],max_env[i]))
  X <- mapply(spline_trans_for_pred, grid, attrib = object$dissim_dat_params,
              SIMPLIFY = FALSE)
  sp_t <- spline_trans_for_pred(grid[,1], attrib = object$dissim_dat_params[[1]], splineInterval = NULL, splineDegree = NULL)
  sp_t <- spline_trans_for_pred(grid[,1], attrib = object$dissim_dat_params[[1]],  splineInterval = c(1,2,3), splineDegree = NULL)
  expect_error(sp_t <- spline_trans_for_pred(grid, attrib = object$dissim_dat_params[[1]], splineInterval = NULL, splineDegree = NULL))
  expect_error(sp_t <- spline_trans_for_pred(grid, attrib = object$dissim_dat_params, splineInterval = NULL, splineDegree = NULL))
  #test outputs
  testthat::expect_equal(class(sp_t), 'matrix')

  object <- fm_bs
  Xold <- data.matrix(object$dissim_dat)
  splineLength <- sapply(object$dissim_dat_params, `[[`, "dim")[2,]
  betas <- object$median.coefs.se[2:length(object$median.coefs.se)]
  betas.quantiles <- object$quantiles.coefs.se[,2:length(object$starting_gdm$coef),drop=FALSE]
  k <- ncol(object$env.dat)
  nr_df<-((nrow(object$sp.dat)^2)-nrow(object$sp.dat))/2
  nc_dt<-ncol(object$env.dat)
  ne<-ncol(object$env.dat)
  diff_table <- diff_table_cpp(as.matrix(object$env.dat))
  colnames(diff_table) <-c(colnames(object$env.dat))
  grid <- matrix(rep(seq(0, 1, length.out = 100), k), ncol = k)
  grid <- t(t(grid) * as.vector(diff(apply(diff_table, 2, range))) + apply(diff_table, 2, min))
  grid <- data.frame(grid)
  min_env <- apply(object$env.dat,2,function(x)min(abs(x)))
  max_env <- apply(object$env.dat,2,function(x)max(abs(x)))
  grid_real <- grid
  for(i in 1:ncol(grid)) grid_real[,i] <-  scales::rescale(grid[,i],to=c(min_env[i],max_env[i]))
  X <- mapply(spline_trans_for_pred, grid, attrib = object$dissim_dat_params,
              SIMPLIFY = FALSE)
  attrib <- object$dissim_dat_params[[1]]
  sp_t <- spline_trans_for_pred(grid[,1], attrib = object$dissim_dat_params[[1]], splineInterval = NULL, splineDegree = NULL)
  expect_error(sp_t <- spline_trans_for_pred(grid[,1], attrib = NULL, splineInterval = NULL, splineDegree = NULL))
  sp_t <- spline_trans_for_pred(grid[,1], attrib = object$dissim_dat_params[[1]],splineInterval = attrib$knots, splineDegree = NULL)
  sp_t <- spline_trans_for_pred(grid[,1], attrib = object$dissim_dat_params[[1]], splineInterval = attrib$knots, splineDegree = attrib$degree)
  expect_error(sp_t <- spline_trans_for_pred(grid, attrib = object$dissim_dat_params[[1]], splineInterval = NULL, splineDegree = NULL))
  expect_error(sp_t <- spline_trans_for_pred(grid, attrib = object$dissim_dat_params, splineInterval = NULL, splineDegree = NULL))
  #test outputs
  testthat::expect_equal(class(sp_t)[1], 'bs')
  testthat::expect_equal(class(sp_t)[2], 'basis')
  testthat::expect_equal(class(sp_t)[3], 'matrix')

})
