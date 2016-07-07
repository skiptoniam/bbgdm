context('spline-class')

#generate some fake data
set.seed(12345)
sp.dat <- matrix(rbinom(200,1,.5),20,10)
x <- simulate_covariates(sp.dat,2)

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

})
