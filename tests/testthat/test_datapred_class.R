context('data.prep-class')

#generate some fake data
set.seed(12345)

sites <- paste0('site',sample(10,100,replace=TRUE))
sp <- paste0('sp',sample(20,100,replace=TRUE))
abn <- rpois(100,3)+1
sp.site.table <- data.frame(id=1:100,sites=sites,sp=sp,abn=abn)


test_that('check table to presence absence matrix works', {

  #test it works
  pa1 <- table2pam(sp.site.table, site.id = "sites", sp.id = "sp", abund = FALSE, abund.col = " No.of.specimens",siteXsp=TRUE)
  pN1 <- table2pam(sp.site.table, site.id = "sites", sp.id = "sp", abund = TRUE, abund.col = "abn",siteXsp=TRUE)

  pa2 <- table2pam(sp.site.table, site.id = "sites", sp.id = "sp", abund = FALSE, abund.col = " No.of.specimens",siteXsp=FALSE)
  pN2 <- table2pam(sp.site.table, site.id = "sites", sp.id = "sp", abund = TRUE, abund.col = "abn",siteXsp=FALSE)

  #throw some errors
  expect_error(table2pam(sp, site.id = "sites", sp.id = "sp", abund = FALSE, abund.col = " No.of.specimens",siteXsp=FALSE))
  expect_error(table2pam(sp.site.table, site.id = "sits", sp.id = "sp", abund = FALSE, abund.col = " No.of.specimens",siteXsp=FALSE))
  expect_error(table2pam(sp.site.table, site.id = "sites", sp.id = "s", abund = FALSE, abund.col = " No.of.specimens",siteXsp=FALSE))
  expect_error(table2pam(sp.site.table, site.id = "sites", sp.id = "sp", abund = TRUE, abund.col = " No.of.specimens",siteXsp=FALSE))
  expect_error(table2pam(sp.site.table, site.id = 1, sp.id = "sp", abund = TRUE, abund.col = " No.of.specimens",siteXsp=FALSE))
  expect_error(table2pam(sp.site.table, site.id = 1, sp.id = "sp", abund = TRUE, abund.col = sites,siteXsp=FALSE))

  testthat::expect_equal(class(pa1), 'matrix')

})

test_that('check pam2dissim works', {

  y <- table2pam(sp.site.table, site.id = "sites", sp.id = "sp")
  y1 <- table2pam(sp.site.table, site.id = "sites", sp.id = "sp", abund = TRUE, abund.col = "abn",siteXsp=TRUE)
  dism <-  pam2dissim(y,dism_metric = 'number_non_shared')
  dism1 <-  pam2dissim(y,dism_metric = 'bray_curtis')

  #can't deal with abundance yet.
  expect_error(dism <-  pam2dissim(c(1,2,2,34,3),dism_metric = 'number_non_shared'))
  expect_error(dism <-  pam2dissim(y1,dism_metric = 'bay_curtis'))
  expect_error(dism <-  pam2dissim(y1,dism_metric = 'number_non_shared'))
  expect_error(dism <-  pam2dissim(y1,dism_metric = 'bray_curtis'))

  expect_equal(class(dism),'list')
  expect_equal(class(dism$no_share),'matrix')
  expect_equal(class(dism$sum_share),'matrix')
  expect_equal(class(dism1),'matrix')

})

test_that('simulate works', {

  set.seed(12345)
  sp.dat <- matrix(rbinom(200,1,.5),20,10)
  env.dat <- simulate_covariates(sp.dat,2)
  expect_error(simulate_covariates(sp,2))
  expect_error(simulate_covariates(sp.dat,'a'))
  expect_error(simulate_covariates(c(1,2),2))

  expect_equal(class(env.dat),'data.frame')

})
