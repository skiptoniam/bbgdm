context('bbgdm.fit-class')

#generate some fake data and transform into model.matrix format for bbgdm.fit
set.seed(12345)
sp.dat <- matrix(rbinom(200,1,.5),20,10)
env.dat <- simulate_covariates(sp.dat,2)
form <- ~ 1 + covar_1 + covar_2
left <- "cbind(nonsharedspp_ij,sumspp_ij)"
env.dat <- model.frame(form, as.data.frame(env.dat))
env.dat <- model.frame(as.data.frame(env.dat))
dissim_dat <- dissim_table(sp.dat,env.dat,dism_metric='number_non_shared',spline_type='ispline')
dissim_dat_table <- as.data.frame(dissim_dat$diff_table)
dissim_dat_params <- dissim_dat$diff_table_params
preds <- colnames(dissim_dat_table[,3:ncol(dissim_dat_table)])
form <- update.formula(form, paste0(left," ~ ",paste(preds, collapse=" + ")))
temp <- model.frame(form, as.data.frame(dissim_dat_table))
y <- model.response(temp)
X <- model.matrix(form, as.data.frame(dissim_dat_table))
Nsite <- nrow(env.dat)
w <- gtools::rdirichlet(Nsite, rep(1/Nsite,Nsite))
wij <- w%*%t(w)
wij <- wij[upper.tri(wij)]
link <- 'logit'

test_that('bbgdm fit works', {
  # simple models with different optimisers
  mod  <- bbgdm.fit(X,y,link=link)
})

test_that('bbgdm.fit options work', {

  # simple fit with different link functions
  mod1  <- bbgdm.fit(X,y,link='logit')
  mod2  <- bbgdm.fit(X,y,link='probit')
  mod3  <- bbgdm.fit(X,y,link='cloglog')
  mod4  <- bbgdm.fit(X,y,link='negexp')

  # simple fit with different optim methods
  mod5  <- bbgdm.fit(X,y,link='logit',optim.meth = 'optim')
  mod6  <- bbgdm.fit(X,y,link='logit',optim.meth = 'nlmnib')

  #add some weights
  mod7  <- bbgdm.fit(X,y,wt=wij,link='logit',optim.meth = 'optim')

  #add erroneous weights
  expect_error(mod8  <- bbgdm.fit(X,y,wt=wij[1:30],link='logit',optim.meth = 'optim'))
  expect_error(mod9  <- bbgdm.fit(X,y,wt=as.character(wij),link='logit',optim.meth = 'optim'))

  # check classes
  expect_equal(class(mod7), 'gdm')
  expect_equal(class(mod7$y), 'matrix')
  expect_equal(class(mod7$X), 'matrix')
  expect_equal(class(mod7$var), 'matrix')

})

