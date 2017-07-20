## ------------------------------------------------------------------------
# install.packages(c('devtools','vegan'))
# devtools::install_github('skiptoniam/bbgdm')

## ----loadlibraries, echo = TRUE, error=FALSE, warning=FALSE, message=FALSE----
library(bbgdm)
library(vegan)

## ----dune data, error=FALSE, warning=FALSE, message=FALSE----------------
data(dune)
data(dune.env)
dune_pa <- ifelse(dune>0,1,0)

## ----fit model, results="hide", error=FALSE, warning=FALSE, message=FALSE----
form <- ~1+A1
fm1 <- bbgdm(form,dune_pa, dune.env,family="binomial",link='logit',
             dism_metric="number_non_shared",spline_type = 'ispline',
             nboot=100, geo=FALSE,optim.meth='nlmnib')

## ----fig.width = 6, fig.height = 6, fig.align='center'-------------------
resids <- diagnostics(fm1)
par(mfrow=c(2,2))
plot(resids)

## ----fig.width = 4, fig.height = 4, fig.align='center'-------------------
response <- as.response(fm1)
par(mfrow=c(1,1))
plot(response)

## ----results='asis'------------------------------------------------------
library(xtable)
wt <- bbgdm.wald.test(fm1)
tab <- xtable(wt)
print(tab, type = "html")

