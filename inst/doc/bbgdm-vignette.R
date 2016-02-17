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
fm1 <- bbgdm(form,dune_pa, dune.env,family="binomial",dism_metric="number_shared",
             nboot=100, scale_covar=F,geo=F,optim.meth='nlmnib')

## ---- fig.show='hold' ,echo = TRUE, fig.height=6,fig.width=8-------------
plotResponse(fm1,plotdim = c(1,1))

## ---- fig.show='hold' ,echo = TRUE, fig.height=6,fig.width=8-------------
bbgdm.check(fm1)

## ------------------------------------------------------------------------
bbgdm.wald.test(fm1)

