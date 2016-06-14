BBGDM
=====

[![Travis-CI Build Status](https://travis-ci.org/skiptoniam/bbgdm.svg?branch=master)](https://travis-ci.org/skiptoniam/bbgdm)

BBGDM is a R package for running Generalized Dissimilarity Models with Bayesian Bootstrap for parameter estimation. To install package run the following command in your R terminal

``` r
install.packages(c('devtools'))
devtools::install_github('skiptoniam/bbgdm')
```

### Load the required libaries, we need vegan for the dune dataset.

``` r
library(bbgdm)
library(vegan)
```

### Run bbgdm on the famous dune meadow data

The dune meadow vegetation data, dune, has cover class values of 30 species on 20 sites. Make the abundance data presence/absence.

``` r
data(dune)
data(dune.env)
dune_pa <- ifelse(dune>0,1,0)
```

### Fit a bbgdm

Now we have a species by sites matrix of simulated data and a set data for a one dimensional gradient.

``` r
form <- ~1+A1
fm1 <- bbgdm(form,dune_pa, dune.env,family="binomial",link='logit',
             dism_metric="number_non_shared",spline_type = 'ispline',
             nboot=100, geo=FALSE,optim.meth='nlmnib')
```

    ## binomial regression is on the way. 
    ## Calculating number_non_shared dissimilarity
    ## Transforming covariates to  ispline  with  2 degrees of freedom.
    ## Bayesian bootstrap  10  iterations
    ## Bayesian bootstrap  20  iterations
    ## Bayesian bootstrap  30  iterations
    ## Bayesian bootstrap  40  iterations
    ## Bayesian bootstrap  50  iterations
    ## Bayesian bootstrap  60  iterations
    ## Bayesian bootstrap  70  iterations
    ## Bayesian bootstrap  80  iterations
    ## Bayesian bootstrap  90  iterations
    ## Bayesian bootstrap  100  iterations

### Plot response curves

``` r
plotResponse(fm1,plotdim = c(1,1))
```

![](readme_files/figure-markdown_github/unnamed-chunk-5-1.png) \#\#\# Plot diagnostics

``` r
bbgdm.check(fm1)
```

![](readme_files/figure-markdown_github/unnamed-chunk-6-1.png) \#\#\# Run 'Wald-like' test on parameters

``` r
library(xtable)
```

    ## Warning: package 'xtable' was built under R version 3.2.5

``` r
wt <- bbgdm.wald.test(fm1)
xtable(wt, type = "html")
```

% latex table generated in R 3.2.2 by xtable 1.8-2 package % Tue Jun 14 11:10:12 2016
\begin{table}[ht]
\centering
\begin{tabular}{rrrr}
  \hline
 & bbgdm\_W & bbgdm\_df & bbgdm\_p-value \\ 
  \hline
intercept & 10.96 & 1.00 & 0.00 \\ 
  A1 & 1.32 & 3.00 & 0.73 \\ 
   \hline
\end{tabular}
\end{table}
