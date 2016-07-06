### BBGDM

[![Travis-CI Build Status](https://travis-ci.org/skiptoniam/bbgdm.svg?branch=master)](https://travis-ci.org/skiptoniam/bbgdm)

BBGDM is a R package for running Generalized Dissimilarity Models with Bayesian Bootstrap for parameter estimation. To install package run the following command in your R terminal

``` r
install.packages(c('devtools'))
devtools::install_github('skiptoniam/bbgdm')
```

##### Load the required libraries, we need vegan for the dune dataset.

``` r
library(bbgdm)
library(vegan)
```

##### Run bbgdm on the famous dune meadow data

The dune meadow vegetation data, dune, has cover class values of 30 species on 20 sites. Make the abundance data presence/absence.

``` r
data(dune)
data(dune.env)
dune_pa <- ifelse(dune>0,1,0)
```

##### Fit a bbgdm

Now we have a species by sites matrix of simulated data and a set data for a one dimensional gradient.

``` r
form <- ~1+A1
fm1 <- bbgdm(form,dune_pa, dune.env,family="binomial",link='logit',
             dism_metric="number_non_shared",spline_type = 'ispline',
             nboot=100, geo=FALSE,optim.meth='nlmnib')
```

##### Print model summary

``` r
print(fm1)
```

    ##  A Bayesian Bootstrap GDM fitted against:
    ##  20 sites,
    ##  30 species and 
    ##  190 dissimilarities used as observations in the model.
    ## 
    ##  A total of 100 Bayesian Bootstraps were run.
    ## 
    ##  Spline base parameter estimates are: 
    ##  (Intercept) 0.6535
    ##  x_1 0
    ##  x_2 0.2768
    ##  x_3 0.5422

##### Plot diagnostics

``` r
resids <- diagnostics(fm1)
par(mfrow=c(2,2))
plot(resids)
```

<img src="readme_files/figure-markdown_github/unnamed-chunk-6-1.png" title="" alt="" style="display: block; margin: auto;" />

##### Plot response curves

``` r
response <- as.response(fm1)
par(mfrow=c(1,1))
plot(response)
```

<img src="readme_files/figure-markdown_github/unnamed-chunk-7-1.png" title="" alt="" style="display: block; margin: auto;" />

##### Run 'Wald-like' test on parameters

``` r
library(xtable)
wt <- bbgdm.wald.test(fm1)
tab <- xtable(wt)
print(tab, type = "html")
```

<!-- html table generated in R 3.2.2 by xtable 1.8-0 package -->
<!-- Thu Jul 07 01:29:08 2016 -->
<table border="1">
<tr>
<th>
</th>
<th>
bbgdm\_W
</th>
<th>
bbgdm\_df
</th>
<th>
bbgdm\_p-value
</th>
</tr>
<tr>
<td align="right">
intercept
</td>
<td align="right">
12.52
</td>
<td align="right">
1.00
</td>
<td align="right">
0.00
</td>
</tr>
<tr>
<td align="right">
A1
</td>
<td align="right">
3.77
</td>
<td align="right">
3.00
</td>
<td align="right">
0.29
</td>
</tr>
</table>
##### Predict BBGDM model

``` r
#generate some random spatial autocorrelated data.
library(raster)
set.seed(123)
xy <- expand.grid(x=seq(145, 150, 0.1), y=seq(-40, -35, 0.1))
d <- as.matrix(dist(xy))
w <- exp(-1/nrow(xy) * d)
ww <- chol(w)
xy$z <- t(ww) %*% rnorm(nrow(xy), 0, 0.1)
xy$z <- scales::rescale(xy$z,range(dune.env$A1))
coordinates(xy) <- ~x+y
r <- rasterize(xy, raster(points2grid(xy)), 'z')
#give it the same name as variable in bbgdm model.
names(r)<- 'A1'
r2 <- raster(r)
res(r2) <- 0.05
r2 <- resample(r, r2)
#use this layer to predict turnover.
pred.dune.sim.dat <- predict(fm1,r2,uncertainty = TRUE)
```

    ## using default three cell neighbourhood to estimate dissimilarity

``` r
#plot the data and the turnover predictions, and error.
colram <- colorRampPalette(c("darkblue","yellow","red"))
colram.se <- colorRampPalette(c('antiquewhite','pink','red'))
par(mfrow=c(1,3),mar=c(3,2,2,6))
plot(r2,main='Simulated A1')
plot(pred.dune.sim.dat[[1]],col=colram(100),main='BBGDM turnover')
plot(pred.dune.sim.dat[[2]],col=colram.se(100),main='BBGDM CV of turnover')
```

<img src="readme_files/figure-markdown_github/unnamed-chunk-9-1.png" title="" alt="" style="display: block; margin: auto;" />
