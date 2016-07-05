### BBGDM

[![Travis-CI Build Status](https://travis-ci.org/skiptoniam/bbgdm.svg?branch=master)](https://travis-ci.org/skiptoniam/bbgdm)

BBGDM is a R package for running Generalized Dissimilarity Models with Bayesian Bootstrap for parameter estimation. To install package run the following command in your R terminal

``` r
install.packages(c('devtools'))
devtools::install_github('skiptoniam/bbgdm')
```

##### Load the required libaries, we need vegan for the dune dataset.

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

##### Plot diagnostics

``` r
resids <- diagnostics(fm1)
par(mfrow=c(2,2))
plot(resids)
```

<img src="readme_files/figure-markdown_github/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

##### Plot response curves

``` r
response <- as.response(fm1)
par(mfrow=c(1,1))
plot(response)
```

<img src="readme_files/figure-markdown_github/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

##### Run 'Wald-like' test on parameters

``` r
library(xtable)
wt <- bbgdm.wald.test(fm1)
tab <- xtable(wt)
print(tab, type = "html")
```

<!-- html table generated in R 3.2.2 by xtable 1.8-2 package -->
<!-- Tue Jul 05 14:32:29 2016 -->
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
8.78
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
1.19
</td>
<td align="right">
3.00
</td>
<td align="right">
0.76
</td>
</tr>
</table>
