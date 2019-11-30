---
title: Estimating $p_2$ models with `p2model`
author: Ruggero Bellio
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Here we present an example to illustrate how to estimate a $p_2$ model using
the `p2model` package.  The data used are the `high-tech managers` data, which were analyzed also in the  Example 1 of the paper illustrating 
the methodology (Bellio and Soriani, 2017). The data can be accessed
from the `NetData` package (Nowak et al., 2012).


```r
library(p2model)
library(NetData)
data(kracknets, package = "NetData")
```

The data are saved in the `friendship_data_frame`. The number of actors is $g=21$, and the object is a  data frame
with $21^2=441$ rows and three columns. The first five rows are as follows

```r
head(friendship_data_frame)
```

```
##   ego alter friendship_tie
## 1   1     1              0
## 2   1     2              1
## 3   1     3              0
## 4   1     4              1
## 5   1     5              0
## 6   1     6              0
```
The response values are contained in the third column.  The first step is to define a $21 \times 21$ adjacency matrix containing the response. This can be done with a simple loop

```r
g <- 21
Y <- matrix(0, g, g)
ind <- 1
for(i in 1:nrow(friendship_data_frame)){
  sele <- friendship_data_frame[i,]
  Y[sele$ego, sele$alter] <- sele$friendship_tie 	
  }
```
We then define two design matrixes for sender and receiver effects, that in this example are both equal to the $21 \times 2$
matrix containing the values of the two actor attributes `AGE` and 
`TENURE`. They  can be found in the `attributes` data frame,
with 21 rows and 4 columns. The required  matrix  can be created with the `model.matrix` command, excluding the first column (which is a 
column of ones)

```r
Xn <- model.matrix( ~ AGE + TENURE, attributes)[,-1]     
```
The next step is to create two 3-dimensional arrays for density and reciprocity effects. For density effects, we need a $21 \times 21 \times 4$ 
array. The first covariate is just a constant one. Then there are two 
dichotomized differences of sending and receiving actor covariate values for `DEPT` and `LEVEL`, and  the absolute difference of sending and receiving actor covariate values
for what concerns `AGE`

```r
XvD <- array(1, dim = c(g, g, 4))
for(i in 1:g)
 for(j in 1:g){ 
    XvD[i, j, 2] <- as.numeric(attributes$DEPT[i]==attributes$DEPT[j])
    XvD[i, j, 3] <- as.numeric(attributes$LEVEL[i]==attributes$LEVEL[j])
    XvD[i, j, 4] <- abs(attributes$AGE[i] - attributes$AGE[j])
   }
```
For reciprocity, there is only the intercept, and the related covariate array has just a  constant component

```r
XvC <- array(1, dim = c(g, g, 1))     
```
Now everything is in place for performing the actual estimation of the model,
that can be made by the `p2.fit` function, implementing approximate maximum
likelihood estimation based on the Laplace approximation

```r
fit <- fit_p2(Y, Xn, Xn, XvD, XvC)    
summary(fit)
```

```
## --------------------------------------------
## Approximate maximum likelihood estimation of p2 model
## Log-likelihood at maximum -169.4082 
## Model selection criteria
## AIC =  362.8164   BIC = 375.3507 
## --------------------------------------------
## Density coefficients
##           Estimate Std. error
## delta1  0.03885399 1.73521103
## delta1  1.57645521 0.34814083
## delta1  1.15354386 0.40545038
## delta1 -0.05488179 0.02458207
## --------------------------------------------
## Reciprocity coefficients
##        Estimate Std. error
## delta2 2.117182  0.6329099
## --------------------------------------------
## Sender coefficients
##          Estimate Std. error
## gamma1 -0.1309324 0.05911481
## gamma1  0.1403820 0.06213760
## --------------------------------------------
## Receiver coefficients
##            Estimate Std. error
## gamma2 -0.002099014 0.03503430
## gamma2  0.045075545 0.04146455
## --------------------------------------------
## Variance matrix of random effects
##           [,1]       [,2]
## [1,]  2.054007 -1.1362520
## [2,] -1.136252  0.9976877
## --------------------------------------------
```
The results are approximated maximum likelihood estimates based on the maximization of $L^*(\theta)$, i.e. 
the value of $\hat\theta^*$. 


We can plot the estimated sender and receiver effects for each node by
the `plot_effects_p2` function

```r
plot_effects_p2(fit)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

Another useful graphical representation is given by the goodness-of-fit plots, which
are obtained by the `gof_p2` function

```r
par(mfrow = c(2, 3))
g2 <- gof_p2(fit, GOF = ~idegree + odegree + distance + espartners + dspartners + triadcensus)
plot(g2, main="")
```

![plot of chunk gof plot](figure/gof plot-1.png)

The plots are obtained by exploiting  the possibility to simulate networks from a given model, 
that can be used directly by applying the 
`simulate_p2` function

```r
obj.sim <- simulate_p2(fit, nsim = 100)
summary(obj.sim)
```

```
##          Length Class        Mode
## networks 100    network.list list
```
Finally, the package makes available maximum likelihood estimation using 
parameter dependent Laplace Importance Sampling, implemented by the `fit_p2_IS`
function, though the results in Bellio and Soriani (2019) suggests that first-order
Laplace approximation usually provides satisfactory results. The method is computationally
demanding, but the `fit_p2_IS` function supports multiple-core parallel computation. The
basic way to call the function in this example would  simply  be `fit_p2_IS(fit, nc = parallel::detectCores())`.

