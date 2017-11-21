## ---- read data, message=FALSE-------------------------------------------
library(p2model)
library(NetData)
data(kracknets, package = "NetData")

## ---- display data-------------------------------------------------------
head(friendship_data_frame)

## ---- define adjacency matrix--------------------------------------------
g <- 21
Y <- matrix(0, g, g)
ind <- 1
for(i in 1:nrow(friendship_data_frame)){
  sele <- friendship_data_frame[i,]
  Y[sele$ego, sele$alter] <- sele$friendship_tie 	
  }

## ---- Define the model matrices I----------------------------------------
Xn <- model.matrix( ~ AGE + TENURE, attributes)[,-1]     

## ---- Define the model matrices II---------------------------------------
XvD <- array(1, dim = c(g, g, 4))
for(i in 1:g)
 for(j in 1:g){ 
    XvD[i, j, 2] <- as.numeric(attributes$DEPT[i]==attributes$DEPT[j])
    XvD[i, j, 3] <- as.numeric(attributes$LEVEL[i]==attributes$LEVEL[j])
    XvD[i, j, 4] <- abs(attributes$AGE[i] - attributes$AGE[j])
   }

## ---- Define the model matrices, III-------------------------------------
XvC <- array(1, dim = c(g, g, 1))     

## ---- Fitting the model--------------------------------------------------
fit <- fit_p2(Y, Xn, Xn, XvD, XvC)    

## ---- fig.show='hold', fig.width=7.5, fig.height=4.5---------------------
plot_effects_p2(fit)

## ---- gof plot, fig.show='hold', fig.width=7.5, fig.height=4.5-----------
par(mfrow = c(2, 3))
g2 <- gof_p2(fit, GOF = ~idegree + odegree + distance + espartners + dspartners + triadcensus)
plot(g2, main="")

## ---- simulation---------------------------------------------------------
obj.sim <- simulate_p2(fit, nsim = 100)
summary(obj.sim)

