# p2model
This package can be used to estimate the parameters of a p2 model for directed binary networks with correlated 
random effects. It implements (approximate) maximum likelihood estimation, following the methodology studied in 
Bellio and Soriani (2021, Statistica Neerlandica https://onlinelibrary.wiley.com/doi/10.1111/stan.12223). 
The internal engine for the estimation is based on the  `TMB` package. A vignette is available with some illustrative examples.

## Installation
``` 
#install.packages("remotes")
remotes::install_github("rugbel/p2model", build = TRUE)
```

