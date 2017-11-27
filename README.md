# p2model
This package can be used to estimate the parameters of a p2 model for directed binary networks with correlated 
random effects. It implements (approximate) maximum likelihood estimation, following the methodology studied in 
Bellio and Soriani (2017, _Submitted Manuscript_). The internal engine for the estimation is based on the 
`TMB` package. A vignette is available with some illustrative examples.

## Installation
``` 
#install.packages("devtools")
devtools::install_github(“rugbel/p2model", build_vignettes = TRUE)
```
