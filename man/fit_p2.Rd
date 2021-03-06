% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p2.R
\name{fit_p2}
\alias{fit_p2}
\title{Fitting function for p2 models}
\usage{
fit_p2(y, XnS, XnR, XvD, XvC, M = 0, seed = NULL, trace = FALSE,
  init = NULL, penalized = FALSE, penSigma = NULL, opt = nlminb,
  singular.ok = TRUE, ...)
}
\arguments{
\item{y}{An adjacency matrix of size \eqn{g \times g}{g x g}.}

\item{XnS}{Matrix of sender effects of size \eqn{g \times k_s}{g x ks}.}

\item{XnR}{Matrix of receiver effects of size \eqn{g \times k_r}{g x kr}.}

\item{XvD}{Density effects, a 3-dim array of size \eqn{g \times g \times
k_d}{g x g x kd}.}

\item{XvC}{Reciprocity effects, a 3-dim array of size \eqn{g \times g \times
k_c}{g x g x kc}.}

\item{M}{Number of replication for TMB-based importance sampling. Default is
0, giving the  Laplace approximation  (strongly recommended).}

\item{seed}{Random seed for importance sampling. Default is NULL.}

\item{trace}{TRUE for tracing information during the estimation. Default is FALSE.}

\item{init}{Optional starting value for model parameters. Default is NULL.}

\item{penalized}{Set TRUE for penalized estimation of variance parameters.
Default is FALSE.}

\item{penSigma}{Optional variance matrix of random effects to be used as
user-defined penalty in penalized estimation. Default is NULL.}

\item{opt}{Name of the optimising function. Default is nlminb.}

\item{singular.ok}{Should singular variance matrix of random effects be
allowed? Default is TRUE.}

\item{...}{Optional arguments passed to the optimiser.}
}
\value{
The returned value is an object of class
\code{"p2"}, a list  containing the following components:
\item{\code{theta}}{the vector of model estimates.}
\item{\code{loglik}}{the log likelihood function at the estimate.}
\item{\code{AIC, BIC}}{model selection criteria at the estimate.}
\item{\code{XnS.null, XnR.null}}{logical, flagging whether there
  are no sender or receiver effects, respectively.}
\item{\code{seed}}{the random seed used for estimation (when \code{M}>0).}
\item{\code{theta.cov}}{Variance matrix of estimates.}
\item{\code{theta.se}}{Standard errors of estimates.}
\item{\code{opt}}{Optimiser employed for estimation.}
\item{\code{opt.details}}{Object returned by the optimiser.}
\item{\code{ADobj}}{Object returned by \code{TMB:MakeADFun}, containing
the log likelihood function and its gradient (among its components).}
\item{\code{model.data}}{List containing  all the model data, fed
to \code{TMB:MakeADFun}.}
\item{\code{sdrep}}{Summary object returned by \code{TMB::sdreport}. This is an
object with many slots, useful for TMB users.}
\item{\code{ranef, ranef.se}}{Estimated random effects and their standard errors.}
\item{\code{Sigma, Sigma.se}}{Estimated variance matrix of random effects
and related standard errors.}
\item{\code{rho, rho.se}}{Estimated correlation of random effects
  and its standard error.}
\item{\code{M}}{Number of importance sampling replications.}
\item{\code{penflag, penSigma}}{Arguments for penalized estimation.}
}
\description{
Estimates a p2 model using the Laplace approximation.
}
\details{
The function allows for penalized estimation, which may be recommendable
to prevent numerical issues, such as estimated variance matrix of random effects
close to singularity. The fit is carried out by means of the TMB package,
and by selecting M>0 it would be possible to employ the importance sampler
made available by that package. However, the function \code{fitIS} provides
a safer alternative for estimation based on importance sampling.
}
\examples{
# Analysis of the  kracknets data from the NetData package
library(NetData)
data(kracknets)
# data preparation
g <- 21
Y <- matrix(0, g, g)
ind <-1
for(i in 1:nrow(friendship_data_frame)){
   sele <- friendship_data_frame[i, ]
   Y[sele$ego, sele$alter] <- sele$friendship_tie
   }
Xn <- model.matrix(~ AGE + TENURE, attributes)[, -1]
XvD <- array(1, dim=c(g, g, 4))
for(i in 1:g)
 for(j in 1:g){
    XvD[i, j, 2] <- as.numeric(attributes$DEPT[i]==attributes$DEPT[j])
    XvD[i, j, 3] <- as.numeric(attributes$LEVEL[i]==attributes$LEVEL[j])
    XvD[i, j, 4] <- abs(attributes$AGE[i] - attributes$AGE[j])
  }
XvC <- array(1, dim=c(g, g, 1))
# Now we are ready to fit the model
mod <- fit_p2(Y, Xn, Xn, XvD, XvC)
print(mod)
}
\references{
Bellio, R. and Soriani, N. (2019). Maximum likelihood estimation based on the
Laplace approximation for \eqn{p_2}{p2} network regression models. \emph{Submitted manuscript}.
}
\author{
Ruggero Bellio
}
