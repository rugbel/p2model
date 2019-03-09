# display version number and date when the package is loaded
#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
    'Package:  p2model\n',
    'Version:  ', desc$Version, '\n',
    'Date:     ', desc$Date, '\n',
    'Authors:  Ruggero Bellio (University of Udine)\n')
  }


#' Fitting function for p2 models
#'
#' Estimates a p2 model using the Laplace approximation.
#'
#' The function allows for penalized estimation, which may be recommendable
#' to prevent numerical issues, such as estimated variance matrix of random effects
#' close to singularity. The fit is carried out by means of the TMB package,
#' and by selecting M>0 it would be possible to employ the importance sampler
#' made available by that package. However, the function \code{fitIS} provides
#' a safer alternative for estimation based on importance sampling.
#'
#' @param y An adjacency matrix of size \eqn{g \times g}{g x g}.
#' @param XnS Matrix of sender effects of size \eqn{g \times k_s}{g x ks}.
#' @param XnR Matrix of receiver effects of size \eqn{g \times k_r}{g x kr}.
#' @param XvD Density effects, a 3-dim array of size \eqn{g \times g \times
#'   k_d}{g x g x kd}.
#' @param XvC Reciprocity effects, a 3-dim array of size \eqn{g \times g \times
#'   k_c}{g x g x kc}.
#' @param M Number of replication for TMB-based importance sampling. Default is
#'   0, giving the  Laplace approximation  (strongly recommended).
#' @param seed Random seed for importance sampling. Default is NULL.
#' @param trace TRUE for tracing information during the estimation. Default is FALSE.
#' @param init Optional starting value for model parameters. Default is NULL.
#' @param penalized Set TRUE for penalized estimation of variance parameters.
#'   Default is FALSE.
#' @param penSigma Optional variance matrix of random effects to be used as
#'   user-defined penalty in penalized estimation. Default is NULL.
#' @param opt Name of the optimising function. Default is nlminb.
#' @param singular.ok Should singular variance matrix of random effects be
#'   allowed? Default is TRUE.
#' @param ... Optional arguments passed to the optimiser.
#' @return The returned value is an object of class
#' \code{"p2"}, a list  containing the following components:
#' \item{\code{theta}}{the vector of model estimates.}
#' \item{\code{loglik}}{the log likelihood function at the estimate.}
#' \item{\code{AIC, BIC}}{model selection criteria at the estimate.}
#' \item{\code{XnS.null, XnR.null}}{logical, flagging whether there
#'   are no sender or receiver effects, respectively.}
#' \item{\code{seed}}{the random seed used for estimation (when \code{M}>0).}
#' \item{\code{theta.cov}}{Variance matrix of estimates.}
#' \item{\code{theta.se}}{Standard errors of estimates.}
#' \item{\code{opt}}{Optimiser employed for estimation.}
#' \item{\code{opt.details}}{Object returned by the optimiser.}
#' \item{\code{ADobj}}{Object returned by \code{TMB:MakeADFun}, containing
#' the log likelihood function and its gradient (among its components).}
#' \item{\code{model.data}}{List containing  all the model data, fed
#' to \code{TMB:MakeADFun}.}
#' \item{\code{sdrep}}{Summary object returned by \code{TMB::sdreport}. This is an
#' object with many slots, useful for TMB users.}
#' \item{\code{ranef, ranef.se}}{Estimated random effects and their standard errors.}
#' \item{\code{Sigma, Sigma.se}}{Estimated variance matrix of random effects
#' and related standard errors.}
#' \item{\code{rho, rho.se}}{Estimated correlation of random effects
#'   and its standard error.}
#' \item{\code{M}}{Number of importance sampling replications.}
#' \item{\code{penflag, penSigma}}{Arguments for penalized estimation.}
#' @export
#' @import NetData
#' @useDynLib p2model
#' @author Ruggero Bellio
#' @references  Bellio, R. and Soriani, N. (2019). Maximum likelihood estimation based on the
#' Laplace approximation for \eqn{p_2}{p2} network regression models. \emph{Submitted manuscript}.
#' @examples
#' # Analysis of the  kracknets data from the NetData package
#' library(NetData)
#' data(kracknets)
#' # data preparation
#' g <- 21
#' Y <- matrix(0, g, g)
#' ind <-1
#' for(i in 1:nrow(friendship_data_frame)){
#'    sele <- friendship_data_frame[i, ]
#'    Y[sele$ego, sele$alter] <- sele$friendship_tie
#'    }
#' Xn <- model.matrix(~ AGE + TENURE, attributes)[, -1]
#' XvD <- array(1, dim=c(g, g, 4))
#' for(i in 1:g)
#'  for(j in 1:g){
#'     XvD[i, j, 2] <- as.numeric(attributes$DEPT[i]==attributes$DEPT[j])
#'     XvD[i, j, 3] <- as.numeric(attributes$LEVEL[i]==attributes$LEVEL[j])
#'     XvD[i, j, 4] <- abs(attributes$AGE[i] - attributes$AGE[j])
#'   }
#' XvC <- array(1, dim=c(g, g, 1))
#' # Now we are ready to fit the model
#' mod <- fit_p2(Y, Xn, Xn, XvD, XvC)
#' print(mod)
fit_p2 <- function(y, XnS, XnR, XvD, XvC, M = 0, seed = NULL, trace = FALSE,
                   init = NULL, penalized = FALSE, penSigma = NULL,
                   opt = nlminb, singular.ok = TRUE, ...)
{
  # Do some argument checking
  if(is.null(y) | is.null(XvD) | is.null(XvC)) stop("y, XvD and XvC must be provided \n")
  if(!is.matrix(y))
    {
     warning("network data should be provided as a matrix\n")
     y <- as.matrix(y)
    }
  if(ncol(y) != nrow(y)) stop("y must be a squared matrix\n")
  g <- ncol(y)
  if(!is.null(XnS) && !is.matrix(XnS)) stop("XnS must be a matrix\n")
  if(!is.null(XnR) && !is.matrix(XnR)) stop("XnR must be a matrix\n")
  if(!is.array(XvD) | length(dim(XvD))!=3) stop("XvD must be a 3-dim array\n")
  if(!is.array(XvC) | length(dim(XvC))!=3) stop("XvC must be a 3-dim array\n")
  if(!is.null(XnS) &&  nrow(XnS)!=g) stop("Wrong dimension of XnS\n")
  if(!is.null(XnR) &&  nrow(XnR)!=g) stop("Wrong dimension of XnR\n")
  if(dim(XvD)[2]!=g | dim(XvD)[1]!=g) stop("Wrong dimension of XvD\n")
  if(dim(XvC)[2]!=g | dim(XvC)[1]!=g) stop("Wrong dimension of XvC\n")
  if(!is.numeric(M) || M<0) M <- 0
  if(!is.matrix(penSigma)) penSigma <- NULL
  if(M>0) warning("IS as implemented by TMB is still experimental --> use the fitIS function instead\n")
  #  Create starting values
  kd <- dim(XvD)[3]
  kc <- dim(XvC)[3]
  map <- list()
  XnS.int <- XnS
  XnR.int <- XnR
  if(is.null(XnS) & !is.null(XnR))
    {
     XnS.int <- matrix(0, nrow = g, ncol = 1)
     map <- list(gamma1 = factor(NA))
   }
  if(!is.null(XnS) & is.null(XnR))
    {
     XnR.int <- matrix(0, nrow = g, ncol = 1)
     map <- list(gamma2 = factor(NA))
    }
  if(is.null(XnS) & is.null(XnR))
    {
     XnS.int <- matrix(0, nrow = g, ncol = 1)
     XnR.int <- matrix(0, nrow=g, ncol = 1)
     map <- list(gamma1 = factor(NA), gamma2 = factor(NA))
   }
  kr <- ncol(XnR.int)
  ks <- ncol(XnS.int)
  model.param <-  list(gamma1 = rep(0, ks) , gamma2 = rep(0, kr), delta1 = rep(0, kd),
                  delta2 = rep(0, kc), alpha=c(1, 0, 1), a=rep(0, g), b = rep(0, g))
  #  Create list of model data for optimization
  if(!penalized) flagpen <- 0
  else
  {
    if(is.null(penSigma)) flagpen <- 1 else flagpen <- 2
  }
  flagSigma <- if(is.null(penSigma)) 0 else  c(penSigma[1,1], penSigma[1,2], penSigma[2,2])
  model.data <- list(XnS = XnS.int, XnR = XnR.int, XvD = XvD, XvC = XvC, y = y,
                penflag = as.integer(flagpen), penSigma = as.vector(flagSigma))
  # Create AD function with data and parameters
  myseed <- if(is.null(seed)) sample(1:10^5, 1) else seed
  obj <- TMB::MakeADFun(data = model.data, parameters = model.param, random = c("a", "b"),
         DLL = "p2model", map = map, silent = !trace,
         MCcontrol = list(doMC = M>0, seed = myseed, n = M))
  if(trace) cat("\nfitting the model\n")
  start <- if(is.null(init)) obj$par else init
  lower.vect <-  if(singular.ok) c(rep(-Inf, length(obj$par)-3), 0, -Inf, 0) else -Inf
  mod <- try(opt(start, obj$fn, obj$gr, lower = lower.vect, ...), silent = TRUE)
  if(!is.character(mod))
  {
    par <- mod$par
    if(trace) cat("\ncomputing Hessian\n")
    sde <- try(TMB::sdreport(obj, par.fixed=par), silent=TRUE)
    if(is.character(sde))
      {
      	warning("\nNumerical issues in the computation of the Hessian\n")
      	hess <- try(numDeriv::jacobian(obj$gr, par, method="simple"), silent=TRUE)
      	if(is.character(hess))
      	 stop("\nNumerical issues in the computation of standard errors are too serious...
      	        Change optimizer, try out a better starting point or switch to penalized estimation\n")
        else  sde <- TMB::sdreport(obj, par.fixed=par, hessian.fixed=hess)
      }
    theta.vcov <-  sde$cov.fixed
    theta.se <- sqrt(diag(theta.vcov))
    Sigma <- matrix(c(sde$value[1], sde$value[2], sde$value[2], sde$value[3]), 2, 2)
    Sigma.se <- matrix(c(sde$sd[1], sde$sd[2], sde$sd[2], sde$sd[3]), 2, 2)
    rho <- sde$value[4]
    rho.se <- sde$sd[4]
    random <- sde$par.random
    random.se <- sqrt(sde$diag.cov.random)
    lnl <- obj$fn(par) + obj$env$report()$pen     #### eliminate the penalty
    lnl <- lnl + g * log(pi * 2)                  #### this constant is introduced by TMB
    res <- list(theta = par, loglik = -lnl,
              AIC = 2 * lnl + 2 * length(par), BIC = 2 * lnl + log(g) * length(par),
              XnS.null = is.null(XnS), XnR.null = is.null(XnR), seed = myseed,
              theta.vcov = theta.vcov, theta.se = theta.se,
              opt = opt, opt.details = mod, ADobj = obj, model.data = model.data,
              sdrep = sde, ranef = random, ranef.se = random.se,
              Sigma = Sigma, Sigma.se = Sigma.se, rho = rho, rho.se = rho.se, M = M,
              penflag = model.data$penflag, penSigma = model.data$penSigma)
   }
   else stop("Optimization did not converge. Change optimizer, try out a better starting point or
              switch to penalized estimation")
  # Assign S3 class values and return
  class(res)=c("p2")
  return(res)
}



#' @export
print.p2 <- function(x, ...) {
  cat("p2 model fit by ML \n")

  cat("Log likelihood at convergence:", x$loglik, "\n")

  cat("Model coefficients: \n")
  p <- length(x$theta)-3
  print.default(x$theta[1:p], ...)

  cat("Variance matrix of random effects: \n")

  print.default(x$Sigma)
  invisible(x)
}



#' @export
summary.p2 <- function(object,... ) {
  est <- object$theta
  se <- object$theta.se
  resultD <- resultC <- resultS <- resultR <- NULL
  condD <-  attr(est,"names")=="delta1"
  if(sum(condD)>0)  resultD <- cbind("Estimate" = est[condD], "Std. error" = se[condD])
  condC <-  attr(est,"names")=="delta2"
  if(sum(condC)>0) resultC <- cbind("Estimate" = est[condC], "Std. error" = se[condC])
  condS <-  attr(est,"names")=="gamma1"
  if(sum(condS)>0) resultS <- cbind("Estimate" = est[condS], "Std. error" = se[condS])
  condR <-  attr(est,"names")=="gamma2"
  if(sum(condR)>0) resultR <- cbind("Estimate" = est[condR], "Std. error" = se[condR])
  summary <- list(resultC=resultC,
                  resultD=resultD,
                  resultS=resultS,
                  resultR=resultR,
                  loglik=object$loglik,
                  AIC=object$AIC,
                  BIC=object$BIC,
                  vcov=object$theta.vcov,
                  Sigma=object$Sigma)
  class(summary) <- "summary.p2"
  return(summary)
}


#' @export
print.summary.p2 <- function(x, ... ) {
  cat("--------------------------------------------\n")
  cat("Approximate maximum likelihood estimation of p2 model\n") 
  cat("Log-likelihood at maximum", x$loglik,"\n")
  cat("Model selection criteria\n")
  cat("AIC = ", x$AIC, "  BIC =", x$BIC,"\n")
  cat("--------------------------------------------\n")
  cat("Density coefficients\n")
  print(x$resultD)
  cat("--------------------------------------------\n")
  cat("Reciprocity coefficients\n")
  print(x$resultC)
  cat("--------------------------------------------\n")
  if(!is.null(x$resultS)){
    cat("Sender coefficients\n")
    print(x$resultS)
    cat("--------------------------------------------\n")
  }
  if(!is.null(x$resultR)){
    cat("Receiver coefficients\n")
    print(x$resultR)
    cat("--------------------------------------------\n")
  }
  cat("Variance matrix of random effects\n")
  print(x$Sigma)
  cat("--------------------------------------------\n")
}



#' Fitting function for p1 models
#'
#' Estimates a p1 model via (exact) maximum likelihood.
#'
#' The function is included for compatibility with \code{fit_p2}, and
#' it might also be useful to obtain starting values.
#' @param y An adjacency matrix of size \eqn{g \times g}{g x g}.
#' @param XnS Matrix of sender effects of size \eqn{g \times k_s}{g x ks}.
#' @param XnR Matrix of receiver effects of size \eqn{g \times k_r}{g x kr}.
#' @param XvD Density effects, a 3-dim array of size \eqn{g \times g \times
#'   k_d}{g x g x kd}.
#' @param XvC Reciprocity effects, a 3-dim array of size \eqn{g \times g \times
#'   k_c}{g x g x kc}.
#' @param trace TRUE for tracing information during the estimation. Default is FALSE.
#' @param init Optional starting value for model parameters. Default is NULL.
#' @param opt Name of the optimising function. Default is nlminb.
#' @param ... Optional arguments passed to the optimiser.
#' @return The returned value is an object of class
#' \code{"p1"}, a list  containing the following components:
#' \item{\code{theta}}{the vector of model estimates.}
#' \item{\code{loglik}}{the log likelihood function at the estimate.}
#' \item{\code{AIC, BIC}}{model selection criteria at the estimate.}
#' \item{\code{XnS.null, XnR.null}}{logical, flagging whether there
#'   are no sender or receiver effects, respectively.}
#' \item{\code{seed}}{the random seed used for estimation (when \code{M}>0).}
#' \item{\code{theta.cov}}{Variance matrix of estimates.}
#' \item{\code{theta.se}}{Standard errors of estimates.}
#' \item{\code{opt}}{Optimiser employed for estimation.}
#' \item{\code{opt.details}}{Object returned by the optimiser.}
#' \item{\code{ADobj}}{Object returned by \code{TMB:MakeADFun}, containing
#' the log likelihood function and its gradient (among its components).}
#' \item{\code{model.data}}{List containing  all the model data, fed
#' to \code{TMB:MakeADFun}.}
#' \item{\code{sdrep}}{Summary object returned by \code{TMB::sdreport}. This is an
#' object with many slots, useful for TMB users.}
#' @import NetData
#' @export
#' @author Ruggero Bellio
#' @examples
#' # Analysis of the  kracknets data from the NetData package
#' library(NetData)
#' data(kracknets)
#' # data preparation
#' g <- 21
#' Y <- matrix(0, g, g)
#' ind <-1
#' for(i in 1:nrow(friendship_data_frame)){
#'    sele <- friendship_data_frame[i, ]
#'    Y[sele$ego, sele$alter] <- sele$friendship_tie
#'    }
#' Xn <- model.matrix(~ AGE + TENURE, attributes)[, -1]
#' XvD <- array(1, dim=c(g, g, 4))
#' for(i in 1:g)
#'  for(j in 1:g){
#'     XvD[i, j, 2] <- as.numeric(attributes$DEPT[i]==attributes$DEPT[j])
#'     XvD[i, j, 3] <- as.numeric(attributes$LEVEL[i]==attributes$LEVEL[j])
#'     XvD[i, j, 4] <- abs(attributes$AGE[i] - attributes$AGE[j])
#'   }
#' XvC <- array(1, dim=c(g, g, 1))
#' # Now we are ready to fit the model
#' mod <- fit_p1(Y, Xn, Xn, XvD, XvC)
#' summary(mod)
fit_p1 <- function(y, XnS, XnR, XvD, XvC, trace=FALSE, init=NULL,  opt=nlminb,...)
{
  # Do some argument checking
  if(is.null(y) | is.null(XvD) | is.null(XvC)) stop("y, XvD and XvC must be provided \n")
  if(!is.matrix(y))
  {
    warning("network data should be provided as a matrix\n")
    y <- as.matrix(y)
  }
  if(ncol(y) != nrow(y)) stop("y must be a squared matrix\n")
  g <- ncol(y)
  if(!is.null(XnS) && !is.matrix(XnS)) stop("XnS must be a matrix\n")
  if(!is.null(XnR) && !is.matrix(XnR)) stop("XnR must be a matrix\n")
  if(!is.array(XvD) | length(dim(XvD))!=3) stop("XvD must be a 3-dim array\n")
  if(!is.array(XvC) | length(dim(XvC))!=3) stop("XvC must be a 3-dim array\n")
  if(!is.null(XnS) &&  nrow(XnS)!=g) stop("Wrong dimension of XnS\n")
  if(!is.null(XnR) &&  nrow(XnR)!=g) stop("Wrong dimension of XnR\n")
  if(dim(XvD)[2]!=g | dim(XvD)[1]!=g) stop("Wrong dimension of XvD\n")
  if(dim(XvC)[2]!=g | dim(XvC)[1]!=g) stop("Wrong dimension of XvC\n")
  #  Create starting values
  kd <- dim(XvD)[3]
  kc <- dim(XvC)[3]
  map <- list(alpha=factor(rep(NA, 3)), a=factor(rep(NA, g)), b=factor(rep(NA, g)))
  XnS.int <- XnS
  XnR.int <- XnR
  if(is.null(XnS) & !is.null(XnR))
  {
    XnS.int <- matrix(0, nrow = g, ncol = 1)
    map <- list(gamma1 = factor(NA), alpha=factor(rep(NA, 3)),
                a=factor(rep(NA, g)), b=factor(rep(NA, g)))
  }
  if(!is.null(XnS) & is.null(XnR))
  {
    XnR.int <- matrix(0, nrow = g, ncol = 1)
    map <- list(gamma2 = factor(NA), alpha = factor(rep(NA, 3)),
                a = factor(rep(NA, g)), b = factor(rep(NA, g)))
  }
  if(is.null(XnS) & is.null(XnR))
  {
    XnS.int <- matrix(0, nrow = g, ncol = 1)
    XnR.int <- matrix(0, nrow=g, ncol = 1)
    map <- list(gamma1 = factor(NA), gamma2 = factor(NA), alpha=factor(rep(NA, 3)),
                a=factor(rep(NA, g)), b=factor(rep(NA, g)))
  }
  kr <- ncol(XnR.int)
  ks <- ncol(XnS.int)
  model.param <-  list(gamma1 = rep(0, ks) , gamma2 = rep(0, kr), delta1 = rep(0, kd),
                       delta2 = rep(0, kc), alpha=c(1, 0, 1), a=rep(0, g), b = rep(0, g))
  #  Create list of model data for optimization
  model.data <- list(XnS = XnS.int, XnR = XnR.int, XvD = XvD, XvC = XvC, y = y,
                     penflag = as.integer(0), penSigma = as.vector(0))
  # Create AD function with data and parameters
  obj <- TMB::MakeADFun(data = model.data, parameters = model.param, DLL="p2model",
                        map=map, silent=!trace)
  if(trace) cat("\nrunning TMB program\n")
  start <- if(is.null(init)) obj$par else init
  mod <- opt(start, obj$fn, obj$gr, hessian=obj$he,...)
  par <- mod$par
  # Compute report
  sde <- TMB::sdreport(obj, par.fixed=mod$par)
  theta.vcov <-  sde$cov.fixed
  theta.se <- sqrt(diag(theta.vcov))
  lnl <- obj$fn(par)
  res <- list(theta = par, loglik = -lnl, AIC = 2 * lnl + 2 * length(par),
              BIC = 2*lnl + log(g) * length(par), XnS.null = is.null(XnS),
              XnR.null = is.null(XnR), theta.vcov = theta.vcov, theta.se = theta.se,
              opt.details = mod, ADobj = obj, model.data = model.data, sdrep = sde)
  # Assign S3 class values and return
  class(res) = c("p1")
  return(res)
}


#' @export
print.p1 <- function(x, ...) {
  cat("p1 model fit by ML \n")

  cat("Log likelihood at convergence:", x$loglik, "\n")

  cat("Model coefficients: \n")
  p <- length(x$theta)
  print.default(x$theta[1:p], ...)
  invisible(x)
}


#' @export
print.summary.p1 <- function(x, ... ) {
   cat("--------------------------------------------\n")
   cat("Maximum likelihood estimation of p1 model\n") 
   cat("Log-likelihood at maximum", x$loglik,"\n")
   cat("Model selection criteria\n")
   cat("AIC = ", x$AIC, "  BIC =", x$BIC,"\n")
   cat("--------------------------------------------\n")
   cat("Density coefficients\n")
   print(x$resultD)
   cat("--------------------------------------------\n")
   cat("Reciprocity coefficients\n")
   print(x$resultC)
   cat("--------------------------------------------\n")
   if(!is.null(x$resultS)){
     cat("Sender coefficients\n")
     print(x$resultS)
     cat("--------------------------------------------\n")
   }
   if(!is.null(x$resultR)){
     cat("Receiver coefficients\n")
     print(x$resultR)
     cat("--------------------------------------------\n")
   }
}


#' @export
summary.p1 <- function(object,... ) {
  est <- object$theta
  se <- object$theta.se
  resultD <- resultC <- resultS <- resultR <- NULL
  condD <-  attr(est,"names")=="delta1"
  if(sum(condD)>0)  resultD <- cbind("Estimate" = est[condD], "Std. error" = se[condD])
  condC <-  attr(est,"names")=="delta2"
  if(sum(condC)>0) resultC <- cbind("Estimate" = est[condC], "Std. error" = se[condC])
  condS <-  attr(est,"names")=="gamma1"
  if(sum(condS)>0) resultS <- cbind("Estimate" = est[condS], "Std. error" = se[condS])
  condR <-  attr(est,"names")=="gamma2"
  if(sum(condR)>0) resultR <- cbind("Estimate" = est[condR], "Std. error" = se[condR])
  summary <- list(resultC=resultC,
                  resultD=resultD,
                  resultS=resultS,
                  resultR=resultR,
                  loglik=object$loglik,
                  AIC=object$AIC,
                  BIC=object$BIC,
                  vcov=object$theta.vcov)
  class(summary) <- "summary.p1"
  return(summary)
}



#' Plot receiver and sender estimated random effects
#'
#' @param fit An object returned by \code{fit_p2}.
#' @param resolution Graphical resolution.
#' @param seed Random seed.
#' @param col Type of colors for the output.
#' @param mai Margins for plots.
#' @param arrow.size Arrow size.
#' @param ... Optional graphical arguments.
#' @export
#' @importFrom graphics par plot
plot_effects_p2 <- function(fit, resolution = round(length(fit$ranef)) / 2, seed = 2,
                    col = rev(grDevices::heat.colors(10)),
                    mai = c(1, 1, 1, 1) * 0.5, arrow.size = 0.5,...)
{
  y <- fit$model.data$y
  ranef <- fit$ranef
  old <- options(stringAsFactors = TRUE)
  z <- ranef
  g.cd <- igraph::graph.adjacency(y)
  mypalette  <- grDevices::colorRampPalette(col)
  breaks <- seq(min(z), max(z), l = resolution + 1)
  facz <- as.numeric(cut(z, breaks = breaks, include.lowest=TRUE))
  cond <-  is.element(1:resolution, facz)
  ind <- (1:resolution)[!cond] + 1
  breaks <- breaks[-ind]
  colors <-  mypalette(resolution)[facz]
  cola <- colors[1:ncol(y)]
  par(mfrow = c(1, 2), mai = mai)
  par(...)
  set.seed(seed)
  plot(g.cd, vertex.color = cola, main = "Sender effects", edge.arrow.size = arrow.size)
  fields::image.plot(legend.only = TRUE, zlim = range(z), breaks = breaks,
             col = unique(colors[order(z)]))
  set.seed(seed)
  colb <- colors[1:ncol(y) + ncol(y)]
  par(...)
  plot(g.cd, vertex.color = colb, main = "Receiver effects", edge.arrow.size = arrow.size)
  on.exit(options(old), add = TRUE)
  invisible(NULL)
}


#' @keywords internal
.recTheta <- function(theta, condS, condR, ks, kr, kd, kc)
{
   if(condS & condR)
     {
     	gamma1 <- theta[1:ks]
     	gamma2 <- theta[(ks+1):(ks+kr)]
     	delta1 <- theta[(kr+ks+1):(kr+ks+kd)]
     	delta2 <- theta[(kr+ks+kd+1):(kr+ks+kd+kc)]
     	}
   if(condS & !condR)
     {
     	gamma1 <- theta[1:ks]
     	gamma2 <- 0
     	delta1 <- theta[(ks+1):(ks+kd)]
     	delta2 <- theta[(ks+kd+1):(ks+kd+kc)]
     	}
   if(!condS & condR)
     {
     	gamma1 <- 0
     	gamma2 <- theta[1:kr]
     	delta1 <- theta[(kr+1):(kr+kd)]
     	delta2 <- theta[(kr+kd+1):(kr+kd+kc)]
     	}
   if(!condS & !condR)
     {
     	gamma1 <- 0
     	gamma2 <- 0
     	delta1 <- theta[1:kd]
     	delta2 <- theta[(kd+1):(kd+kc)]
     	}
   return(list(gamma1=gamma1, gamma2=gamma2, delta1=delta1, delta2=delta2))
}


#' Simulate networks from a p2 model
#'
#' The function can be used to simulate from an estimated model or from a given
#' parameter value. In the latter case, the model settings are included in a fitted model
#' object.
#'
#' @param objfit An object produced by a call to \code{fit_p2}.
#' @param nsim Number of simulated networks.
#' @param conditional Should the random effects be set at the estimated values,
#' or else simulated from the normal distribution? Default is FALSE.
#' @param theta Parameter value to be used for the simulation. In case it is \code{NULL},
#' the \code{theta} slot of \code{objfit} is used.
#' @param Sigma Variance matrix of random effects to be used for the simulation. In case
#' it is \code{NULL}, the \code{Sigma} slot of \code{objfit} is used.
#' @return Either a single network (when \code{nsim}=1) or a list of \code{nsim} networks.
#' @importFrom stats rbinom simulate
#' @importFrom network network
#' @export
#'
simulate_p2 <- function(objfit, nsim = 100, conditional = FALSE, theta = NULL, Sigma = NULL)
{
   random <- !is.null(objfit$ranef)
   # allocate all the networks using the ergm package
   y <-  network(objfit$model.data$y, directed = TRUE)
   nw.sim <- simulate(y ~ edges, coef = 0, nsim = nsim)
   g <- ncol(y[,])
   if(random & conditional)
    {
     a0 <- objfit$ranef[1:g]
     b0 <- objfit$ranef[g+1:g]
    }
   else a0 <- b0 <- rep(0,g)
   # recover matrices
   XnS <- objfit$model.data$XnS; ks <- ncol(XnS)
   XnR <- objfit$model.data$XnR; kr <- ncol(XnR)
   XvC <- objfit$model.data$XvC; kc <- dim(XvC)[3]
   XvD <- objfit$model.data$XvD; kd <- dim(XvD)[3]
   # recover parameters
   condS <- !(objfit$XnS.null)
   condR <- !(objfit$XnR.null)
   thetasim <- if(is.null(theta)) objfit$theta else theta
   Sigmasim <- if(is.null(Sigma) & random) objfit$Sigma else Sigma
   recT <- .recTheta(thetasim, condS, condR, ks, kr, kd, kc)
   gamma1 <- recT$gamma1
   gamma2 <- recT$gamma2
   delta1 <- recT$delta1
   delta2 <- recT$delta2
   # simulate all the random effects
   xs1 <- as.vector(XnS %*% gamma1)
   xr2 <- as.vector(XnR %*% gamma2)
   # simulate all the networks
   for(s in 1:nsim)
   {
   	Y <- matrix(0, g, g)
    a <- a0
    b <- b0
    if(random & !conditional)
   	  {
   	   v <- MASS::mvrnorm(n = g, mu = rep(0,2), Sigma = Sigmasim, tol = 1e-6, empirical = FALSE)
   	   a <- v[,1]
   	   b <- v[,2]
   	  }
    for(i in 1:(g-1))
     for(j in (i+1):g)
   	  {
   	    alpha.i <-  xs1[i] + a[i]
   	    alpha.j <-  xs1[j] + a[j]
   	    beta.i <-   xr2[i] + b[i]
   	    beta.j <-   xr2[j] + b[j]
   	    mu.ij <- sum(XvD[i,j,] * delta1)
   	    mu.ji <- sum(XvD[j,i,] * delta1)
   	    rho.ij <- sum(XvC[i,j,] * delta2)
   	    xi1 <- mu.ij + alpha.i + beta.j
   	    xi2 <-  mu.ji + alpha.j + beta.i
        xi3 <-  mu.ij + mu.ji + alpha.i + beta.j + alpha.j + beta.i + rho.ij
   	  	den <- 1.0 + exp(xi1) + exp(xi2) + exp(xi3)
   	  	p1 <- (exp(xi1) + exp(xi3)) / den
   	  	y1 <- rbinom(1, prob=p1, size=1)
   	  	Y[i,j] <- y1
   	  	p2.1 <- exp(y1 * xi1 +xi2 + y1 * rho.ij) / (exp(xi1*y1) + exp(y1 * xi1 +xi2 + y1 * rho.ij))
   	  	Y[j,i] <- rbinom(1, prob=p2.1, size=1)
   	   }
   if(nsim>1) nw.sim[[s]] <- network(Y, directed=TRUE)
   }
   if(nsim == 1)
   {
   	 out <- network(Y, directed=TRUE)
    }
    else
   {
    out <- list(networks = nw.sim)
    }
  return(out)
}


#' Computes GOF statistics for a p2 model fit
#'
#' This is a slight modification of the function \code{gof.ergmm} from the
#' package \code{latentnet}.
#'
#' The function makes only a minor modification to the function \code{gof.ergmm} from the
#' package \code{latentnet}. The only change is actually that the function \code{simulate.p2}
#' is employed to generate the networks.
#' This means that this function is based on 'statnet' project software
#' (\url{http://statnet.org}). For license and citation information see \url{http://statnet.org/attribution}.
#' It would have been much preferable to pass the simulated networks to \code{gof.ergmm}, but
#' unfortunately such functionality is not supported, therefore it has been necessary to duplicate
#' the function.
#'
#' @param objfit An object produced by a call to \code{fit_p2}.
#' @param nsim Number of simulated networks.
#' @param GOF A formula object to specify the statistics to use to diagnosis the goodness-of-fit
#' of the model, among those suitable for directed networks. See the help file of \code{gof.ergmm}
#' for more details.
#' @param verbose	Provides verbose information on the progress of the simulation. Default is FALSE.
#' @param conditional Should the random effects be set at the estimated values,
#' or else simulated from the normal distribution? Default is FALSE.
#' @return Like \code{gof.ergmm}, it returns an object of class \code{gofobject}.
#' @importFrom stats as.formula quantile
#' @importFrom network network network.size
#' @import ergm statnet.common
#' @export
gof_p2 <-
function(objfit, nsim=100, GOF = ~idegree + odegree + distance, verbose = FALSE,
         conditional = FALSE)
 {
   SimNetworkSeriesObj <- simulate_p2(objfit, nsim = nsim, conditional = conditional)
   nw <- network(objfit$model.data$y, directed=TRUE)
   all.gof.vars <- all.vars(GOF)
   for (i in seq(along = all.gof.vars))
       {
        all.gof.vars[i] <- match.arg(all.gof.vars[i], c("distance",
            "espartners", "dspartners", "odegree", "idegree",  "triadcensus"))
       }
   GOF <- as.formula(paste("~", paste(all.gof.vars, collapse = "+")))
   pval.triadcensus <- pval.dist <- pval.deg <- pval.espart <- pval.espart <- NULL
   obs.triadcensus <- pobs.triadcensus <- sim.triadcensus <- psim.triadcensus <- pval.triadcensus <- bds.triadcensus <- NULL
   obs.dist <- pobs.dist <- sim.dist <- psim.dist <- pval.dist <- bds.dist <- NULL
   obs.espart <- pobs.espart <- sim.espart <- psim.espart <- pval.espart <- bds.espart <- NULL
   obs.dspart <- pobs.dspart <- sim.dspart <- psim.dspart <- pval.dspart <- bds.dspart <- NULL
   obs.ideg <- pobs.ideg <- sim.ideg <- psim.ideg <- pval.ideg <- bds.ideg <- pval.ideg <- NULL
   obs.odeg <- pobs.odeg <- sim.odeg <- psim.odeg <- pval.odeg <- bds.odeg <- pval.odeg <- NULL
   n <- network.size(nw)
   if("distance" %in% all.gof.vars) {
        obs.dist <- ergm::ergm.geodistdist(nw)
        obs.dist[obs.dist == Inf] <- n
        sim.dist <- array(0, dim = c(nsim, n))
        dimnames(sim.dist) <- list(paste(c(1:nsim)), paste(1:n))
    }
   if("odegree" %in% all.gof.vars) {
        mesp <- paste("c(", paste(0:(n - 1), collapse = ","),
            ")", sep = "")
        obs.odeg <- summary(as.formula(paste("nw ~ odegree(",
            mesp, ")", sep = "")), drop = FALSE)
        sim.odeg <- array(0, dim = c(nsim, n))
        dimnames(sim.odeg) <- list(paste(c(1:nsim)), paste(0:(n -
            1)))
        names(obs.odeg) <- dimnames(sim.odeg)[[2]]
    }
   if("idegree" %in% all.gof.vars) {
        mesp <- paste("c(", paste(0:(n - 1), collapse = ","),
            ")", sep = "")
        obs.ideg <- summary(as.formula(paste("nw ~ idegree(",
            mesp, ")", sep = "")), drop = FALSE)
        sim.ideg <- array(0, dim = c(nsim, n))
        dimnames(sim.ideg) <- list(paste(c(1:nsim)), paste(0:(n -
            1)))
        names(obs.ideg) <- dimnames(sim.ideg)[[2]]
    }
   if("espartners" %in% all.gof.vars) {
        mesp <- paste("c(", paste(0:(network.size(nw) - 2), collapse = ","),
            ")", sep = "")
        obs.espart <- summary(as.formula(paste("nw ~ esp(", mesp,
            ")", sep = "")), drop = FALSE)
        sim.espart <- array(0, dim = c(nsim, n - 1))
        dimnames(sim.espart) <- list(paste(c(1:nsim)), paste(0:(n -
            2)))
    }
   if("dspartners" %in% all.gof.vars) {
        mesp <- paste("c(", paste(0:(network.size(nw) - 2), collapse = ","),
            ")", sep = "")
        obs.dspart <- summary(as.formula(paste("nw ~ dsp(", mesp,
            ")", sep = "")), drop = FALSE)
        sim.dspart <- array(0, dim = c(nsim, n - 1))
        dimnames(sim.dspart) <- list(paste(c(1:nsim)), paste(0:(n -
            2)))
    }
   if("triadcensus" %in% all.gof.vars) {
        triadcensus <- 0:15
        namestriadcensus <- c("003", "012", "102", "021D",
                "021U", "021C", "111D", "111U", "030T", "030C",
                "201", "120D", "120U", "120C", "210", "300")
        triadcensus.formula <- "~ triadcensus(0:15)"
        obs.triadcensus <- summary(as.formula(paste("nw", triadcensus.formula,
            sep = "")), drop = FALSE)
        sim.triadcensus <- array(0, dim = c(nsim, length(triadcensus)))
        dimnames(sim.triadcensus) <- list(paste(c(1:nsim)), namestriadcensus)
        names(obs.triadcensus) <- namestriadcensus
    }
  if(verbose) {
        cat("\nCollating simulations\n")
    }
    for(i in 1:nsim) {
        if(verbose) {
            cat("\nCalculating statistics for simulation", i,
                "\n")
        }
        if("distance" %in% all.gof.vars) {
            sim.dist[i, ] <- ergm::ergm.geodistdist(SimNetworkSeriesObj[["networks"]][[i]])
        }
        if("idegree" %in% all.gof.vars) {
            mesp <- paste("c(", paste(0:(n - 1), collapse = ","),
                ")", sep = "")
            gi <- SimNetworkSeriesObj[["networks"]][[i]]
            sim.ideg[i, ] <- summary(as.formula(paste("gi ~ idegree(",
                mesp, ")", sep = "")), drop = FALSE)
        }
        if("odegree" %in% all.gof.vars) {
            mesp <- paste("c(", paste(0:(n - 1), collapse = ","),
                ")", sep = "")
            gi <- SimNetworkSeriesObj[["networks"]][[i]]
            sim.odeg[i, ] <- summary(as.formula(paste("gi ~ odegree(",
                mesp, ")", sep = "")), drop = FALSE)
        }
       if("espartners" %in% all.gof.vars) {
            gi <- SimNetworkSeriesObj[["networks"]][[i]]
            mesp <- paste("c(", paste(0:(network::network.size(gi) - 2),
                collapse = ","), ")", sep = "")
            sim.espart[i, ] <- summary(as.formula(paste("gi ~ esp(",
                mesp, ")", sep = "")), drop = FALSE)
        }
       if("dspartners" %in% all.gof.vars) {
            gi <- SimNetworkSeriesObj[["networks"]][[i]]
            mesp <- paste("c(", paste(0:(network.size(gi) - 2),
                collapse = ","), ")", sep = "")
            sim.dspart[i, ] <- summary(as.formula(paste("gi ~ dsp(",
                mesp, ")", sep = "")), drop = FALSE)
        }
        if("triadcensus" %in% all.gof.vars) {
            gi <- SimNetworkSeriesObj[["networks"]][[i]]
            sim.triadcensus[i, ] <- summary(as.formula(paste("gi",
                triadcensus.formula, sep = "")), drop = FALSE)
        }
    }

    if("distance" %in% all.gof.vars) {
        pval.dist <- apply(sim.dist <= obs.dist[col(sim.dist)],
            2, mean)
        pval.dist.top <- apply(sim.dist >= obs.dist[col(sim.dist)],
            2, mean)
        pval.dist <- cbind(obs.dist, apply(sim.dist, 2, min),
            apply(sim.dist, 2, mean), apply(sim.dist, 2, max),
            pmin(1, 2 * pmin(pval.dist, pval.dist.top)))
        dimnames(pval.dist)[[2]] <- c("obs", "min", "mean", "max",
            "MC p-value")
        pobs.dist <- obs.dist/sum(obs.dist)
        psim.dist <- sweep(sim.dist, 1, apply(sim.dist, 1, sum),
            "/")
        psim.dist[is.na(psim.dist)] <- 1
        bds.dist <- apply(psim.dist, 2, quantile, probs = c(0.025,
            0.975))
    }
    if("idegree" %in% all.gof.vars) {
        pval.ideg <- apply(sim.ideg <= obs.ideg[col(sim.ideg)],
            2, mean)
        pval.ideg.top <- apply(sim.ideg >= obs.ideg[col(sim.ideg)],
            2, mean)
        pval.ideg <- cbind(obs.ideg, apply(sim.ideg, 2, min),
            apply(sim.ideg, 2, mean), apply(sim.ideg, 2, max),
            pmin(1, 2 * pmin(pval.ideg, pval.ideg.top)))
        dimnames(pval.ideg)[[2]] <- c("obs", "min", "mean", "max",
            "MC p-value")
        pobs.ideg <- obs.ideg/sum(obs.ideg)
        psim.ideg <- sweep(sim.ideg, 1, apply(sim.ideg, 1, sum),
            "/")
        psim.ideg[is.na(psim.ideg)] <- 1
        bds.ideg <- apply(psim.ideg, 2, quantile, probs = c(0.025,
            0.975))
    }
    if("odegree" %in% all.gof.vars) {
        pval.odeg <- apply(sim.odeg <= obs.odeg[col(sim.odeg)],
            2, mean)
        pval.odeg.top <- apply(sim.odeg >= obs.odeg[col(sim.odeg)],
            2, mean)
        pval.odeg <- cbind(obs.odeg, apply(sim.odeg, 2, min),
            apply(sim.odeg, 2, mean), apply(sim.odeg, 2, max),
            pmin(1, 2 * pmin(pval.odeg, pval.odeg.top)))
        dimnames(pval.odeg)[[2]] <- c("obs", "min", "mean", "max",
            "MC p-value")
        pobs.odeg <- obs.odeg/sum(obs.odeg)
        psim.odeg <- sweep(sim.odeg, 1, apply(sim.odeg, 1, sum),
            "/")
        psim.odeg[is.na(psim.odeg)] <- 1
        bds.odeg <- apply(psim.odeg, 2, quantile, probs = c(0.025,
            0.975))
    }
        if("espartners" %in% all.gof.vars) {
        pval.espart <- apply(sim.espart <= obs.espart[col(sim.espart)],
            2, mean)
        pval.espart.top <- apply(sim.espart >= obs.espart[col(sim.espart)],
            2, mean)
        pval.espart <- cbind(obs.espart, apply(sim.espart, 2,
            min), apply(sim.espart, 2, mean), apply(sim.espart,
            2, max), pmin(1, 2 * pmin(pval.espart, pval.espart.top)))
        dimnames(pval.espart)[[2]] <- c("obs", "min", "mean",
            "max", "MC p-value")
        pobs.espart <- obs.espart/sum(obs.espart)
        psim.espart <- sweep(sim.espart, 1, apply(sim.espart,
            1, sum), "/")
        psim.espart[is.na(psim.espart)] <- 1
        bds.espart <- apply(psim.espart, 2, quantile, probs = c(0.025,
            0.975))
    }
    if("dspartners" %in% all.gof.vars) {
        pval.dspart <- apply(sim.dspart <= obs.dspart[col(sim.dspart)],
            2, mean)
        pval.dspart.top <- apply(sim.dspart >= obs.dspart[col(sim.dspart)],
            2, mean)
        pval.dspart <- cbind(obs.dspart, apply(sim.dspart, 2,
            min), apply(sim.dspart, 2, mean), apply(sim.dspart,
            2, max), pmin(1, 2 * pmin(pval.dspart, pval.dspart.top)))
        dimnames(pval.dspart)[[2]] <- c("obs", "min", "mean",
            "max", "MC p-value")
        pobs.dspart <- obs.dspart/sum(obs.dspart)
        psim.dspart <- sweep(sim.dspart, 1, apply(sim.dspart,
            1, sum), "/")
        psim.dspart[is.na(psim.dspart)] <- 1
        bds.dspart <- apply(psim.dspart, 2, quantile, probs = c(0.025,
            0.975))
    }
    if("triadcensus" %in% all.gof.vars) {
        pval.triadcensus <- apply(sim.triadcensus <= obs.triadcensus[col(sim.triadcensus)],
            2, mean)
        pval.triadcensus.top <- apply(sim.triadcensus >= obs.triadcensus[col(sim.triadcensus)],
            2, mean)
        pval.triadcensus <- cbind(obs.triadcensus, apply(sim.triadcensus,
            2, min), apply(sim.triadcensus, 2, mean), apply(sim.triadcensus,
            2, max), pmin(1, 2 * pmin(pval.triadcensus, pval.triadcensus.top)))
        dimnames(pval.triadcensus)[[2]] <- c("obs", "min", "mean",
            "max", "MC p-value")
        pobs.triadcensus <- obs.triadcensus/sum(obs.triadcensus)
        psim.triadcensus <- sweep(sim.triadcensus, 1, apply(sim.triadcensus,
            1, sum), "/")
        psim.triadcensus[is.na(psim.triadcensus)] <- 1
        bds.triadcensus <- apply(psim.triadcensus, 2, quantile,
            probs = c(0.025, 0.975))
    }
    returnlist <- list(n, pval.triadcensus, pval.dist,
        pval.ideg, pval.odeg, pval.deg, pval.espart, pval.dspart,
        obs.triadcensus, pobs.triadcensus, sim.triadcensus,
        psim.triadcensus, pval.triadcensus, bds.triadcensus,
        obs.dist, pobs.dist, sim.dist, psim.dist, pval.dist,
        bds.dist, obs.ideg, pobs.ideg, sim.ideg, psim.ideg, pval.ideg,
        bds.ideg, obs.odeg, pobs.odeg, sim.odeg, psim.odeg, pval.odeg,
        bds.odeg, obs.espart, pobs.espart, sim.espart, psim.espart,
        pval.espart, bds.espart, obs.dspart, pobs.dspart, sim.dspart,
        psim.dspart, pval.dspart, bds.dspart, GOF)
    names(returnlist) <- c("network.size", "summary.triadcensus",
        "summary.dist", "summary.ideg", "summary.odeg", "summary.deg",
        "summary.espart", "summary.dspart",
       "obs.triadcensus", "pobs.triadcensus", "sim.triadcensus",
        "psim.triadcensus", "pval.triadcensus", "bds.triadcensus",
        "obs.dist", "pobs.dist", "sim.dist", "psim.dist", "pval.dist",
        "bds.dist", "obs.ideg", "pobs.ideg", "sim.ideg", "psim.ideg",
        "pval.ideg", "bds.ideg", "obs.odeg", "pobs.odeg", "sim.odeg",
        "psim.odeg", "pval.odeg", "bds.odeg",  "obs.espart",
        "pobs.espart", "sim.espart", "psim.espart", "pval.espart",
        "bds.espart", "obs.dspart", "pobs.dspart", "sim.dspart",
        "psim.dspart", "pval.dspart", "bds.dspart", "GOF")
    class(returnlist) <- "gofobject"
    returnlist
}










