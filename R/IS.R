###logl is the loglik (not minus!)
#' @keywords internal
.crit <- function(mle, logl, lapfit, penalized = FALSE)
{
  k <- length(mle)
  if(penalized) {
    alpha <- mle[(k-2):k]
    Sigma <- a2S(alpha)
    DS <-  Sigma[1,1] * Sigma[2,2] - Sigma[2,1]^2
    logl <- logl  - 0.5 * log(DS)
   }
  g <- length(lapfit$ranef) / 2
  return(list(loglik = logl, AIC = -2 * logl + 2 * k, BIC = - 2 * logl + log(g) * k))
}


#' @importFrom stats nlminb
#' @importFrom TMB MakeADFun
#' @keywords internal
.likIStidy <- function(theta, lapfit, listM, ncores = 1)
{
  obj.par <- .para.fun(theta, lapfit)
  data <- lapfit$model.data
  objh <- MakeADFun(data = data, parameters = obj.par$para, DLL = "p2model",
                         map = obj.par$map, silent = TRUE)
  u.fit <- nlminb(objh$par, objh$fn, objh$gr, objh$he)
  u.theta<- u.fit$par ###here the penalty is already included!!!
  h0 <- (-1) * u.fit$obj
  G <-  objh$he(u.fit$par)
  Hinv <- solve(G)
  C.theta <- t(chol(Hinv))
  myf <- function(v) {
      x <- u.theta + C.theta %*% v
      hx <- (-1) * objh$fn(x)
      return(exp(hx - h0 + 0.5 * sum(v^2)))
     }
  ogg <-  unlist(parallel::mclapply(listM, myf, mc.cores = ncores))
  out <- log(mean(ogg)) + h0 + log(det(C.theta))
  return(-out)
}


#' Improve the fit of p2 models
#'
#' Estimates a p2 model using \emph{parameter dependent Laplace Importance Sampling}.
#'
#' The function can be called only after fitting a model by \code{fit_p2}, so that
#' all the model settings are inherited from those employed in the call to\code{fit_p2}.
#' The function can take advantage of multiple cores for parallel evaluation of the
#' samples simulated by the importance sampler. The optimiser used is \code{ucminf}, which
#' can exploit the estimated Hessian from  \code{fit_p2} in a Quasi-Newton optimisation.
#'
#' @param objfit An object produced by a call to \code{fit_p2}.
#' @param M number of importance sampling simulations. Default is \code{M=5000}.
#' @param ncores Number of cores to employ.
#' @param init  Starting point for the model parameter value. If it is
#' \code{NULL}, the \code{theta} slot of \code{objfit} is used.
#' @param init.hess  Should the hessian matrix from \code{fit_p2} be used
#' as the initial Hessian by the \code{ucminf} optimiser? Default is \code{TRUE}.
#' @param seed Random seed for the importance sampler.
#' @param trace,hessian,grtol,xtol Arguments for the \code{ucminf} optimiser. For
#' their meaning, see the help file for that function.
#' @return The returned value is an object of class
#' \code{"p2"}, a list  containing the same components listed for \code{fit_p2}.
#' @export
#' @importFrom stats rnorm
#' @useDynLib p2model
#' @author Ruggero Bellio
#' @seealso \code{\link{fit_p2}}, \code{\link{ucminf}}
#' @examples
#' \dontrun{
#' mod <- fit_p2(Y, Xn, Xn, XvD, XvC)
#' nc <- parallel::detectCores()
#' modIS <- fit_p2_IS(mod, ncores=nc)}
#'
fit_p2_IS <- function(objfit, M = 5000, ncores = 1, init = NULL, init.hess = TRUE, seed = NULL,
                  trace = 0, hessian = 1, grtol = 1e-6, xtol = 1e-12)
{
  start <- if(is.null(init)) objfit$theta else init
  k <- length(start)
  H0 <- diag(k)
  if(init.hess) {
   H0 <- qr.solve(objfit$sdrep$cov.fixed)
  }
  invhessian.lt  <- solve(H0)[lower.tri(H0, diag = TRUE)]
  myseed <- if(is.null(seed)) sample(1:10^5, 1) else seed
  g <-  ncol(objfit$model.data$y)
  set.seed(myseed)
  matV <- matrix(rnorm(M * 2 * g), ncol = 2 * g)
  listM <- split(matV, c(row(matV)))
  mod <- ucminf::ucminf(start, .likIStidy, listM = listM, ncores = ncores,  lapfit = objfit,
                        control = list(invhessian.lt = invhessian.lt, trace = trace,
                        grad = "central"), hessian = hessian)
  theta.vcov <- theta.se <- NULL
  if(hessian){
    theta.vcov <- solve(mod$hessian)
    theta.se <- sqrt(diag(theta.vcov))
  }
  ind <- (k-2):k
  alpha <- mod$par[ind]
  Sigma <- a2S(alpha)
  rho <- Sigma[1,2] / sqrt(Sigma[1,1] * Sigma[2,2])
  Sigma.se <- rho.se <- NULL
  if(hessian){
      deltaM <- msm::deltamethod(list( ~ x1 * x1, ~ x1 * x2, ~ x2 * x2 + x3 * x3), mod$par[ind],
                                 theta.vcov[ind, ind])
    Sigma.se <- matrix(c(deltaM[1], deltaM[2], deltaM[2], deltaM[3]), 2, 2)
      rho.se <- msm::deltamethod(list( ~ x1 * x2 / sqrt(x1 * x1 * (x2 * x2 + x3 * x3))), mod$par[ind],
                                 theta.vcov[ind, ind])
   }
 sdrep <- TMB::sdreport(objfit$ADobj, par.fixed = mod$par)
 random <- sdrep$par.random
 random.se <- sqrt(sdrep$diag.cov.random)
 objcrit <- .crit(mod$par, -mod$value, objfit, objfit$model.data$penflag)
 res <- list(theta = mod$par, loglik = objcrit$loglik, AIC = objcrit$AIC, BIC = objcrit$BIC,
             XnS.null = objfit$XnS.null, XnR.null = objfit$XnR.null, lower = NULL,
             seed = myseed,  theta.vcov = theta.vcov, theta.se = theta.se, opt = "ucminf",
             opt.details = mod, ADobj = NULL, model.data = objfit$model.data, sdrep = NULL,
             ranef = random, ranef.se = random.se, Sigma = Sigma, rho = rho,
             Sigma.se = Sigma.se, rho.se = rho.se, M = M,
             penflag = objfit$penflag, penSigma = objfit$penSigma)
 class(res)=c("p2")
 return(res)
}


#' @keywords internal
.recThetaIS <- function(theta, condS, condR, ks, kr, kd, kc, g)
{
   if(condS & condR)
     {
     	gamma1 <- theta[1:ks]
     	gamma2 <- theta[(ks+1):(ks+kr)]
     	delta1 <- theta[(kr+ks+1):(kr+ks+kd)]
     	delta2 <- theta[(kr+ks+kd+1):(kr+ks+kd+kc)]
     	alpha <- theta[(kr+ks+kd+kc+1):(kr+ks+kd+kc+3)]
     	}
   if(condS & !condR)
     {
     	gamma1 <- theta[1:ks]
     	gamma2 <- 0
     	delta1 <- theta[(ks+1):(ks+kd)]
     	delta2 <- theta[(ks+kd+1):(ks+kd+kc)]
     	alpha <- theta[(ks+kd+kc+1):(ks+kd+kc+3)]
     	}
   if(!condS & condR)
     {
     	gamma1 <- 0
     	gamma2 <- theta[1:kr]
     	delta1 <- theta[(kr+1):(kr+kd)]
     	delta2 <- theta[(kr+kd+1):(kr+kd+kc)]
     	alpha <- theta[(kr+kd+kc+1):(kr+kd+kc+3)]
     	}
   if(!condS & !condR)
     {
     	gamma1 <- 0
     	gamma2 <- 0
     	delta1 <- theta[1:kd]
     	delta2 <- theta[(kd+1):(kd+kc)]
     	alpha <- theta[(kd+kc+1):(kd+kc+3)]
     	}
   return(list(gamma1 = gamma1, gamma2 = gamma2, delta1 = delta1, delta2 = delta2,
               alpha = alpha, a = rep(0, g), b = rep(0, g)))
}


#' @keywords internal
.para.fun <- function(theta, lapfit)
{
  data <- lapfit$model.data
  g <- ncol(data$y)
  condS <- !(lapfit$XnS.null)
  condR <- !(lapfit$XnR.null)
  kr <- ncol(data$XnR)
  ks <- ncol(data$XnS)
  kc <- dim(data$XvC)[3]
  kd <- dim(data$XvD)[3]
  para <- .recThetaIS(theta, condS, condR, ks, kr, kd, kc, g)
  map <- list(gamma1 = factor(rep(NA, ks)), gamma2 = factor(rep(NA, kr)), delta1 = factor(rep(NA, kd)),
              delta2 = factor(rep(NA, kc)), alpha = factor(rep(NA, 3)))
  return(list(para = para, map = map))
}


#' @keywords internal
a2S <- function(alpha)
{
  out <- matrix(0,2,2)
  out[1,1] <- alpha[1]^2
  out[1,2] <- out[2,1] <- alpha[1] * alpha[2]
  out[2,2] <- alpha[2]^2 + alpha[3]^2
  return(out)
}

