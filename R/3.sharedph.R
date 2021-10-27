#' Shared phase type frailty model
#'
#' Class of objects for shared phase type frailty models.
#'
#' @slot name Name of the phase type distribution.
#' @slot bhaz1 A list comprising of the parameters.
#' @slot bhaz2 A list comprising of the parameters.
#' @slot coefs Regression parameters.
#'
#' @return Class object.
#' @export
#'
setClass("shared",
  contains = c("phasetype"),
  slots = list(
    bhaz1 = "list",
    bhaz2 = "list",
    coefs = "list"
  )
)



#' Constructor function for shared phase type frailty models
#'
#' @param ph An object of class \linkS4class{phasetype}.
#' @param alpha A probability vector.
#' @param S A sub-intensity matrix.
#' @param structure A valid phase-type structure.
#' @param dimension The dimension of the phase-type structure (if provided).
#' @param bhaz1 Baseline hazard function of first marginal.
#' @param bhaz_pars1 The parameters of the baseline hazard function of first marginal.
#' @param bhaz2 Baseline hazard function of first marginal.
#' @param bhaz_pars2 The parameters of the baseline hazard function of first marginal.
#' @param B Regression parameters.
#'
#' @return An object of class \linkS4class{shared}.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "coxian")
#' shared(obj, bhaz1 = "weibull", bhaz_pars1 = 3, bhaz2 = "weibull", bhaz_pars2 = 2)
shared <- function(ph = NULL, bhaz1 = NULL, bhaz_pars1 = NULL, bhaz2 = NULL, bhaz_pars2 = NULL, B = numeric(0), alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if (all(is.null(c(bhaz1, bhaz_pars1, bhaz2, bhaz_pars2)))) {
    stop("input baseline hazard functions and parameters")
  }
  if (is.null(ph)) {
    ph <- phasetype(alpha = alpha, S = S, structure = structure, dimension = dimension)
  }
  if ((!bhaz1 %in% c("exponential", "weibull", "gompertz")) | (!bhaz2 %in% c("exponential", "weibull", "gompertz"))) {
    stop("invalid bhaz")
  }
  if (bhaz1 %in% c("exponential", "weibull")) {
    if (is.null(bhaz_pars1)) bhaz_pars1 <- 1
    if (length(bhaz_pars1) != 1 | sum(bhaz_pars1 <= 0) > 0) {
      stop("bhaz parameter should be positive and of length one")
    } else {
      names(bhaz_pars1) <- "theta"
    }
  } else if (bhaz1 %in% c("gompertz")) {
    if (is.null(bhaz_pars1)) bhaz_pars1 <- c(0.1, 1)
    if (length(bhaz_pars1) != 2 | (bhaz_pars1[1] <= 0) | (bhaz_pars1[2] <= 0)) {
      stop("bhaz parameter should be positive and of length two: a, b > 0")
    } else {
      names(bhaz_pars1) <- c("b", "c")
    }
  }

  if (bhaz2 %in% c("exponential", "weibull")) {
    if (is.null(bhaz_pars2)) bhaz_pars2 <- 1
    if (length(bhaz_pars2) != 1 | sum(bhaz_pars2 <= 0) > 0) {
      stop("bhaz parameter should be positive and of length one")
    } else {
      names(bhaz_pars2) <- "theta"
    }
  } else if (bhaz2 %in% c("gompertz")) {
    if (is.null(bhaz_pars2)) bhaz_pars2 <- c(0.1, 1)
    if (length(bhaz_pars2) != 2 | (bhaz_pars2[1] <= 0) | (bhaz_pars2[2] <= 0)) {
      stop("bhaz parameter should be positive and of length two: a, b > 0")
    } else {
      names(bhaz_pars2) <- c("b", "c")
    }
  }

  f1 <- function(theta, t) theta
  f2 <- function(theta, t) theta * t^(theta - 1)
  f3 <- function(theta, t) theta[1] * exp(theta[2] * t)
  nb1 <- which(bhaz1 == c("exponential", "weibull", "gompertz"))
  nb2 <- which(bhaz2 == c("exponential", "weibull", "gompertz"))
  hz1 <- base::eval(parse(text = paste("f", nb1, sep = "")))
  hz2 <- base::eval(parse(text = paste("f", nb2, sep = "")))

  f1 <- function(theta, t) theta * t
  f2 <- function(theta, t) t^theta
  f3 <- function(theta, t) theta[1] * (exp(theta[2] * t) - 1) / theta[2]
  nb1 <- which(bhaz1 == c("exponential", "weibull", "gompertz"))
  nb2 <- which(bhaz2 == c("exponential", "weibull", "gompertz"))
  c_hz1 <- base::eval(parse(text = paste("f", nb1, sep = "")))
  c_hz2 <- base::eval(parse(text = paste("f", nb2, sep = "")))

  f1 <- function(theta, t) t / theta
  f2 <- function(theta, t) t^(1 / theta)
  f3 <- function(theta, t) log(theta[2] * t / theta[1] + 1) / theta[2]
  nb1 <- which(bhaz1 == c("exponential", "weibull", "gompertz"))
  nb2 <- which(bhaz2 == c("exponential", "weibull", "gompertz"))
  c_hz_inv1 <- base::eval(parse(text = paste("f", nb1, sep = "")))
  c_hz_inv2 <- base::eval(parse(text = paste("f", nb2, sep = "")))

  name <- if (methods::is(ph, "frailty")) ph@name else paste("shared ", ph@name, sep = "")

  methods::new("shared",
    name = name,
    pars = ph@pars,
    bhaz1 = list(
      name = bhaz1, pars = bhaz_pars1, hazard = hz1,
      cum_hazard = c_hz1, cum_hazard_inv = c_hz_inv1
    ),
    bhaz2 = list(
      name = bhaz2, pars = bhaz_pars2, hazard = hz2,
      cum_hazard = c_hz2, cum_hazard_inv = c_hz_inv2
    ),
    coefs = list(B = B),
    fit = ph@fit
  )
}



#' Show method for shared phase type frailty models
#'
#' @param object An object of class \linkS4class{shared}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "shared", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("b-hazard1 name: ", object@bhaz1$name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@bhaz1$pars)
  cat("b-hazard2 name: ", object@bhaz2$name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@bhaz2$pars)
  cat("coefficients: ", "\n", sep = "")
  print(object@coefs)
})


#' Simulation method for shared phase type frailty models
#'
#' @param x An object of class \linkS4class{shared}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed shared phase-type frailty random vectors.
#' @export
#'
#' @examples
#' ph_obj <- phasetype(structure = "coxian")
#' shared_obj <- shared(ph_obj, bhaz1 = "weibull", bhaz_pars1 = 3, bhaz2 = "weibull", bhaz_pars2 = 2)
#' sim(shared_obj, n = 100)
setMethod("sim", c(x = "shared"), function(x, n = 1000) {
  theta1 <- x@bhaz1$pars
  c_hz_inv1 <- x@bhaz1$cum_hazard_inv
  theta2 <- x@bhaz2$pars
  c_hz_inv2 <- x@bhaz2$cum_hazard_inv
  rph <- rphasetype(n, x@pars$alpha, x@pars$S)
  U1 <- c_hz_inv1(theta1, -log(stats::runif(n)) / rph)
  U2 <- c_hz_inv2(theta2, -log(stats::runif(n)) / rph)
  return(matrix(c(U1, U2), ncol = 2))
})

#' Density method for shared phase type frailty models
#'
#' @param x An object of class \linkS4class{shared}.
#' @param y A matrix of locations.
#' @param X A matrix of covariates.
#'
#' @return A vector containing the corresponding joint density evaluations.
#' @export
#'
#' @examples
#' ph_obj <- phasetype(structure = "coxian")
#' shared_obj <- shared(ph_obj, bhaz1 = "weibull", bhaz_pars1 = 3, bhaz2 = "weibull", bhaz_pars2 = 2)
#' dens(shared_obj, matrix(c(1, 2), ncol = 2))
setMethod("dens", c(x = "shared"), function(x, y, X = numeric(0)) {
  theta1 <- x@bhaz1$pars
  fn1 <- x@bhaz1$cum_hazard
  fn1_der <- x@bhaz1$hazard
  theta2 <- x@bhaz2$pars
  fn2 <- x@bhaz2$cum_hazard
  fn2_der <- x@bhaz2$hazard
  B0 <- x@coefs$B
  y <- as.matrix(y)
  X <- as.matrix(X)

  if (dim(y)[2] != 2) {
    stop("the matrix of locations must have two columns")
  }

  if (any(dim(X) == 0)) {
    den <- 2 * fn1_der(theta1, y[, 1]) * fn2_der(theta2, y[, 2]) * ph_laplace_der_nocons(fn1(theta1, y[, 1]) + fn2(theta2, y[, 2]), 3, x@pars$alpha, x@pars$S)
    return(den)
  } else {
    if (length(B0) == 0) {
      den <- 2 * fn1_der(theta1, y[, 1]) * fn2_der(theta2, y[, 2]) * ph_laplace_der_nocons(fn1(theta1, y[, 1]) + fn2(theta2, y[, 2]), 3, x@pars$alpha, x@pars$S)
      warning("empty regression parameter, returns joint density without covariate information")
      return(den)
    } else {
      if (length(B0) != dim(X)[2]) {
        stop("dimension of covariates different from regression parameter")
      } else if (dim(y)[1] != dim(X)[1]) {
        stop("dimension of observations different from covariates")
      } else {
        ex <- exp(X %*% B0)
        den <- 2 * ex * fn1_der(theta1, y[, 1]) * fn2_der(theta2, y[, 2]) * ph_laplace_der_nocons((fn1(theta1, y[, 1]) + fn2(theta2, y[, 2])) * ex, 3, x@pars$alpha, x@pars$S)
        return(den)
      }
    }
  }
})


#' Survival method for shared phase type frailty models
#'
#' @param x An object of class \linkS4class{shared}.
#' @param q A matrix of locations.
#' @param X A matrix of covariates.
#'
#' @return A vector containing the corresponding joint survival evaluations.
#' @export
#'
#' @examples
#' ph_obj <- phasetype(structure = "coxian")
#' shared_obj <- shared(ph_obj, bhaz1 = "weibull", bhaz_pars1 = 3, bhaz2 = "weibull", bhaz_pars2 = 2)
#' surv(shared_obj, matrix(c(1, 2), ncol = 2))
setMethod("surv", c(x = "shared"), function(x, q, X = numeric(0)) {
  theta1 <- x@bhaz1$pars
  fn1 <- x@bhaz1$cum_hazard
  theta2 <- x@bhaz2$pars
  fn2 <- x@bhaz2$cum_hazard
  B0 <- x@coefs$B
  q <- as.matrix(q)
  X <- as.matrix(X)

  if (dim(q)[2] != 2) {
    stop("the matrix of locations must have two columns")
  }

  if (any(dim(X) == 0)) {
    sur <- ph_laplace(fn1(theta1, q[, 1]) + fn2(theta2, q[, 2]), x@pars$alpha, x@pars$S)
    return(sur)
  } else {
    if (length(B0) == 0) {
      sur <- ph_laplace(fn1(theta1, q[, 1]) + fn2(theta2, q[, 2]), x@pars$alpha, x@pars$S)
      warning("empty regression parameter, returns joint survival without covariate information")
      return(sur)
    } else {
      if (length(B0) != dim(X)[2]) {
        stop("dimension of covariates different from regression parameter")
      } else if (dim(q)[1] != dim(X)[1]) {
        stop("dimension of observations different from covariates")
      } else {
        ex <- exp(X %*% B0)
        sur <- ph_laplace((fn1(theta1, q[, 1]) + fn2(theta2, q[, 2])) * ex, x@pars$alpha, x@pars$S)
        return(sur)
      }
    }
  }
})


#' Fit method for shared phase type frailty models
#'
#' @param x An object of class \linkS4class{shared}.
#' @param y Matrix or data.
#' @param rcen <atrix of indicators of right-censored observations
#' @param X A matrix of covariates.
#' @param initialpoint Initial value for discretization of density.
#' @param truncationpoint Ultimate value for discretization of density.
#' @param maxprobability <ax probability allowed for an interval in the discretization.
#' @param maxdelta Max size of interval allowed for the discretization.
#' @param stepsEM Number of EM steps to be performed.
#' @param stepsPH Number of EM steps for the phase-type component at each iteration of the global EM.
#' @param maxit Maximum number of iterations when optimizing g function.
#' @param reltol Relative tolerance when optimizing g function.
#' @param every Number of iterations between likelihood display updates.
#'
#' @return An object of class \linkS4class{shared}.
#'
#' @export
#'
#' @examples
#' ph_obj <- phasetype(structure = "coxian")
#' shared_obj <- shared(ph_obj, bhaz1 = "weibull", bhaz_pars1 = 3, bhaz2 = "weibull", bhaz_pars2 = 2)
#' data <- sim(shared_obj, n = 100)
#' fit(shared_obj, data, stepsEM = 50, every = 10)
setMethod(
  "fit", c(x = "shared", y = "ANY"),
  function(x,
           y,
           rcen = numeric(0),
           X = numeric(0),
           stepsEM = 1000,
           stepsPH = 50,
           initialpoint = 0.0001,
           truncationpoint = 10,
           maxprobability = 0.01,
           maxdelta = 0.05,
           maxit = 100,
           reltol = 1e-8,
           every = 100) {
    if (!all(y > 0)) {
      stop("data should be positive")
    }

    rcen <- as.matrix(rcen)

    if (dim(rcen)[1] > 0) {
      if (!all((sum(rcen == 0) + sum(rcen == 1)) != (dim(y)[1] * dim(y)[2]))) {
        stop("right censoring indicator should contain only zeroes and ones")
      }
      if (dim(y) != dim(rcen)) {
        stop("data and right censoring indicator should have the same dimensions")
      }

      rowcheck <- rowSums(rcen)
      yaux <- y[(rowcheck == 1), ]
      rcenaux <- rcen[(rowcheck == 1), ]
      rcens10 <- yaux[(rcenaux[, 1] == 1), ]
      rcens01 <- yaux[(rcenaux[, 1] == 0), ]
      rcens11 <- y[(rowcheck == 2), ]
      y <- y[(rowcheck == 0), ]
    } else {
      rcens01 <- matrix(numeric(0), ncol = 2)
      rcens10 <- matrix(numeric(0), ncol = 2)
      rcens11 <- matrix(numeric(0), ncol = 2)
    }
    # 1 indicates censoring  and 0 noncensoring

    par_haz1 <- x@bhaz1$pars
    chaz1 <- x@bhaz1$cum_hazard
    haz1 <- x@bhaz1$hazard

    par_haz2 <- x@bhaz2$pars
    chaz2 <- x@bhaz2$cum_hazard
    haz2 <- x@bhaz2$hazard

    X <- as.matrix(X)

    if (any(dim(X) == 0)) {
      LL <- function(alphafn, Sfn, theta1, theta2, obs, cens01, cens10, cens11) {
        sum(log(2 * ph_laplace_der_nocons(chaz1(theta1, obs[, 1]) + chaz2(theta2, obs[, 2]), 3, alphafn, Sfn) * haz1(theta1, obs[, 1]) * haz2(theta2, obs[, 2]))) + sum(log(ph_laplace_der_nocons(chaz1(theta1, cens01[, 1]) + chaz2(theta2, cens01[, 2]), 2, alphafn, Sfn) * haz1(theta1, cens01[, 1]))) + sum(log(ph_laplace_der_nocons(chaz1(theta1, cens10[, 1]) + chaz2(theta2, cens10[, 2]), 2, alphafn, Sfn) * haz2(theta2, cens10[, 2]))) + sum(log(ph_laplace(chaz1(theta1, cens11[, 1]) + chaz2(theta2, cens11[, 2]), alphafn, Sfn)))
      }

      conditional_density <- function(z, alphafn, Sfn, theta1, theta2, obs, cens01, cens10, cens11) {
        (sum(0.5 * z^(2) * exp(-z * (chaz1(theta1, obs[, 1]) + chaz2(theta2, obs[, 2]))) * ph_density(z, alphafn, Sfn) / ph_laplace_der_nocons(chaz1(theta1, obs[, 1]) + chaz2(theta2, obs[, 2]), 3, alphafn, Sfn))
        + sum(z * exp(-z * (chaz1(theta1, cens01[, 1]) + chaz2(theta2, cens01[, 2]))) * ph_density(z, alphafn, Sfn) / ph_laplace_der_nocons(chaz1(theta1, cens01[, 1]) + chaz2(theta2, cens01[, 2]), 2, alphafn, Sfn))
          + sum(z * exp(-z * (chaz1(theta1, cens10[, 1]) + chaz2(theta2, cens10[, 2]))) * ph_density(z, alphafn, Sfn) / ph_laplace_der_nocons(chaz1(theta1, cens10[, 1]) + chaz2(theta2, cens10[, 2]), 2, alphafn, Sfn))
          + sum(exp(-z * (chaz1(theta1, cens11[, 1]) + chaz2(theta2, cens11[, 2]))) * ph_density(z, alphafn, Sfn) / ph_laplace(chaz1(theta1, cens11[, 1]) + chaz2(theta2, cens11[, 2]), alphafn, Sfn))) / (dim(obs)[1] + dim(cens01)[1] + dim(cens10)[1] + dim(cens11)[1])
      }

      Ezgiveny <- function(parmax, alphafn, Sfn, theta1, theta2, obs, cens01, cens10, cens11) {
        theta1max <- parmax[1:length(theta1)]
        theta2max <- utils::tail(parmax, length(theta2))
        return(-sum(log(haz1(theta1max, obs[, 1])) + log(haz2(theta2max, obs[, 2])) - (chaz1(theta1max, obs[, 1]) + chaz2(theta2max, obs[, 2])) * 3 * ph_laplace_der_nocons(chaz1(theta1, obs[, 1]) + chaz2(theta2, obs[, 2]), 4, alphafn, Sfn) / ph_laplace_der_nocons(chaz1(theta1, obs[, 1]) + chaz2(theta2, obs[, 2]), 3, alphafn, Sfn))
        - sum(log(haz1(theta1max, cens01[, 1])) - (chaz1(theta1max, cens01[, 1]) + chaz2(theta2max, cens01[, 2])) * 2 * ph_laplace_der_nocons(chaz1(theta1, cens01[, 1]) + chaz2(theta2, cens01[, 2]), 3, alphafn, Sfn) / ph_laplace_der_nocons(chaz1(theta1, cens01[, 1]) + chaz2(theta2, cens01[, 2]), 2, alphafn, Sfn))
          - sum(log(haz2(theta2max, cens10[, 2])) - (chaz1(theta1max, cens10[, 1]) + chaz2(theta2max, cens10[, 2])) * 2 * ph_laplace_der_nocons(chaz1(theta1, cens10[, 1]) + chaz2(theta2, cens10[, 2]), 3, alphafn, Sfn) / ph_laplace_der_nocons(chaz1(theta1, cens10[, 1]) + chaz2(theta2, cens10[, 2]), 2, alphafn, Sfn))
          + sum((chaz1(theta1max, cens11[, 1]) + chaz2(theta2max, cens11[, 2])) * ph_laplace_der_nocons(chaz1(theta1, cens11[, 1]) + chaz2(theta2, cens11[, 2]), 2, alphafn, Sfn) / ph_laplace(chaz1(theta1, cens11[, 1]) + chaz2(theta2, cens11[, 2]), alphafn, Sfn)))
      }

      ph_par <- x@pars
      alpha_fit <- clone_vector(ph_par$alpha)
      S_fit <- clone_matrix(ph_par$S)

      for (k in 1:stepsEM) {
        par_haz_fit <- suppressWarnings(
          stats::optim(
            par = c(par_haz1, par_haz2),
            fn = Ezgiveny,
            theta1 = par_haz1,
            theta2 = par_haz2,
            alphafn = alpha_fit,
            Sfn = S_fit,
            obs = y,
            cens01 = rcens01,
            cens10 = rcens10,
            cens11 = rcens11,
            hessian = FALSE,
            control = list(
              maxit = maxit,
              reltol = reltol
            )
          )$par
        )

        # Discretization of density
        deltat <- 0
        t <- initialpoint

        prob <- numeric(0)
        value <- numeric(0)

        j <- 1

        while (t < truncationpoint) {
          if (conditional_density(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11) < maxprobability / maxdelta) {
            deltat <- maxdelta
          } else {
            deltat <- maxprobability / conditional_density(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11)
          }
          proba_aux <- deltat / 6 * (conditional_density(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11) + conditional_density(t + deltat, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11))
          while (proba_aux > maxprobability) {
            deltat <- deltat * 0.9
            proba_aux <- deltat / 6 * (conditional_density(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11) + conditional_density(t + deltat, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11))
          }
          if (proba_aux > 0) {
            value[j] <- (t * conditional_density(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11) + 4 * (t + deltat / 2) * conditional_density(t + deltat / 2, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11) + (t + deltat) * conditional_density(t + deltat, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11)) / (conditional_density(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11) + conditional_density(t + deltat, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11))
            prob[j] <- proba_aux
            j <- j + 1
          }
          t <- t + deltat
        }

        # PH fitting
        for (l in 1:stepsPH) {
          EMstep(alpha_fit, S_fit, value, prob)
        }

        par_haz1 <- par_haz_fit[1:length(par_haz1)]
        par_haz2 <- utils::tail(par_haz_fit, length(par_haz2))

        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", LL(alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11),
            sep = " "
          )
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@fit <- list(
        logLik = LL(alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11),
        nobs = sum(prob)
      )
      x <- shared(x, bhaz1 = x@bhaz1$name, bhaz_pars1 = par_haz1, bhaz2 = x@bhaz2$name, bhaz_pars2 = par_haz2)
      return(x)
    } else {

      B0 <- x@coefs$B
      h <- dim(X)[2]
      n <- length(y[, 1]) + length(rcens01[, 1]) + length(rcens10[, 1]) + length(rcens11[, 1])

      if (n != dim(X)[1]) {
        stop("Number of observations different from number of covariates")
      }


      LL_cov <- function(alphafn, Sfn, theta1, theta2, obs, cens01, cens10, cens11, scal, scal01, scal10, scal11) {
        (sum(log(2 * scal^2 * ph_laplace_der_nocons(scal * (chaz1(theta1, obs[, 1]) + chaz2(theta2, obs[, 2])), 3, alphafn, Sfn) * haz1(theta1, obs[, 1]) * haz2(theta2, obs[, 2])))
         + sum(log(scal01 * ph_laplace_der_nocons(scal01 * (chaz1(theta1, cens01[, 1]) + chaz2(theta2, cens01[, 2])), 2, alphafn, Sfn) * haz1(theta1, cens01[, 1])))
         + sum(log(scal10 * ph_laplace_der_nocons(scal10 * (chaz1(theta1, cens10[, 1]) + chaz2(theta2, cens10[, 2])), 2, alphafn, Sfn) * haz2(theta2, cens10[, 2])))
         + sum(log(ph_laplace(scal11 * (chaz1(theta1, cens11[, 1]) + chaz2(theta2, cens11[, 2])), alphafn, Sfn))))
      }

      conditional_density_cov <- function(z, alphafn, Sfn, theta1, theta2, obs, cens01, cens10, cens11, scal, scal01, scal10, scal11) {
        (sum(0.5 * z^2 * exp(-z * scal * (chaz1(theta1, obs[, 1]) + chaz2(theta2, obs[, 2]))) * ph_density(z, alphafn, Sfn) / ph_laplace_der_nocons(scal * (chaz1(theta1, obs[, 1]) + chaz2(theta2, obs[, 2])), 3, alphafn, Sfn))
         + sum(z * exp(-z * scal01 * (chaz1(theta1, cens01[, 1]) + chaz2(theta2, cens01[, 2]))) * ph_density(z, alphafn, Sfn) / ph_laplace_der_nocons(scal01 * (chaz1(theta1, cens01[, 1]) + chaz2(theta2, cens01[, 2])), 2, alphafn, Sfn))
         + sum(z * exp(-z * scal10 * (chaz1(theta1, cens10[, 1]) + chaz2(theta2, cens10[, 2]))) * ph_density(z, alphafn, Sfn) / ph_laplace_der_nocons(scal10 * (chaz1(theta1, cens10[, 1]) + chaz2(theta2, cens10[, 2])), 2, alphafn, Sfn))
         + sum(exp(-z * scal11 * (chaz1(theta1, cens11[, 1]) + chaz2(theta2, cens11[, 2]))) * ph_density(z, alphafn, Sfn) / ph_laplace(scal11 * (chaz1(theta1, cens11[, 1]) + chaz2(theta2, cens11[, 2])), alphafn, Sfn))) / (dim(obs)[1] + dim(cens01)[1] + dim(cens10)[1] + dim(cens11)[1])
      }

      Ezgiveny_cov <- function(parmax, alphafn, Sfn, theta1, theta2, obs, cens01, cens10, cens11, scal, scal01, scal10, scal11, cinf00, cinf01, cinf10, cinf11) {
        theta1max <- parmax[1:length(theta1)]
        theta2max <- parmax[(length(theta1) + 1):(length(theta1) + length(theta2))]
        Bmax <- parmax[(length(theta1) + length(theta2) + 1):length(parmax)]
        scale00max <- exp(cinf00 %*% Bmax)
        scale01max <- if (length(cinf01) == 0) numeric(0) else exp(cinf01 %*% Bmax)
        scale10max <- if (length(cinf10) == 0) numeric(0) else exp(cinf10 %*% Bmax)
        scale11max <- if (length(cinf11) == 0) numeric(0) else exp(cinf11 %*% Bmax)
        return(-sum(log(scale00max * haz1(theta1max, obs[, 1])) + log(scale00max * haz2(theta2max, obs[, 2])) - scale00max * (chaz1(theta1max, obs[, 1]) + chaz2(theta2max, obs[, 2])) * 3 * ph_laplace_der_nocons(scal * (chaz1(theta1, obs[, 1]) + chaz2(theta2, obs[, 2])), 4, alphafn, Sfn) / ph_laplace_der_nocons(scal * (chaz1(theta1, obs[, 1]) + chaz2(theta2, obs[, 2])), 3, alphafn, Sfn))
               - sum(log(scale01max * haz1(theta1max, cens01[, 1])) - scale01max * (chaz1(theta1max, cens01[, 1]) + chaz2(theta2max, cens01[, 2])) * 2 * ph_laplace_der_nocons(scal01 * (chaz1(theta1, cens01[, 1]) + chaz2(theta2, cens01[, 2])), 3, alphafn, Sfn) / ph_laplace_der_nocons(scal01 * (chaz1(theta1, cens01[, 1]) + chaz2(theta2, cens01[, 2])), 2, alphafn, Sfn))
               - sum(log(scale10max * haz2(theta2max, cens10[, 2])) - scale10max * (chaz1(theta1max, cens10[, 1]) + chaz2(theta2max, cens10[, 2])) * 2 * ph_laplace_der_nocons(scal10 * (chaz1(theta1, cens10[, 1]) + chaz2(theta2, cens10[, 2])), 3, alphafn, Sfn) / ph_laplace_der_nocons(scal10 * (chaz1(theta1, cens10[, 1]) + chaz2(theta2, cens10[, 2])), 2, alphafn, Sfn))
               + sum(scale11max * (chaz1(theta1max, cens11[, 1]) + chaz2(theta2max, cens11[, 2])) * ph_laplace_der_nocons(scal11 * (chaz1(theta1, cens11[, 1]) + chaz2(theta2, cens11[, 2])), 2, alphafn, Sfn) / ph_laplace(scal11 * (chaz1(theta1, cens11[, 1]) + chaz2(theta2, cens11[, 2])), alphafn, Sfn)))
      }

      ph_par <- x@pars
      alpha_fit <- clone_vector(ph_par$alpha)
      S_fit <- clone_matrix(ph_par$S)

      if (length(B0) == 0) {
        B_fit <- rep(0, h)
      } else if (length(B0) != h) {
        B_fit <- rep(0, h)
        warning("Dimension of covariates different from regression parameter. Vector of zeroes used as initial value")
      } else {
        B_fit <- B0
      }


      if (dim(rcen)[1] > 0) {
        rowcheck <- rowSums(rcen)
        Xaux <- X[(rowcheck == 1), ]
        rcenaux <- rcen[(rowcheck == 1), ]
        Xrcens10 <- Xaux[(rcenaux[, 1] == 1), ]
        Xrcens01 <- Xaux[(rcenaux[, 1] == 0), ]
        Xrcens11 <- X[(rowcheck == 2), ]
        X <- X[(rowcheck == 0), ]
      } else {
        Xrcens01 <- matrix(numeric(0), ncol = 2)
        Xrcens10 <- matrix(numeric(0), ncol = 2)
        Xrcens11 <- matrix(numeric(0), ncol = 2)
      }

      scale00 <- exp(X %*% B_fit)
      scale01 <- if (length(Xrcens01) == 0) numeric(0) else exp(Xrcens01 %*% B_fit)
      scale10 <- if (length(Xrcens10) == 0) numeric(0) else exp(Xrcens10 %*% B_fit)
      scale11 <- if (length(Xrcens11) == 0) numeric(0) else exp(Xrcens11 %*% B_fit)

      for (k in 1:stepsEM) {
        par_haz_fit <- suppressWarnings(
          stats::optim(
            par = c(par_haz1, par_haz2, B_fit),
            fn = Ezgiveny_cov,
            theta1 = par_haz1,
            theta2 = par_haz2,
            alphafn = alpha_fit,
            Sfn = S_fit,
            obs = y,
            cens01 = rcens01,
            cens10 = rcens10,
            cens11 = rcens11,
            scal = scale00,
            scal01 = scale01,
            scal10 = scale10,
            scal11 = scale11,
            cinf00 = X,
            cinf01 = Xrcens01,
            cinf10 = Xrcens10,
            cinf11 = Xrcens11,
            hessian = FALSE,
            control = list(
              maxit = maxit,
              reltol = reltol
            )
          )$par
        )

        # Discretization of density
        deltat <- 0
        t <- initialpoint

        prob <- numeric(0)
        value <- numeric(0)

        j <- 1

        while (t < truncationpoint) {
          if (conditional_density_cov(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11) < maxprobability / maxdelta) {
            deltat <- maxdelta
          } else {
            deltat <- maxprobability / conditional_density_cov(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11)
          }
          proba_aux <- deltat / 6 * (conditional_density_cov(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11) + 4 * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11) + conditional_density_cov(t + deltat, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11))
          while (proba_aux > maxprobability) {
            deltat <- deltat * 0.9
            proba_aux <- deltat / 6 * (conditional_density_cov(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11) + 4 * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11) + conditional_density_cov(t + deltat, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11))
          }
          if (proba_aux > 0) {
            value[j] <- (t * conditional_density_cov(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11) + 4 * (t + deltat / 2) * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11) + (t + deltat) * conditional_density_cov(t + deltat, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11)) / (conditional_density_cov(t, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11) + 4 * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11) + conditional_density_cov(t + deltat, alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11))
            prob[j] <- proba_aux
            j <- j + 1
          }
          t <- t + deltat
        }

        # PH fitting
        for (l in 1:stepsPH) {
          EMstep(alpha_fit, S_fit, value, prob)
        }

        par_haz1 <- par_haz_fit[1:length(par_haz1)]
        par_haz2 <- par_haz_fit[(length(par_haz1) + 1):(length(par_haz1) + length(par_haz2))]
        B_fit <- par_haz_fit[(length(par_haz1) + length(par_haz2) + 1):length(par_haz_fit)]

        scale00 <- exp(X %*% B_fit)
        scale01 <- if (length(Xrcens01) == 0) numeric(0) else exp(Xrcens01 %*% B_fit)
        scale10 <- if (length(Xrcens10) == 0) numeric(0) else exp(Xrcens10 %*% B_fit)
        scale11 <- if (length(Xrcens11) == 0) numeric(0) else exp(Xrcens11 %*% B_fit)

        if (k %% every == 0) {
          cat("\r", "iteration:", k,
              ", logLik:", LL_cov(alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11),
              sep = " "
          )
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@fit <- list(
        logLik = LL_cov(alpha_fit, S_fit, par_haz1, par_haz2, y, rcens01, rcens10, rcens11, scale00, scale01, scale10, scale11),
        nobs = sum(prob)
      )
      x <- shared(x, bhaz1 = x@bhaz1$name, bhaz_pars1 = par_haz1, bhaz2 = x@bhaz2$name, bhaz_pars2 = par_haz2, B = B_fit)
      return(x)
    }
  }
)




#' Coef method for shared frailty class
#'
#' @param object An object of class \linkS4class{shared}.
#'
#' @return Parameters of shared phase type frailty model.
#' @export
#'
#' @examples
#' ph_obj <- phasetype(structure = "coxian")
#' shared_obj <- shared(ph_obj, bhaz1 = "weibull", bhaz_pars1 = 3, bhaz2 = "weibull", bhaz_pars2 = 2)
#' coef(shared_obj)
setMethod("coef", c(object = "shared"), function(object) {
  L <- object@pars
  L$bhazpars1 <- object@bhaz1$pars
  L$bhazpars2 <- object@bhaz2$pars
  L$B <- object@coefs$B
  L
})
