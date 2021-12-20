#' Phase-type distributions
#'
#' Class of objects for phase-type distributions.
#'
#' @slot name Name of the phase type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("phasetype",
  slots = list(
    name = "character",
    pars = "list",
    fit = "list"
  ),
  prototype = list(
    name = NA_character_,
    pars = list(),
    fit = list()
  )
)

#' Constructor function for phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S A sub-intensity matrix.
#' @param structure A valid phase-type structure ("general", "coxian", "hyperexponential", "gcoxian", "gerlang").
#' @param dimension The dimension of the phase-type structure (if structure is provided).
#'
#' @return An object of class \linkS4class{phasetype}.
#' @export
#'
#' @examples
#' phasetype(structure = "gcoxian", dimension = 5)
#' phasetype(alpha = c(.5, .5), S = matrix(c(-1, .5, .5, -1), 2, 2))
phasetype <- function(alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }
  if (!is.null(structure)) {
    rs <- random_structure(dimension, structure = structure)
    alpha <- rs[[1]]
    S <- rs[[2]]
    name <- structure
  } else {
    if (dim(S)[1] != dim(S)[2]) {
      stop("matrix S should be square")
    }
    if (length(alpha) != dim(S)[1]) {
      stop("incompatible dimensions")
    }
    name <- "custom"
  }
  methods::new("phasetype",
    name = paste(name, " phasetype(", length(alpha), ")", sep = ""),
    pars = list(alpha = alpha, S = S)
  )
}

#' Moment method for phase-type distributions
#'
#' @param x An object of class \linkS4class{phasetype}.
#' @param k A positive integer (moment order).
#'
#' @return The raw moment of the \linkS4class{phasetype} (or underlying \linkS4class{phasetype}) object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- phasetype(structure = "general", dimension = 3)
#' moment(obj, 2)
setMethod(
  "moment", signature(x = "phasetype"),
  function(x, k = 1) {
    if (k <= 0) {
      return("k should be positive")
    }
    if ((k %% 1) != 0) {
      return("k should be an integer")
    }
    if (methods::is(x, "frailty")) {
      warning("moment of undelying phase-type structure is provided for frailty objects")
    }
    m <- solve(-x@pars$S)
    prod <- diag(nrow(m))
    for (i in 1:k) {
      prod <- prod %*% m
    }
    return(factorial(k) * sum(x@pars$alpha %*% prod))
  }
)

#' Show method for phase type distributions
#'
#' @param object An object of class \linkS4class{phasetype}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "phasetype", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
})

#' Simulation method for phase-type distributions
#'
#' @param x An object of class \linkS4class{phasetype}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed phase-type variables.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "general")
#' sim(obj, n = 100)
setMethod("sim", c(x = "phasetype"), function(x, n = 1000) {
  U <- rphasetype(n, x@pars$alpha, x@pars$S)
  return(U)
})

#' Density method for phase-type distributions
#'
#' @param x An object of class \linkS4class{phasetype}.
#' @param y A vector of locations.
#'
#' @return A vector containing the density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "general")
#' dens(obj, c(1, 2, 3))
setMethod("dens", c(x = "phasetype"), function(x, y) {
  y_inf <- (y == Inf)
  dens <- y
  dens[!y_inf] <- ph_density(y, x@pars$alpha, x@pars$S)
  dens[y_inf] <- 0
  return(dens)
})

#' Distribution method for phase-type distributions
#'
#' @param x An object of class \linkS4class{phasetype}.
#' @param q A vector of locations.
#' @param lower.tail Logical parameter specifying whether lower tail (cdf) or upper tail is computed.
#'
#' @return A vector containing the CDF (or tail) evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "general")
#' cdf(obj, c(1, 2, 3))
setMethod("cdf", c(x = "phasetype"), function(x, q, lower.tail = TRUE) {
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- ph_cdf(q[!q_inf], x@pars$alpha, x@pars$S, lower.tail)
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  return(cdf)
})

#' Survival method for phase-type distributions
#'
#' @param x An object of class \linkS4class{phasetype}.
#' @param q A vector of locations.
#'
#' @return A vector containing the survival function evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "general")
#' surv(obj, c(1, 2, 3))
setMethod("surv", c(x = "phasetype"), function(x, q) {
  cdf(x, q, lower.tail = FALSE)
})

#' Hazard rate method for phase-type distributions
#'
#' @param x An object of class \linkS4class{phasetype}.
#' @param y A vector of locations.
#'
#' @return A vector containing the hazard rate evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "general")
#' haz(obj, c(1, 2, 3))
setMethod("haz", c(x = "phasetype"), function(x, y) {
  d <- dens(x, y)
  s <- cdf(x, y, lower.tail = FALSE)
  return(d / s)
})

#' Quantile method for phase-type distributions
#'
#' @param x An object of class \linkS4class{phasetype}.
#' @param p A vector of probabilities.
#'
#' @return A vector containing the quantile evaluations at the given probabilities.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "general")
#' quan(obj, c(0.5, 0.9, 0.99))
setMethod("quan", c(x = "phasetype"), function(x, p) {
  quan <- numeric(length(p))
  for (i in seq_along(p)) {
    quan[i] <- stats::uniroot(f = function(q) p[i] - cdf(x, 1 / (1 - q) - 1), interval = c(0, 1))$root
  }
  return(1 / (1 - quan) - 1)
})


#' Fit method for phase-type frailty models
#'
#' @param x An object of class \linkS4class{phasetype}.
#' @param y Vector or data.
#' @param rcen Vector of right-censored observations.
#' @param weight Vector of weights for observations.
#' @param X A matrix of covariates.
#' @param initialpoint Initial value for discretization of density.
#' @param truncationpoint Ultimate value for discretization of density.
#' @param maxprobability Max probability allowed for an interval in the discretization.
#' @param maxdelta Max size of interval allowed for the discretization.
#' @param stepsEM Number of EM steps to be performed.
#' @param stepsPH Number of EM steps for the phase-type component at each iteration of the global EM.
#' @param maxit Maximum number of iterations when optimizing g function.
#' @param reltol Relative tolerance when optimizing g function.
#' @param every Number of iterations between likelihood display updates.
#'
#' @return An object of class \linkS4class{phasetype}.
#'
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general", dimension = 2), bhaz = "weibull", bhaz_pars = 2)
#' data <- sim(obj, n = 100)
#' fit(obj, data, stepsEM = 50, every = 10)
setMethod(
  "fit", c(x = "phasetype", y = "ANY"),
  function(x,
           y,
           rcen = numeric(0),
           weight = numeric(0),
           X = numeric(0),
           stepsEM = 100,
           stepsPH = 50,
           initialpoint = 0.0001,
           truncationpoint = 10,
           maxprobability = 0.01,
           maxdelta = 0.05,
           maxit = 100,
           reltol = 1e-8,
           every = 100) {
    if (!all(c(y) > 0)) {
      stop("data should be positive")
    }
    rcenweight <- numeric(0)
    if (length(weight) == 0) {
      weight <- rep(1, length(y))
    } else if (!all(c(weight) >= 0)) {
      stop("weights should be non-negative")
    }
    rcenind <- rcen
    if (length(rcenind) > 0) {
      if (!all((sum(rcenind == 0) + sum(rcenind == 1)) == (length(y)))) {
        stop("right censoring indicator should contain only zeroes and ones")
      }
      if (length(y) != length(rcenind)) {
        stop("data and right censoring indicator should have the same dimensions")
      }
      rcen <- y[(rcenind == 1)]
      y <- y[(rcenind == 0)]
      rcenweight <- weight[(rcenind == 1)]
      weight <- weight[(rcenind == 0)]
    }

    par_haz <- x@bhaz$pars
    chaz <- x@bhaz$cum_hazard
    haz <- x@bhaz$hazard
    X <- as.matrix(X)

    if (any(dim(X) == 0)) {
      LL <- function(alphafn, Sfn, theta, obs, cens, weightobs, weightcens) {
        sum(weightobs * log(ph_laplace_der_nocons(chaz(theta, obs), 2, alphafn, Sfn) * haz(theta, obs))) + sum(weightcens * log(ph_laplace(chaz(theta, cens), alphafn, Sfn)))
      }

      conditional_density <- function(z, alphafn, Sfn, theta, obs, cens, weightobs, weightcens) {
        (sum(z * weightobs * exp(-z * chaz(theta, obs)) * ph_density(z, alphafn, Sfn) / ph_laplace_der_nocons(chaz(theta, obs), 2, alphafn, Sfn)) + sum(weightcens * exp(-z * chaz(theta, cens)) * ph_density(z, alphafn, Sfn) / ph_laplace(chaz(theta, cens), alphafn, Sfn))) / (sum(weight) + sum(rcenweight))
      }

      Ezgiveny <- function(thetamax, alphafn, Sfn, theta, obs, cens, weightobs, weightcens) {
        -sum(weightobs * (log(haz(thetamax, obs)) - chaz(thetamax, obs) * 2 * ph_laplace_der_nocons(chaz(theta, obs), 3, alphafn, Sfn) / ph_laplace_der_nocons(chaz(theta, obs), 2, alphafn, Sfn))) + sum(weightcens * chaz(thetamax, cens) * ph_laplace_der_nocons(chaz(theta, cens), 2, alphafn, Sfn) / ph_laplace(chaz(theta, cens), alphafn, Sfn))
      }

      ph_par <- x@pars
      alpha_fit <- clone_vector(ph_par$alpha)
      S_fit <- clone_matrix(ph_par$S)

      par_haz_fit <- par_haz

      for (k in 1:stepsEM) {
        par_haz_fit <- suppressWarnings(
          stats::optim(
            par = par_haz,
            fn = Ezgiveny,
            theta = par_haz,
            alphafn = alpha_fit,
            Sfn = S_fit,
            obs = y,
            cens = rcen,
            weightobs = weight,
            weightcens = rcenweight,
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
          if (conditional_density(t, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight) < maxprobability / maxdelta) {
            deltat <- maxdelta
          } else {
            deltat <- maxprobability / conditional_density(t, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight)
          }
          proba_aux <- deltat / 6 * (conditional_density(t, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight) + conditional_density(t + deltat, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight))
          while (proba_aux > maxprobability) {
            deltat <- deltat * 0.9
            proba_aux <- deltat / 6 * (conditional_density(t, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight) + conditional_density(t + deltat, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight))
          }
          if (proba_aux > 0) {
            value[j] <- (t * conditional_density(t, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight) + 4 * (t + deltat / 2) * conditional_density(t + deltat / 2, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight) + (t + deltat) * conditional_density(t + deltat, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight)) / (conditional_density(t, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight) + conditional_density(t + deltat, alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight))
            prob[j] <- proba_aux
            j <- j + 1
          }
          t <- t + deltat
        }

        # PH fitting
        for (l in 1:stepsPH) {
          EMstep(alpha_fit, S_fit, value, prob)
        }

        par_haz <- par_haz_fit

        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", LL(alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight),
            sep = " "
          )
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@fit <- list(
        logLik = LL(alpha_fit, S_fit, par_haz, y, rcen, weight, rcenweight),
        nobs = sum(prob)
      )
      x <- frailty(x, bhaz = x@bhaz$name, bhaz_pars = par_haz)
      return(x)
    } else {
      B0 <- x@coefs$B
      h <- dim(X)[2]
      n1 <- length(y)
      n2 <- length(rcen)

      if (n1 + n2 != dim(X)[1]) {
        stop("Number of observations different from number of covariates")
      }

      LL_cov <- function(alphafn, Sfn, theta, obs, cens, scaleobs, scalecens, weightobs, weightcens) {
        sum(weightobs * log(ph_laplace_der_nocons(chaz(theta, obs) * scaleobs, 2, alphafn, Sfn) * haz(theta, obs) * scaleobs)) + sum(weightcens * log(ph_laplace(chaz(theta, cens) * scalecens, alphafn, Sfn)))
      }

      conditional_density_cov <- function(z, alphafn, Sfn, theta, obs, cens, scaleobs, scalecens, weightobs, weightcens) {
        (sum(z * weightobs * exp(-z * chaz(theta, obs) * scaleobs) * ph_density(z, alphafn, Sfn) / ph_laplace_der_nocons(chaz(theta, obs) * scaleobs, 2, alphafn, Sfn)) + sum(weightcens * exp(-z * chaz(theta, cens) * scalecens) * ph_density(z, alphafn, Sfn) / ph_laplace(chaz(theta, cens) * scalecens, alphafn, Sfn))) / (sum(weightobs) + sum(weightcens))
      }

      Ezgiveny_cov <- function(parmax, alphafn, Sfn, theta, obs, cens, scaleobs, scalecens, covinf, weightobs, weightcens) {
        thetamax <- parmax[1:length(theta)]
        Bmax <- parmax[(length(theta) + 1):length(parmax)]
        exmax <- exp(covinf %*% Bmax)
        scaleobsmax <- exmax[1:length(obs)]
        scalecensmax <- utils::tail(exmax, length(cens))
        -sum(weightobs * (log(haz(thetamax, obs) * scaleobsmax) - scaleobsmax * chaz(thetamax, obs) * 2 * ph_laplace_der_nocons(chaz(theta, obs) * scaleobs, 3, alphafn, Sfn) / ph_laplace_der_nocons(chaz(theta, obs) * scaleobs, 2, alphafn, Sfn))) + sum(weightcens * scalecensmax * chaz(thetamax, cens) * ph_laplace_der_nocons(chaz(theta, cens) * scalecens, 2, alphafn, Sfn) / ph_laplace(chaz(theta, cens) * scalecens, alphafn, Sfn))
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

      if (length(rcenind) > 0) {
        X <- base::rbin(X[(rcenind == 0), ], X[(rcenind == 1), ])
      }

      ex <- exp(X %*% B_fit)
      scale1 <- ex[1:length(y)]
      scale2 <- utils::tail(ex, length(rcen))

      for (k in 1:stepsEM) {
        par_fit <- suppressWarnings(
          stats::optim(
            par = c(par_haz, B_fit),
            fn = Ezgiveny_cov,
            theta = par_haz,
            alphafn = alpha_fit,
            Sfn = S_fit,
            obs = y,
            cens = rcen,
            scaleobs = scale1,
            scalecens = scale2,
            covinf = X,
            weightobs = weight,
            weightcens = rcenweight,
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
          if (conditional_density_cov(t, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight) < maxprobability / maxdelta) {
            deltat <- maxdelta
          } else {
            deltat <- maxprobability / conditional_density_cov(t, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight)
          }
          proba_aux <- deltat / 6 * (conditional_density_cov(t, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight) + 4 * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight) + conditional_density_cov(t + deltat, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight))
          while (proba_aux > maxprobability) {
            deltat <- deltat * 0.9
            proba_aux <- deltat / 6 * (conditional_density_cov(t, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight) + 4 * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight) + conditional_density_cov(t + deltat, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight))
          }
          if (proba_aux > 0) {
            value[j] <- (t * conditional_density_cov(t, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight) + 4 * (t + deltat / 2) * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight) + (t + deltat) * conditional_density_cov(t + deltat, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight)) / (conditional_density_cov(t, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight) + 4 * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight) + conditional_density_cov(t + deltat, alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight))
            prob[j] <- proba_aux
            j <- j + 1
          }
          t <- t + deltat
        }

        # PH fitting
        for (l in 1:stepsPH) {
          EMstep(alpha_fit, S_fit, value, prob)
        }

        par_haz <- par_fit[1:length(par_haz)]
        B_fit <- par_fit[(length(par_haz) + 1):length(par_fit)]

        ex <- exp(X %*% B_fit)
        scale1 <- ex[1:length(y)]
        scale2 <- utils::tail(ex, length(rcen))

        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", LL_cov(alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight),
            sep = " "
          )
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@fit <- list(
        logLik = LL_cov(alpha_fit, S_fit, par_haz, y, rcen, scale1, scale2, weight, rcenweight),
        nobs = sum(prob)
      )
      x <- frailty(x, bhaz = x@bhaz$name, bhaz_pars = par_haz, B = B_fit)
      return(x)
    }
  }
)


#' Coef method for phasetype class
#'
#' @param object An object of class \linkS4class{phasetype}.
#'
#' @return Parameters of phase type model.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "general")
#' coef(obj)
setMethod("coef", c(object = "phasetype"), function(object) {
  object@pars
})
