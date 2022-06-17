#' Phase type mixed Poisson model
#'
#' Class of objects for univariate phase type mixed Poisson models.
#'
#' @slot name Name of the phase type distribution.
#' @slot coefs Regression parameters.
#'
#' @return Class object.
#' @export
#'
setClass("mp",
  contains = c("phasetype"),
  slots = list(
    coefs = "list"
  )
)

#' Constructor function for univariate phase type mixed Poisson models
#'
#' @param ph An object of class \linkS4class{phasetype}.
#' @param alpha A probability vector.
#' @param S A sub-intensity matrix.
#' @param structure A valid phase-type structure.
#' @param dimension The dimension of the phase-type structure (if provided).
#' @param B Regression parameters.
#'
#' @return An object of class \linkS4class{mp}.
#' @export
#'
#' @examples
#' mp(phasetype(structure = "coxian", dimension = 4))
mp <- function(ph = NULL, B = numeric(0), alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if (is.null(ph)) {
    ph <- phasetype(alpha = alpha, S = S, structure = structure, dimension = dimension)
  }

  name <- if (methods::is(ph, "mp")) ph@name else paste("Mixed Poisson ", ph@name, sep = "")

  methods::new("mp",
    name = name,
    pars = ph@pars,
    coefs = list(B = B),
    fit = ph@fit
  )
}

#' Show method for univariate phase type mixed Poisson models
#'
#' @param object An object of class \linkS4class{mp}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "mp", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("coefficients: ", "\n", sep = "")
  print(object@coefs)
})

#' Simulation method for univariate phase type mixed Poisson models
#'
#' @param x An object of class \linkS4class{mp}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed phase-type mixed Poisson variables.
#' @export
#'
#' @examples
#' obj <- mp(phasetype(structure = "general"))
#' sim(obj, n = 100)
setMethod("sim", c(x = "mp"), function(x, n = 1000) {
  vec_rpois <- Vectorize(stats::rpois, vectorize.args = "lambda")
  U <- vec_rpois(1, rphasetype(n, x@pars$alpha, x@pars$S))
  return(U)
})


#' Density method for univariate phase type mixed Poisson models
#'
#' @param x An object of class \linkS4class{mp}.
#' @param y A vector of locations.
#' @param X A matrix of covariates.
#'
#' @return A vector containing the density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- mp(phasetype(structure = "general"))
#' dens(obj, c(1, 2, 3))
setMethod("dens", c(x = "mp"), function(x, y, X = numeric(0)) {
  B0 <- x@coefs$B
  y_inf <- (y == Inf)
  dens <- y
  X <- as.matrix(X)
  if (any(dim(X) == 0)) {
    dens[!y_inf] <- mp_density(y[!y_inf], x@pars$alpha, x@pars$S)
    dens[y_inf] <- 0
    return(dens)
  } else {
    if (length(B0) == 0) {
      dens[!y_inf] <- mp_density(y[!y_inf], x@pars$alpha, x@pars$S)
      dens[y_inf] <- 0
      warning("empty regression parameter, returns density without covariate information")
      return(dens)
    } else {
      if (length(B0) != dim(X)[2]) {
        stop("dimension of covariates different from regression parameter")
      } else if (length(y) != dim(X)[1]) {
        stop("dimension of observations different from covariates")
      } else {
        ex <- exp(X %*% B0)
        dens[!y_inf] <- ex^y[!y_inf] * mp_density_cov(y[!y_inf], ex[!y_inf], x@pars$alpha, x@pars$S)
        dens[y_inf] <- 0
        return(dens)
      }
    }
  }
})



#' Distribution method for univariate phase type mixed Poisson models
#'
#' @param x An object of class \linkS4class{mp}.
#' @param q A vector of locations.
#' @param X A matrix of covariates.
#' @param lower.tail Logical parameter specifying whether lower tail (cdf) or upper tail is computed.
#'
#' @return A vector containing the CDF evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- mp(phasetype(structure = "general"))
#' cdf(obj, c(1, 2, 3))
setMethod("cdf", c(x = "mp"), function(x, q, X = numeric(0), lower.tail = TRUE) {
  aux_fn <- function(n) {
    sum(dens(x, 0:n, X))
  }
  vec_fn <- Vectorize(aux_fn, "n")
  q_inf <- (q == Inf)
  cdf <- q
  X <- as.matrix(X)
  if (lower.tail) {
      cdf[!q_inf] <- vec_fn(q[!q_inf])
  } else {
      cdf[!q_inf] <- 1 - vec_fn(q[!q_inf])
  }
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  return(cdf)
})

#' Fit method for phase-type mixed Poisson models
#'
#' @param x An object of class \linkS4class{mp}.
#' @param y Vector or data.
#' @param exposure Vector of exposures for observations.
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
#' @return An object of class \linkS4class{mp}.
#'
#' @export
#'
#' @examples
#' obj <- mp(phasetype(structure = "general", dimension = 2))
#' data <- sim(obj, n = 100)
#' fit(obj, data, stepsEM = 50, every = 10)
setMethod(
  "fit", c(x = "mp", y = "ANY"),
  function(x,
           y,
           exposure = numeric(0),
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
    if (!all(c(y) >= 0)) {
      stop("data should be non-negative")
    }
    rcenweight <- numeric(0)
    if (length(weight) == 0) {
      weight <- rep(1, length(y))
    } else if (!all(c(weight) >= 0)) {
      stop("weights should be non-negative")
    }

    X <- as.matrix(X)

    if (any(dim(X) == 0)) {
      LL <- function(alphafn, Sfn, obs, weightobs) {
        sum(weightobs * log(mp_density(obs, alphafn, Sfn)))
      }

      conditional_density <- function(z, alphafn, Sfn, obs, weightobs) {
        (sum(weightobs * z^obs * exp(-z) * ph_density(z, alphafn, Sfn) / (factorial(obs) * mp_density(obs, alphafn, Sfn)))) / (sum(weightobs))
      }

      ph_par <- x@pars
      alpha_fit <- clone_vector(ph_par$alpha)
      S_fit <- clone_matrix(ph_par$S)

      for (k in 1:stepsEM) {

        # Discretization of density
        deltat <- 0
        t <- initialpoint

        prob <- numeric(0)
        value <- numeric(0)

        j <- 1

        while (t < truncationpoint) {
          if (conditional_density(t, alpha_fit, S_fit, y, weight) < maxprobability / maxdelta) {
            deltat <- maxdelta
          } else {
            deltat <- maxprobability / conditional_density(t, alpha_fit, S_fit, y, weight)
          }
          proba_aux <- deltat / 6 * (conditional_density(t, alpha_fit, S_fit, y, weight) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, y, weight) + conditional_density(t + deltat, alpha_fit, S_fit, y, weight))
          while (proba_aux > maxprobability) {
            deltat <- deltat * 0.9
            proba_aux <- deltat / 6 * (conditional_density(t, alpha_fit, S_fit, y, weight) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, y, weight) + conditional_density(t + deltat, alpha_fit, S_fit, y, weight))
          }
          if (proba_aux > 0) {
            value[j] <- (t * conditional_density(t, alpha_fit, S_fit, y, weight) + 4 * (t + deltat / 2) * conditional_density(t + deltat / 2, alpha_fit, S_fit, y, weight) + (t + deltat) * conditional_density(t + deltat, alpha_fit, S_fit, y, weight)) / (conditional_density(t, alpha_fit, S_fit, y, weight) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, y, weight) + conditional_density(t + deltat, alpha_fit, S_fit, y, weight))
            prob[j] <- proba_aux
            j <- j + 1
          }
          t <- t + deltat
        }

        # PH fitting
        for (l in 1:stepsPH) {
          EMstep(alpha_fit, S_fit, value, prob)
        }

        if (k %% every == 0) {
          cat("\r", "iteration:", k,
              ", logLik:", LL(alpha_fit, S_fit, y, weight),
              sep = " "
          )
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@fit <- list(
        logLik = LL(alpha_fit, S_fit, y, weight),
        nobs = sum(prob)
      )
      return(x)
    } else {
      B0 <- x@coefs$B
      h <- dim(X)[2]
      n <- length(y)

      if (n != dim(X)[1]) {
        stop("Number of observations different from number of covariates")
      }

      if (length(exposure) == 0) {
        exposure <- rep(1, n)
      } else if (length(exposure) != n) {
        stop("Number of observations different from number of exposures")
      }

      dm <- data.frame(y, X, exposure)

      if (length(B0) == (h + 1)) {
        B_fit <- B0
      } else {
        B_fit <- coef(stats::glm(y ~ . - exposure, offset = log(exposure), weights = weight, family = stats::poisson(), data = dm))
      }

      X0 <- cbind(rep(1, n), X)

      ex <- exp(X0 %*% B_fit)
      mult_fact <- ex * exposure

      LL_cov <- function(alphafn, Sfn, obs, weightobs, m_fact) {
        sum(weightobs * log(m_fact^obs * mp_density_cov(obs, m_fact, alphafn, Sfn)))
      }

      conditional_density_cov <- function(z, alphafn, Sfn, obs, weightobs, m_fact) {
        (sum(weightobs * z^obs * exp(-z * m_fact) * ph_density(z, alphafn, Sfn) / (factorial(obs) * mp_density_cov(obs, m_fact, alphafn, Sfn)))) / (sum(weightobs))
      }


      ph_par <- x@pars
      alpha_fit <- clone_vector(ph_par$alpha)
      S_fit <- clone_matrix(ph_par$S)

      for (k in 1:stepsEM) {

        # Discretization of density
        deltat <- 0
        t <- initialpoint

        prob <- numeric(0)
        value <- numeric(0)

        j <- 1

        while (t < truncationpoint) {
          if (conditional_density_cov(t, alpha_fit, S_fit, y, weight, mult_fact) < maxprobability / maxdelta) {
            deltat <- maxdelta
          } else {
            deltat <- maxprobability / conditional_density_cov(t, alpha_fit, S_fit, y, weight, mult_fact)
          }
          proba_aux <- deltat / 6 * (conditional_density_cov(t, alpha_fit, S_fit, y, weight, mult_fact) + 4 * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, y, weight, mult_fact) + conditional_density_cov(t + deltat, alpha_fit, S_fit, y, weight, mult_fact))
          while (proba_aux > maxprobability) {
            deltat <- deltat * 0.9
            proba_aux <- deltat / 6 * (conditional_density_cov(t, alpha_fit, S_fit, y, weight, mult_fact) + 4 * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, y, weight, mult_fact) + conditional_density_cov(t + deltat, alpha_fit, S_fit, y, weight, mult_fact))
          }
          if (proba_aux > 0) {
            value[j] <- (t * conditional_density_cov(t, alpha_fit, S_fit, y, weight, mult_fact) + 4 * (t + deltat / 2) * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, y, weight, mult_fact) + (t + deltat) * conditional_density_cov(t + deltat, alpha_fit, S_fit, y, weight, mult_fact)) / (conditional_density_cov(t, alpha_fit, S_fit, y, weight, mult_fact) + 4 * conditional_density_cov(t + deltat / 2, alpha_fit, S_fit, y, weight, mult_fact) + conditional_density_cov(t + deltat, alpha_fit, S_fit, y, weight, mult_fact))
            prob[j] <- proba_aux
            j <- j + 1
          }
          t <- t + deltat
        }

        L <- mp_aux_density(y, mult_fact, alpha_fit, S_fit)
        e_theta <- (y + 1) * L$dens_aux / L$density

        dm$e_theta <- e_theta

        B_fit <- coef(stats::glm(y ~ . - exposure - e_theta, offset = log(exposure * e_theta), weights = weight, family = stats::poisson(), data = dm))

        # PH fitting
        for (l in 1:stepsPH) {
          EMstep(alpha_fit, S_fit, value, prob)
        }

        ex <- exp(X0 %*% B_fit)
        mult_fact <- ex * exposure

        if (k %% every == 0) {
          cat("\r", "iteration:", k,
              ", logLik:", LL_cov(alpha_fit, S_fit, y, weight, mult_fact), sum(prob),
              sep = " "
          )
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@coefs$B <- B_fit
      x@fit <- list(
        logLik = LL_cov(alpha_fit, S_fit, y, weight, mult_fact),
        nobs = sum(prob)
      )
      return(x)
    }
  }
)



#' Coef method for univariate phase-type mixed Poisson class
#'
#' @param object An object of class \linkS4class{mp}.
#'
#' @return Parameters of the univariate phase type mixed Poisson model.
#' @export
#'
#' @examples
#' obj <- mp(phasetype(structure = "general", dimension = 2))
#' coef(obj)
setMethod("coef", c(object = "mp"), function(object) {
  L <- object@pars
  L$B <- object@coefs$B
  L
})

