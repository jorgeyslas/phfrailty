#' Phase Type distributions
#'
#' Class of objects for phase type distributions
#'
#' @slot name name of the phase type distribution.
#' @slot pars a list comprising of the parameters.
#' @slot fit a list containing estimation information.
#'
#' @return Class object
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

#' Constructor function for phase type distributions
#'
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param structure a valid phase-type structure ("general", "coxian", "hyperexponential", "gcoxian", "gerlang").
#' @param dimension the dimension of the phase-type structure (if structure is provided).
#'
#' @return An object of class \linkS4class{phasetype}.
#' @export
#'
#' @examples
#' phasetype(structure = "gcoxian", dim = 5)
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
#' @param x an object of class \linkS4class{phasetype}.
#' @param k a positive integer (moment order).
#'
#' @return The raw moment of the \linkS4class{phasetype} (or undelying \linkS4class{phasetype}) object.
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- phasetype(structure = "general", dimension = 3)
#' moment(obj, 2)
setMethod("moment", signature(x = "phasetype"),
          function (x, k = 1){
            if(k <= 0) return("k should be positive")
            if((k%%1) != 0) return("k should be an integer")
            if(methods::is(x, "frailty")) warning("moment of undelying phase-type structure is provided for frailty objects")
            m <- solve(-x@pars$S)
            prod <- diag(nrow(m))
            for(i in 1:k){prod <- prod %*% m}
            return(factorial(k)*sum(x@pars$alpha %*% prod))
          }
)

#' Show method for phase type distributions
#'
#' @param object an object of class \linkS4class{phasetype}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "phasetype", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
})

#' Simulation method for phase type distributions
#'
#' @param x an object of class \linkS4class{phasetype}.
#' @param n an integer of length of realization.
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

#' Density method for phase type distributions
#'
#' @param x an object of class \linkS4class{phasetype}.
#' @param y a vector of locations.
#'
#' @return A list containing the locations and corresponding density evaluations.
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

#' Distribution method for phase type distributions
#'
#' @param x an object of class \linkS4class{phasetype}.
#' @param q a vector of locations.
#' @param lower.tail logical parameter specifying whether lower tail (cdf) or upper tail is computed.
#'
#' @return A list containing the locations and corresponding CDF evaluations.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "general")
#' cdf(obj, c(1, 2, 3))
setMethod("cdf", c(x = "phasetype"), function(x,
                                       q,
                                       lower.tail = TRUE) {
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- ph_cdf(q[!q_inf], x@pars$alpha, x@pars$S, lower.tail)
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  return(cdf)
})

#' Hazard rate method for phase type distributions
#'
#' @param x an object of class \linkS4class{phasetype}.
#' @param y a vector of locations.
#'
#' @return A list containing the locations and corresponding hazard rate evaluations.
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

#' Quantile method for phase type distributions
#'
#' @param x an object of class \linkS4class{phasetype}.
#' @param p a vector of probabilities.
#'
#' @return A list containing the probabilities and corresponding quantile evaluations.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "general")
#' quan(obj, c(0.5, 0.9, 0.99))
setMethod("quan", c(x = "phasetype"), function(x,
                                        p) {
  quan <- numeric(length(p))
  for (i in seq_along(p)) {
    quan[i] <- stats::uniroot(f = function(q) p[i] - cdf(x, 1 / (1 - q) - 1), interval = c(0, 1))$root
  }
  return(1 / (1 - quan) - 1)
})


#' Fit method for phase-type frailty models
#'
#' @param x an object of class \linkS4class{phasetype}.
#' @param y vector or data.
#' @param weight vector of weights.
#' @param rcen vector of right-censored observations
#' @param rcenweight vector of weights for right-censored observations.
#' @param initialpoint initial value for discretization of density
#' @param truncationpoint ultimate value for discretization of density
#' @param maxprobability max probability allowed for an interval in the discretization
#' @param maxdelta max size of interval allowed for the discretization
#' @param stepsEM number of EM steps to be performed.
#' @param stepsPH number of EM steps for the phase-type component at each iteration of the global EM.
#' @param maxit maximum number of iterations when optimizing g function.
#' @param reltol relative tolerance when optimizing g function.
#' @param every number of iterations between likelihood display updates.
#'
#' @return An object of class \linkS4class{phasetype}.
#'
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general", dimension = 2), bhaz = "weibull", bhaz_pars = 2)
#' data <- sim(obj, n = 100)
#' fit(obj, data, stepsEM = 1000, every = 200)
setMethod(
  "fit", c(x = "phasetype", y = "ANY"),
  function(x,
           y,
           weight = numeric(0),
           rcen = numeric(0),
           rcenweight = numeric(0),
           stepsEM = 1000,
           stepsPH = 50,
           initialpoint = 0.001,
           truncationpoint = 10,
           maxprobability = 0.01,
           maxdelta = 0.05,
           maxit = 100,
           reltol = 1e-8,
           every = 100) {
    if(!all(c(y, rcen) > 0)) stop("data should be positive")

    par_haz <- x@bhaz$pars
    chaz <- x@bhaz$cum_hazard
    haz <- x@bhaz$hazard

    LL <- function(alpha, S, beta, obs, cens) {
      log(ph_laplace_der_nocons(chaz(beta, obs), 2, alpha, S) * haz(beta, obs)) + log(ph_laplace(chaz(beta, cens), alpha, S))
    }

    conditional_density <- function(z, y, beta, alpha, S) {
      (sum(z * theta * y^(theta-1) * exp(- z * y^theta) * phdensity(z, alpha, S) / mp3density(y, theta, alpha, S)) + sum(exp(- z * cens^theta) * phdensity(z, alpha, S) / (1 - mp3cdf(cens, theta, alpha, S)))) / (length(y) + length(cens))
    }

    ph_par <- x@pars
    alpha_fit <- clone_vector(ph_par$alpha)
    S_fit <- clone_matrix(ph_par$S)

    for (k in 1:stepsEM) {
      A <- discretizate_density(conditional_density, parameters, initialpoint, truncationpoint, maxdelta, maxprobability)
      EMstep(alpha_fit, S_fit, A$values, A$prob)
      if (k %% every == 0) {
        cat("\r", "iteration:", k,
            ", logLik:", LL(alpha_fit, S_fit, par_haz, y, rcen),
            sep = " "
        )
      }
    }
    cat("\n", sep = "")
    x@pars$alpha <- alpha_fit
    x@pars$S <- S_fit
    x@fit <- list(
      logLik = LL(alpha_fit, S_fit, y, rcen),
      nobs = sum(A$prob)
    )

    return(x)
  }
)


#' Coef method for phasetype class
#'
#' @param object an object of class \linkS4class{phasetype}.
#'
#' @return Parameters of phasetype model.
#' @export
#'
#' @examples
#' obj <- phasetype(structure = "general")
#' coef(obj)
setMethod("coef", c(object = "phasetype"), function(object) {
  object@pars
})

