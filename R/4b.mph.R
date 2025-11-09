#' Multivariate phase-type distributions
#'
#' Class of objects for multivariate phase-type distributions.
#'
#' @slot name Name of the phase type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("mph",
  slots = list(
    name = "character",
    pars = "list",
    fit = "list"
  )
)

#' Constructor function for multivariate phase-type distributions
#'
#' @param alpha A probability vector.
#' @param S A list of sub-intensity matrices.
#' @param structure A vector of valid ph structures.
#' @param dimension The dimension of the ph structure (if provided).
#' @param variables The dimension of the multivariate phase-type.
#'
#' @return An object of class \linkS4class{mph}.
#' @export
#'
#' @examples
#' mph(structure = c("gcoxian", "general"), dimension = 5)
mph <- function(alpha = NULL, S = NULL, structure = NULL, dimension = 3, variables = NULL) {
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("input a vector and matrix, or a structure")
  }
  if (is.null(variables)) {
    variables <- length(structure)
  }
  if (!any(is.null(structure))) {
    rs <- random_structure(dimension, structure = structure[1])
    alpha <- rs[[1]]
    S <- list()
    S[[1]] <- rs[[2]]
    for (i in 2:variables) {
      S[[i]] <- random_structure(dimension, structure = structure[i])[[2]]
    }
    name <- structure
  } else {
    name <- "custom"
  }
  methods::new("mph",
    name = paste(name, " mph(", length(alpha), ")", sep = " "),
    pars = list(alpha = alpha, S = S)
  )
}

#' Show method for multivariate phase-type distributions
#'
#' @param object An object of class \linkS4class{mph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "mph", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
  cat("number of variables: ", length(object@pars$S), "\n", sep = "")
})


#' Simulation method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param n Length of realization.
#' @param equal_marginals Non-negative integer. If positive, it specifies
#' the number of marginals to simulate from, all from the first matrix.
#'
#' @return A realization of a multivariate phase-type distribution.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' sim(obj, 100)
setMethod("sim", c(x = "mph"), function(x, n = 1000, equal_marginals = 0) {
  if (is.vector(x@pars$alpha)) p <- length(x@pars$alpha)
  if (is.matrix(x@pars$alpha)) p <- ncol(x@pars$alpha)

  if (equal_marginals == 0) {
    d <- length(x@pars$S)
    states <- 1:p
    result <- matrix(NA, n, d)
    for (i in 1:n) {
      state <- sample(states, 1, prob = x@pars$alpha)
      in_vect <- rep(0, p)
      in_vect[state] <- 1
      for (j in 1:d) {
        result[i, j] <- rphasetype(1, in_vect, x@pars$S[[j]])
      }
    }
  } else {
    d <- equal_marginals
    states <- 1:p
    result <- matrix(NA, n, d)
    for (i in 1:n) {
      state <- sample(states, 1, prob = x@pars$alpha)
      in_vect <- rep(0, p)
      in_vect[state] <- 1
      for (j in 1:d) {
        result[i, j] <- rphasetype(1, in_vect, x@pars$S[[1]])
      }
    }
  }
  result
})

#' Marginal method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{ph}.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' marginal(obj, 1)
setMethod("marginal", c(x = "mph"), function(x, mar = 1) {
  if (!(mar %in% 1:length(x@pars$S))) {
    stop("maringal provided not available")
  }
  S <- x@pars$S
  phasetype(alpha = x@pars$alpha, S = S[[mar]])
})

#' Density method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param delta Matrix with right-censoring indicators (1 uncensored, 0 right censored).
#' @param y A matrix of observations.
#'
#' @return A list containing the locations and corresponding density evaluations.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' dens(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("dens", c(x = "mph"), function(x, y) {
  alpha <- x@pars$alpha
  S <- x@pars$S
  d <- length(x@pars$S)

  if (is.matrix(y)) {
    n <- nrow(y)
  }
  if (is.vector(y)) {
    n <- 1
    y <- t(y)
  }

  res <- numeric(n)

  p <- length(x@pars$alpha)

  for (j in 1:p) {
    in_vect <- rep(0, p)
    in_vect[j] <- 1
    aux <- matrix(NA, n, d)
    for (i in 1:d) {
      for (m in 1:n) {
        aux[m, i] <- ph_density(y[m, i], in_vect, S[[i]])
      }
    }
    res <- res + alpha[j] * apply(aux, 1, prod)
  }

  res
})

#' Distribution method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param y A matrix of observations.
#' @param lower.tail Logical parameter specifying whether lower tail (CDF) or
#'  upper tail is computed.
#'
#' @return A list containing the locations and corresponding CDF evaluations.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' cdf(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("cdf", c(x = "mph"), function(x,
                                        y,
                                        lower.tail = TRUE) {
  alpha <- x@pars$alpha
  S <- x@pars$S
  d <- length(x@pars$S)

  if (is.matrix(y)) {
    n <- nrow(y)
  }
  if (is.vector(y)) {
    n <- 1
    y <- t(y)
  }

  res <- numeric(n)

  p <- length(x@pars$alpha)

  for (j in 1:p) {
    in_vect <- rep(0, p)
    in_vect[j] <- 1
    aux <- matrix(NA, n, d)
    for (i in 1:d) {
      aux[, i] <- ph_cdf(y[, i], in_vect, S[[i]], lower.tail)
    }
    res <- res + alpha[j] * apply(aux, 1, prod)
  }

  res
})

#' Moment method for multivariate phase-type distributions
#'
#' @param x An object of class \linkS4class{mph}.
#' @param k A vector of non-negative integer values.
#'
#' @return The corresponding joint moment evaluation.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' moment(obj, c(2, 1))
setMethod("moment", c(x = "mph"), function(x, k) {
  if (all(k == 0) | any(k < 0)) {
    stop("k should be non-negative and not zero")
  }
  if (any((k %% 1) != 0)) {
    stop("k should be an integer")
  }
  if (methods::is(x, "miph")) {
    warning("moment of undelying mph structure is provided for miph objects")
  }
  alpha <- x@pars$alpha
  S <- x@pars$S
  d <- length(x@pars$S)

  p <- length(x@pars$alpha)
  res <- 0

  for (j in 1:p) {
    in_vect <- rep(0, p)
    in_vect[j] <- 1
    aux <- rep(0, d)
    for (i in 1:d) {
      aux[i] <- factorial(k[i]) * sum(in_vect %*% matrix_power(k[i], solve(-S[[i]])))
    }
    res <- res + alpha[j] * prod(aux)
  }

  res
})

#' Correlation method for mph class
#'
#' @param x An object of class \linkS4class{mph}.
#' @param k A vector with the location.
#' @return An real value.
#' @export
#'
#' @examples
#' obj <- mph(structure = c("general", "general"))
#' corr(obj)
setMethod("corr", c(x = "mph"), function(x) {
  m11 <- moment(x, c(1, 0))
  m12 <- moment(x, c(0, 1))
  m21 <- moment(x, c(2, 0))
  m22 <- moment(x, c(0, 2))
  cm <- moment(x, c(1, 1))
  return((cm - m11 * m12) / (sqrt(m21 - m11^2) * sqrt(m22 - m12^2)))
})
