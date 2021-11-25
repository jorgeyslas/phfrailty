#' Correlated phase type frailty model
#'
#' Class of objects for correlated phase type frailty models.
#'
#' @slot name Name of the bivariate phase type distribution.
#' @slot bhaz1 A list comprising of the parameters.
#' @slot bhaz2 A list comprising of the parameters.
#' @slot coefs Regression parameters.
#'
#' @return Class object.
#' @export
#'
setClass("correlated",
  contains = c("bivphasetype"),
  slots = list(
    bhaz1 = "list",
    bhaz2 = "list",
    coefs = "list"
  )
)


#' Constructor function for shared phase type frailty models
#'
#' @param bivph An object of class \linkS4class{bivphasetype}.
#' @param alpha A probability vector.
#' @param S11 A sub-intensity matrix.
#' @param S12 A matrix.
#' @param S22 A sub-intensity matrix.
#' @param dimensions The dimension of the bivariate phase-type structure (if provided).
#' @param bhaz1 Baseline hazard function of first marginal.
#' @param bhaz_pars1 The parameters of the baseline hazard function of first marginal.
#' @param bhaz2 Baseline hazard function of first marginal.
#' @param bhaz_pars2 The parameters of the baseline hazard function of first marginal.
#' @param B Regression parameters.
#'
#' @return An object of class \linkS4class{correlated}.
#' @export
#'
#' @examples
#' obj <- bivphasetype(dimensions = c(3, 3))
#' correlated(obj, bhaz1 = "weibull", bhaz_pars1 = 3, bhaz2 = "weibull", bhaz_pars2 = 2)
correlated <- function(bivph = NULL, bhaz1 = NULL, bhaz_pars1 = NULL, bhaz2 = NULL, bhaz_pars2 = NULL, B = numeric(0), alpha = NULL, S11 = NULL, S12 = NULL, S22 = NULL, dimensions = c(3, 3)) {
  if (all(is.null(c(bhaz1, bhaz_pars1, bhaz2, bhaz_pars2)))) {
    stop("input baseline hazard functions and parameters")
  }
  if (is.null(bivph)) {
    ph <- bivphasetype(alpha = alpha, S11 = S11, S12 = S12, S22 = S22, dimensions = dimensions)
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

  name <- if (methods::is(bivph, "frailty")) bivph@name else paste("correlated ", bivph@name, sep = "")

  methods::new("correlated",
    name = name,
    pars = bivph@pars,
    bhaz1 = list(
      name = bhaz1, pars = bhaz_pars1, hazard = hz1,
      cum_hazard = c_hz1, cum_hazard_inv = c_hz_inv1
    ),
    bhaz2 = list(
      name = bhaz2, pars = bhaz_pars2, hazard = hz2,
      cum_hazard = c_hz2, cum_hazard_inv = c_hz_inv2
    ),
    coefs = list(B = B),
    fit = bivph@fit
  )
}


#' Show method for correlated phase type frailty models
#'
#' @param object An object of class \linkS4class{correlated}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "correlated", function(object) {
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


#' Simulation method for correlated phase type frailty models
#'
#' @param x An object of class \linkS4class{correlated}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed correlated phase-type frailty random vectors.
#' @export
#'
#' @examples
#' bivph_obj <- bivphasetype(dimensions = c(3, 3))
#' cobj <- correlated(bivph_obj, bhaz1 = "weibull", bhaz_pars1 = 3, bhaz2 = "weibull", bhaz_pars2 = 2)
#' sim(cobj, n = 100)
setMethod("sim", c(x = "correlated"), function(x, n = 1000) {
  theta1 <- x@bhaz1$pars
  c_hz_inv1 <- x@bhaz1$cum_hazard_inv
  theta2 <- x@bhaz2$pars
  c_hz_inv2 <- x@bhaz2$cum_hazard_inv
  p1_aux <- dim(x@pars$S11)[1]
  p2_aux <- dim(x@pars$S22)[1]
  alpha_aux <- c(x@pars$alpha, rep(0, p2_aux))
  S_aux <- merge_matrices(x@pars$S11, x@pars$S12, x@pars$S22)
  R_aux <- matrix(c(c(rep(1, p1_aux), rep(0, p2_aux)), c(rep(0, p1_aux), rep(1, p2_aux))), ncol = 2)
  rbivph <- rmph(n, alpha_aux, S_aux, R_aux)
  U1 <- c_hz_inv1(theta1, -log(stats::runif(n)) / rbivph[, 1])
  U2 <- c_hz_inv2(theta2, -log(stats::runif(n)) / rbivph[, 2])
  return(matrix(c(U1, U2), ncol = 2))
})


#' Density method for correlated phase type frailty models
#'
#' @param x An object of class \linkS4class{correlated}.
#' @param y A matrix of locations.
#' @param X A matrix of covariates.
#'
#' @return A vector containing the joint density evaluations at the given locations.
#' @export
#'
#' @examples
#' bivph_obj <- bivphasetype(dimensions = c(3,3))
#' cobj <- correlated(bivph_obj, bhaz1 = "weibull", bhaz_pars1 = 3, bhaz2 = "weibull", bhaz_pars2 = 2)
#' dens(cobj, matrix(c(1, 2), ncol = 2))
setMethod("dens", c(x = "correlated"), function(x, y, X = numeric(0)) {
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
    den <- fn1_der(theta1, y[, 1]) * fn2_der(theta2, y[, 2]) * bivph_laplace_der_nocons(matrix(c(fn1(theta1, y[, 1]), fn2(theta2, y[, 2])), ncol = 2), 2, 2, x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22)
    names(den) <- NULL
    return(den)
  } else {
    if (length(B0) == 0) {
      den <- fn1_der(theta1, y[, 1]) * fn2_der(theta2, y[, 2]) * bivph_laplace_der_nocons(matrix(c(fn1(theta1, y[, 1]), fn2(theta2, y[, 2])), ncol = 2), 2, 2, x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22)
      names(den) <- NULL
      warning("empty regression parameter, returns joint density without covariate information")
      return(den)
    } else {
      if (length(B0) != dim(X)[2]) {
        stop("dimension of covariates different from regression parameter")
      } else if (dim(y)[1] != dim(X)[1]) {
        stop("dimension of observations different from covariates")
      } else {
        ex <- exp(X %*% B0)
        den <- ex * fn1_der(theta1, y[, 1]) * fn2_der(theta2, y[, 2]) * bivph_laplace_der_nocons(matrix(c(fn1(theta1, y[, 1]) * ex, fn2(theta2, y[, 2]) * ex), ncol = 2), 2, 2, x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22)
        names(den) <- NULL
        return(den)
      }
    }
  }
})
