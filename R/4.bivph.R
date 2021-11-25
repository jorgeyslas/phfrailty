#' Bivariate phase type distributions
#'
#' Class of objects for bivariate phase type distributions.
#'
#' @slot name Name of the phase type distribution.
#' @slot pars A list comprising of the parameters.
#' @slot fit A list containing estimation information.
#'
#' @return Class object.
#' @export
#'
setClass("bivphasetype",
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

#' Constructor function for bivariate phase type distributions
#'
#' @param alpha A probability vector.
#' @param S11 A sub-intensity matrix.
#' @param S12 A matrix.
#' @param S22 A sub-intensity matrix.
#' @param dimensions The dimensions of the bivariate phase-type (if no parameters are provided).
#'
#' @return An object of class \linkS4class{bivphasetype}.
#' @export
#'
#' @examples
#' bivphasetype(dimensions = c(3, 3))
#' S11 <- matrix(c(-1, .5, .5, -1), 2, 2)
#' S12 <- matrix(c(.2, .4, .3, .1), 2, 2)
#' S22 <- matrix(c(-2, 0, 1, -1), 2, 2)
#' bivphasetype(alpha = c(.5, .5), S11, S12, S22)
bivphasetype <- function(alpha = NULL, S11 = NULL, S12 = NULL, S22 = NULL, dimensions = c(3, 3)) {
  if (any(is.null(alpha)) & any(is.null(S11))) {
    rs <- random_structure_bivph(dimensions[1], dimensions[2])
    alpha <- rs[[1]]
    S11 <- rs[[2]]
    S12 <- rs[[3]]
    S22 <- rs[[4]]
    name <- "random"
  } else {
    if (dim(S11)[1] != dim(S11)[2]) {
      stop("matrix S11 should be square")
    } else if (dim(S22)[1] != dim(S22)[2]) {
      stop("matrix S22 should be square")
    } else if (dim(S11)[1] != dim(S12)[1]) {
      stop("incompatible dimensions of S11 and S12")
    } else if (dim(S22)[1] != dim(S12)[2]) {
      stop("incompatible dimensions of S12 and S22")
    } else if (length(alpha) != dim(S11)[1]) {
      stop("incompatible dimensions of alpha and S11")
    }
    name <- "custom"
  }
  methods::new("bivphasetype",
    name = paste(name, " bivphasetype(", length(alpha), ")", sep = ""),
    pars = list(alpha = alpha, S11 = S11, S12 = S12, S22 = S22)
  )
}

#' Show method for bivariate phase type distributions
#'
#' @param object An object of class \linkS4class{bivphasetype}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "bivphasetype", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@pars)
})

#' Density method for bivariate phase type distributions
#'
#' @param x An object of class \linkS4class{bivphasetype}.
#' @param y A matrix of locations.
#'
#' @return A vector containing the joint density evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- bivphasetype(dimensions = c(3, 3))
#' dens(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("dens", c(x = "bivphasetype"), function(x, y) {
  dens <- bivph_density(y, x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22)
  return(dens)
})

#' Survival method for bivariate phase type distributions
#'
#' @param x An object of class \linkS4class{bivphasetype}.
#' @param y A matrix of locations.
#'
#' @return A vector containing the joint survival evaluations at the given locations.
#' @export
#'
#' @examples
#' obj <- bivphasetype(dimensions = c(3, 3))
#' surv(obj, matrix(c(0.5, 1), ncol = 2))
setMethod("surv", c(x = "bivphasetype"), function(x, y) {
  sur <- bivph_tail(y, x@pars$alpha, x@pars$S11, x@pars$S12, x@pars$S22)
  return(sur)
})

#' Simulation method for bivariate phase type distributions
#'
#' @param x An object of class \linkS4class{bivphasetype}.
#' @param n An integer of length of realization.
#'
#' @return A realization of independent and identically distributed bivariate
#'  phase-type vector.
#' @export
#'
#' @examples
#' obj <- bivphasetype(dimensions = c(3, 3))
#' sim(obj, n = 100)
setMethod("sim", c(x = "bivphasetype"), function(x, n = 1000) {
  p1_aux <- dim(x@pars$S11)[1]
  p2_aux <- dim(x@pars$S22)[1]
  alpha_aux <- c(x@pars$alpha, rep(0, p2_aux))
  S_aux <- merge_matrices(x@pars$S11, x@pars$S12, x@pars$S22)
  R_aux <- matrix(c(c(rep(1, p1_aux), rep(0, p2_aux)), c(rep(0, p1_aux), rep(1, p2_aux))), ncol = 2)
  U <- rmph(n, alpha_aux, S_aux, R_aux)
  return(U)
})

#' Coef method for bivariate phasetype class
#'
#' @param object An object of class \linkS4class{bivphasetype}.
#'
#' @return Parameters of bivariate phase type model.
#' @export
#'
#' @examples
#' obj <- bivphasetype(dimensions = c(3, 3))
#' coef(obj)
setMethod("coef", c(object = "bivphasetype"), function(object) {
  object@pars
})

#' Marginal method for bivariate phasetype class
#'
#' @param x An object of class \linkS4class{bivphasetype}.
#' @param mar Indicator of which marginal.
#' @return An object of the of class \linkS4class{phasetype}.
#' @export
#'
#' @examples
#' obj <- bivphasetype(dimensions = c(3, 3))
#' marginal(obj, 1)
setMethod("marginal", c(x = "bivphasetype"), function(x, mar = 1) {
  if (mar == 1) {
    xn <- phasetype(alpha = x@pars$alpha, S = x@pars$S11)
  } else {
    alpha0 <- x@pars$alpha %*% base::solve(-x@pars$S11) %*% x@pars$S12
    xn <- phasetype(alpha = alpha0, S = x@pars$S22)
  }
  return(xn)
})
