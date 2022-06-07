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



