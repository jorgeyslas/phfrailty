#' Phase type frailty model
#'
#' Class of objects for phase type frailty models
#'
#' @slot name name of the phase type distribution.
#' @slot bhaz a list comprising of the parameters.
#' @slot scale scale.
#'
#' @return Class object
#' @export
#'
setClass("frailty",
  contains = c("phasetype"),
  slots = list(
    bhaz = "list",
    scale = "numeric"
  )
)

#' Constructor function for phase type frailty models
#'
#' @param ph an object of class \linkS4class{phasetype}.
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param structure a valid phase-type structure
#' @param dimension the dimension of the phase-type structure (if provided)
#' @param bhaz baseline hazard function
#' @param bhaz_pars the parameters of the baseline hazard function
#' @param scale scale
#'
#' @return An object of class \linkS4class{frailty}.
#' @export
#'
#' @examples
#' frailty(phasetype(structure = "coxian", dimension = 4), bhaz = "pareto", bhaz_pars = 3)
frailty <- function(ph = NULL, bhaz = NULL, bhaz_pars = NULL, alpha = NULL, S = NULL, structure = NULL, dimension = 3, scale = 1) {
  if (all(is.null(c(bhaz, bhaz_pars)))) {
    stop("input inhomogeneity function and parameters")
  }
  if (is.null(ph)) {
    ph <- phasetype(alpha = alpha, S = S, structure = structure, dimension = dimension)
  }
  if (!bhaz %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz", "gev", "identity")) {
    stop("invalid bhaz")
  }
  if (bhaz %in% c("pareto", "weibull", "lognormal", "gompertz")) {
    if(is.null(bhaz_pars))bhaz_pars <- 1
    if (length(bhaz_pars) != 1 | sum(bhaz_pars <= 0) > 0) {
      stop("bhaz parameter should be positive and of length one")
    } else {
      names(bhaz_pars) <- "beta"
    }
  }
  if (bhaz %in% c("gev")) {
    if(is.null(bhaz_pars))bhaz_pars <- c(0, 1, 1)
    if (length(bhaz_pars) != 3 | (bhaz_pars[2] > 0) == FALSE) {
      stop("bhaz parameter should be of length three: mu, sigma, xi, and sigma > 0")
    } else {
      names(bhaz_pars) <- c("mu", "sigma", "xi")
    }
  }
  if (bhaz %in% c("loglogistic")) {
    if(is.null(bhaz_pars))bhaz_pars <- c(1, 1)
    if (length(bhaz_pars) != 2 | (bhaz_pars[1] <= 0) | (bhaz_pars[2] <= 0)) {
      stop("bhaz parameter should be positive and of length two: alpha, theta > 0")
    } else {
      names(bhaz_pars) <- c("alpha", "theta")
    }
  }
  f1 <- function(beta, t) t^{beta}
  f2 <- function(beta, t) log(t / beta + 1)
  f3 <- function(beta, t) log(t + 1)^{beta}
  f4 <- function(beta, t) log((t / beta[1])^{beta[2]} + 1)
  f5 <- function(beta, t) (exp(t * beta) - 1) / beta
  f6 <- function(beta, t, w) reversTransformData(t, w, beta)
  nb <- which(bhaz == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  ginv <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(beta, t) t^{beta} * log(t)
  f2 <- function(beta, t) -t/(beta * t + beta^2)
  f3 <- function(beta, t) log(t + 1)^{beta} * log(log(t + 1))
  f4 <- NA
  f5 <- function(beta, t) exp(t * beta) * (t * beta - 1)  / beta^2
  f6 <- NA
  nb <- which(bhaz == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  ginv_prime <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(beta, t) beta * t^{beta - 1}
  f2 <- function(beta, t) (t + beta)^{-1}
  f3 <- function(beta, t) beta * log(t + 1)^{beta - 1}/(t + 1)
  f4 <- NA
  f5 <- function(beta, t) exp(t * beta)
  f6 <- NA
  nb <- which(bhaz == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  lambda <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(beta, t) t^{beta - 1} + beta * t^{beta - 1} * log(t)
  f2 <- function(beta, t) -(t + beta)^{-2}
  f3 <- function(beta, t) log(t + 1)^{beta - 1}/(t + 1) + beta * log(t + 1)^{beta - 1} * log(log(t + 1))/(t + 1)
  f4 <- NA
  f5 <- function(beta, t) t * exp(t * beta)
  f6 <- NA
  nb <- which(bhaz == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  lambda_prime <- base::eval(parse(text = paste("f", nb, sep = "")))

  name <- if(is(ph, "frailty")){ph@name}else{paste("frailty ", ph@name, sep = "")}

  methods::new("frailty",
    name = name,
    pars = ph@pars,
    bhaz = list(name = bhaz, pars = bhaz_pars, inverse = ginv,
                inverse_prime = ginv_prime, intensity = lambda, intensity_prime = lambda_prime),
    scale = scale,
    fit = ph@fit
  )
}


#' Show method for phase type frailty models
#'
#' @param object an object of class \linkS4class{frailty}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "frailty", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("b-hazard name: ", object@bhaz$name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@bhaz$pars)
})

#' Simulation method for phase type frailty models
#'
#' @param x an object of class \linkS4class{frailty}.
#' @param n an integer of length of realization.
#'
#' @return A realization of independent and identically distributed phase-type frailty variables.
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general"), bhaz = "lognormal", bhaz_pars = 2)
#' sim(obj, n = 100)
setMethod("sim", c(x = "frailty"), function(x, n = 1000) {
  name <- x@bhaz$name
  pars <- x@bhaz$pars
  scale <- x@scale
  if (name %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz")) {
    U <- scale * riph(n, name, x@pars$alpha, x@pars$S, pars)
  }
  if (name %in% c("gev")) {
    U <- scale * rmatrixgev(n, x@pars$alpha, x@pars$S, pars[1], pars[2], pars[3])
  }
  return(U)
})

#' Density method for phase type frailty models
#'
#' @param x an object of class \linkS4class{frailty}.
#' @param y a vector of locations.
#'
#' @return A list containing the locations and corresponding density evaluations.
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general"), bhaz = "weibull", bhaz_pars = 2)
#' dens(obj, c(1, 2, 3))
setMethod("dens", c(x = "frailty"), function(x, y) {
  fn <- base::eval(parse(text = paste("m", x@bhaz$name, "den", sep = "")))
  scale <- x@scale
  y_inf <- (y == Inf)
  dens <- y
  dens[!y_inf] <- fn(y / scale, x@pars$alpha, x@pars$S, x@bhaz$pars) / scale
  dens[y_inf] <- 0
  return(dens)
})

#' Distribution Method for inhomogeneous phase type distributions
#'
#' @param x an object of class \linkS4class{frailty}.
#' @param q a vector of locations.
#' @param lower.tail logical parameter specifying whether lower tail (cdf) or upper tail is computed.
#'
#' @return A list containing the locations and corresponding CDF evaluations.
#' @export
#'
setMethod("cdf", c(x = "frailty"), function(x,
                                        q,
                                        lower.tail = TRUE) {
  fn <- base::eval(parse(text = paste("m", x@bhaz$name, "cdf", sep = "")))
  scale <- x@scale
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- fn(q[!q_inf] / scale, x@pars$alpha, x@pars$S, x@bhaz$pars, lower.tail)
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  return(cdf)
})


#' Coef Method for frailty Class
#'
#' @param object an object of class \linkS4class{frailty}.
#'
#' @return parameters of frailty model.
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general", dimension = 2), bhaz = "lognormal", bhaz_pars = 2)
#' coef(obj)
setMethod("coef", c(object = "frailty"), function(object) {
  L <- object@pars
  L$gpars <- object@bhaz$pars
  L
})
