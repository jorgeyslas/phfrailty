#' Inhomogeneous Phase Type distributions
#'
#' Class of objects for inhomogeneous phase type distributions
#'
#' @slot name name of the phase type distribution.
#' @slot gfun a list comprising of the parameters.
#' @slot scale scale.
#'
#' @return Class object
#' @export
#'
setClass("phfrailty",
  contains = c("ph"),
  slots = list(
    gfun = "list",
    scale = "numeric"
  )
)

#' Constructor Function for phase type frailty models
#'
#' @param ph An object of class \linkS4class{ph}.
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param structure a valid ph structure
#' @param dimension the dimension of the ph structure (if provided)
#' @param gfun inhomogeneity transform
#' @param gfun_pars the parameters of the inhomogeneity function
#' @param scale scale
#'
#' @return An object of class \linkS4class{phfrailty}.
#' @export
#'
#' @examples
#' phfrailty(ph(structure = "coxian", dimension = 4), gfun = "pareto", gfun_pars = 3)
phfrailty <- function(ph = NULL, gfun = NULL, gfun_pars = NULL, alpha = NULL, S = NULL, structure = NULL, dimension = 3, scale = 1) {
  if (all(is.null(c(gfun, gfun_pars)))) {
    stop("input inhomogeneity function and parameters")
  }
  if (is.null(ph)) {
    ph <- ph(alpha = alpha, S = S, structure = structure, dimension = dimension)
  }
  if (!gfun %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz", "gev", "identity")) {
    stop("invalid gfun")
  }
  if (gfun %in% c("pareto", "weibull", "lognormal", "gompertz")) {
    if(is.null(gfun_pars))gfun_pars <- 1
    if (length(gfun_pars) != 1 | sum(gfun_pars <= 0) > 0) {
      stop("gfun parameter should be positive and of length one")
    } else {
      names(gfun_pars) <- "beta"
    }
  }
  if (gfun %in% c("gev")) {
    if(is.null(gfun_pars))gfun_pars <- c(0, 1, 1)
    if (length(gfun_pars) != 3 | (gfun_pars[2] > 0) == FALSE) {
      stop("gfun parameter should be of length three: mu, sigma, xi, and sigma > 0")
    } else {
      names(gfun_pars) <- c("mu", "sigma", "xi")
    }
  }
  if (gfun %in% c("loglogistic")) {
    if(is.null(gfun_pars))gfun_pars <- c(1, 1)
    if (length(gfun_pars) != 2 | (gfun_pars[1] <= 0) | (gfun_pars[2] <= 0)) {
      stop("gfun parameter should be positive and of length two: alpha, theta > 0")
    } else {
      names(gfun_pars) <- c("alpha", "theta")
    }
  }
  f1 <- function(beta, t) t^{beta}
  f2 <- function(beta, t) log(t / beta + 1)
  f3 <- function(beta, t) log(t + 1)^{beta}
  f4 <- function(beta, t) log((t / beta[1])^{beta[2]} + 1)
  f5 <- function(beta, t) (exp(t * beta) - 1) / beta
  f6 <- function(beta, t, w) reversTransformData(t, w, beta)
  nb <- which(gfun == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  ginv <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(beta, t) t^{beta} * log(t)
  f2 <- function(beta, t) -t/(beta * t + beta^2)
  f3 <- function(beta, t) log(t + 1)^{beta} * log(log(t + 1))
  f4 <- NA
  f5 <- function(beta, t) exp(t * beta) * (t * beta - 1)  / beta^2
  f6 <- NA
  nb <- which(gfun == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  ginv_prime <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(beta, t) beta * t^{beta - 1}
  f2 <- function(beta, t) (t + beta)^{-1}
  f3 <- function(beta, t) beta * log(t + 1)^{beta - 1}/(t + 1)
  f4 <- NA
  f5 <- function(beta, t) exp(t * beta)
  f6 <- NA
  nb <- which(gfun == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  lambda <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(beta, t) t^{beta - 1} + beta * t^{beta - 1} * log(t)
  f2 <- function(beta, t) -(t + beta)^{-2}
  f3 <- function(beta, t) log(t + 1)^{beta - 1}/(t + 1) + beta * log(t + 1)^{beta - 1} * log(log(t + 1))/(t + 1)
  f4 <- NA
  f5 <- function(beta, t) t * exp(t * beta)
  f6 <- NA
  nb <- which(gfun == c("weibull", "pareto", "lognormal", "loglogistic", "gompertz", "gev"))
  lambda_prime <- base::eval(parse(text = paste("f", nb, sep = "")))

  name <- if(is(ph, "iph")){ph@name}else{paste("inhomogeneous ", ph@name, sep = "")}

  methods::new("iph",
    name = name,
    pars = ph@pars,
    gfun = list(name = gfun, pars = gfun_pars, inverse = ginv,
                inverse_prime = ginv_prime, intensity = lambda, intensity_prime = lambda_prime),
    scale = scale,
    fit = ph@fit
  )
}


#' Show Method for inhomogeneous phase type distributions
#'
#' @param object an object of class \linkS4class{iph}.
#' @importFrom methods show
#' @export
#'
setMethod("show", "phfrailty", function(object) {
  cat("object class: ", methods::is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  print(object@pars)
  cat("g-function name: ", object@gfun$name, "\n", sep = "")
  cat("parameters: ", "\n", sep = "")
  methods::show(object@gfun$pars)
})

#' Simulation Method for inhomogeneous phase type distributions
#'
#' @param x an object of class \linkS4class{iph}.
#' @param n an integer of length of realization.
#'
#' @return A realization of independent and identically distributed inhomogeneous phase-type variables.
#' @export
#'
#' @examples
#' obj <- phfrailty(ph(structure = "general"), gfun = "lognormal", gfun_pars = 2)
#' sim(obj, n = 100)
setMethod("sim", c(x = "phfrailty"), function(x, n = 1000) {
  name <- x@gfun$name
  pars <- x@gfun$pars
  scale <- x@scale
  if (name %in% c("pareto", "weibull", "lognormal", "loglogistic", "gompertz")) {
    U <- scale * riph(n, name, x@pars$alpha, x@pars$S, pars)
  }
  if (name %in% c("gev")) {
    U <- scale * rmatrixgev(n, x@pars$alpha, x@pars$S, pars[1], pars[2], pars[3])
  }
  return(U)
})

#' Density Method for inhomogeneous phase type distributions
#'
#' @param x an object of class \linkS4class{phfrailty}.
#' @param y a vector of locations.
#'
#' @return A list containing the locations and corresponding density evaluations.
#' @export
#'
#' @examples
#' obj <- phfrailty(ph(structure = "general"), gfun = "weibull", gfun_pars = 2)
#' dens(obj, c(1, 2, 3))
setMethod("dens", c(x = "phfrailty"), function(x, y) {
  fn <- base::eval(parse(text = paste("m", x@gfun$name, "den", sep = "")))
  scale <- x@scale
  y_inf <- (y == Inf)
  dens <- y
  dens[!y_inf] <- fn(y / scale, x@pars$alpha, x@pars$S, x@gfun$pars) / scale
  dens[y_inf] <- 0
  return(dens)
})

#' Distribution Method for inhomogeneous phase type distributions
#'
#' @param x an object of class \linkS4class{phfrailty}.
#' @param q a vector of locations.
#' @param lower.tail logical parameter specifying whether lower tail (cdf) or upper tail is computed.
#'
#' @return A list containing the locations and corresponding CDF evaluations.
#' @export
#'
setMethod("cdf", c(x = "phfrailty"), function(x,
                                        q,
                                        lower.tail = TRUE) {
  fn <- base::eval(parse(text = paste("m", x@gfun$name, "cdf", sep = "")))
  scale <- x@scale
  q_inf <- (q == Inf)
  cdf <- q
  cdf[!q_inf] <- fn(q[!q_inf] / scale, x@pars$alpha, x@pars$S, x@gfun$pars, lower.tail)
  cdf[q_inf] <- as.numeric(1 * lower.tail)
  return(cdf)
})


#' Coef Method for phfrailty Class
#'
#' @param object an object of class \linkS4class{iph}.
#'
#' @return parameters of phfrailty model.
#' @export
#'
#' @examples
#' obj <- phfrailty(ph(structure = "general", dimension = 2), gfun = "lognormal", gfun_pars = 2)
#' coef(obj)
setMethod("coef", c(object = "phfrailty"), function(object) {
  L <- object@pars
  L$gpars <- object@gfun$pars
  L
})
