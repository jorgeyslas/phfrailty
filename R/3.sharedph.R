#' Shared phase type frailty model
#'
#' Class of objects for shared phase type frailty models
#'
#' @slot name name of the phase type distribution.
#' @slot bhaz1 a list comprising of the parameters.
#' @slot bhaz2 a list comprising of the parameters.
#' @slot coefs regression parameters
#'
#' @return Class object
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
#' @param ph an object of class \linkS4class{phasetype}.
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param structure a valid phase-type structure
#' @param dimension the dimension of the phase-type structure (if provided)
#' @param bhaz1 baseline hazard function of first marginal
#' @param bhaz_pars1 the parameters of the baseline hazard function of first marginal
#' @param bhaz2 baseline hazard function of first marginal
#' @param bhaz_pars2 the parameters of the baseline hazard function of first marginal
#' @param B regression parameters
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
    if(is.null(bhaz_pars1)) bhaz_pars1 <- 1
    if (length(bhaz_pars1) != 1 | sum(bhaz_pars1 <= 0) > 0) {
      stop("bhaz parameter should be positive and of length one")
    } else {
      names(bhaz_pars1) <- "theta"
    }
  } else if (bhaz1 %in% c("gompertz")) {
    if(is.null(bhaz_pars1)) bhaz_pars1 <- c(0.1, 1)
    if (length(bhaz_pars1) != 2 | (bhaz_pars1[1] <= 0) | (bhaz_pars1[2] <= 0)) {
      stop("bhaz parameter should be positive and of length two: a, b > 0")
    } else {
      names(bhaz_pars1) <- c("b", "c")
    }
  }

  if (bhaz2 %in% c("exponential", "weibull")) {
    if(is.null(bhaz_pars2)) bhaz_pars2 <- 1
    if (length(bhaz_pars2) != 1 | sum(bhaz_pars2 <= 0) > 0) {
      stop("bhaz parameter should be positive and of length one")
    } else {
      names(bhaz_pars2) <- "theta"
    }
  } else if (bhaz2 %in% c("gompertz")) {
    if(is.null(bhaz_pars2)) bhaz_pars2 <- c(0.1, 1)
    if (length(bhaz_pars2) != 2 | (bhaz_pars2[1] <= 0) | (bhaz_pars2[2] <= 0)) {
      stop("bhaz parameter should be positive and of length two: a, b > 0")
    } else {
      names(bhaz_pars2) <- c("b", "c")
    }
  }

  f1 <- function(theta, t) theta
  f2 <- function(theta, t)  theta * t^{theta - 1}
  f3 <- function(theta, t) theta[1] * exp(theta[2] * t)
  nb1 <- which(bhaz1 == c("exponential", "weibull", "gompertz"))
  nb2 <- which(bhaz2 == c("exponential", "weibull", "gompertz"))
  hz1 <- base::eval(parse(text = paste("f", nb1, sep = "")))
  hz2 <- base::eval(parse(text = paste("f", nb2, sep = "")))

  f1 <- function(theta, t) theta * t
  f2 <- function(theta, t) t^{theta}
  f3 <- function(theta, t) theta[1] * (exp(theta[2] * t) - 1) / theta[2]
  nb1 <- which(bhaz1 == c("exponential","weibull", "gompertz"))
  nb2 <- which(bhaz2 == c("exponential","weibull", "gompertz"))
  c_hz1 <- base::eval(parse(text = paste("f", nb1, sep = "")))
  c_hz2 <- base::eval(parse(text = paste("f", nb2, sep = "")))

  f1 <- function(theta, t) t / theta
  f2 <- function(theta, t) t^{1 / theta}
  f3 <- function(theta, t)  log(theta[2] * t / theta[1] + 1) / theta[2]
  nb1 <- which(bhaz1 == c("exponential","weibull", "gompertz"))
  nb2 <- which(bhaz2 == c("exponential","weibull", "gompertz"))
  c_hz_inv1 <- base::eval(parse(text = paste("f", nb1, sep = "")))
  c_hz_inv2 <- base::eval(parse(text = paste("f", nb2, sep = "")))

  name <- if(methods::is(ph, "frailty")){ph@name}else{paste("shared ", ph@name, sep = "")}

  methods::new("shared",
               name = name,
               pars = ph@pars,
               bhaz1 = list(name = bhaz1, pars = bhaz_pars1, hazard = hz1,
                           cum_hazard = c_hz1, cum_hazard_inv = c_hz_inv1),
               bhaz2 = list(name = bhaz2, pars = bhaz_pars2, hazard = hz2,
                           cum_hazard = c_hz2, cum_hazard_inv = c_hz_inv2),
               coefs = list(B = B),
               fit = ph@fit
  )
}



#' Show method for shared phase type frailty models
#'
#' @param object an object of class \linkS4class{shared}.
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
#' @param x an object of class \linkS4class{shared}.
#' @param n an integer of length of realization.
#'
#' @return A realization of independent and identically distributed shared phase-type frailty random vectors
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
  return(matrix(c(U1,U2), ncol = 2))
})




#' Coef method for shared frailty class
#'
#' @param object an object of class \linkS4class{shared}.
#'
#' @return parameters of shared phase type frailty model.
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
