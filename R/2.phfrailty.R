#' Phase type frailty model
#'
#' Class of objects for univariate phase type frailty models
#'
#' @slot name name of the phase type distribution.
#' @slot bhaz a list comprising of the parameters.
#' @slot coefs regression parameters
#'
#' @return Class object
#' @export
#'
setClass("frailty",
  contains = c("phasetype"),
  slots = list(
    bhaz = "list",
    coefs = "list"
  )
)

#' Constructor function for univariate phase type frailty models
#'
#' @param ph an object of class \linkS4class{phasetype}.
#' @param alpha a probability vector.
#' @param S a sub-intensity matrix.
#' @param structure a valid phase-type structure
#' @param dimension the dimension of the phase-type structure (if provided)
#' @param bhaz baseline hazard function
#' @param bhaz_pars the parameters of the baseline hazard function
#' @param B regression parameters
#'
#' @return An object of class \linkS4class{frailty}.
#' @export
#'
#' @examples
#' frailty(phasetype(structure = "coxian", dimension = 4), bhaz = "weibull", bhaz_pars = 3)
frailty <- function(ph = NULL, bhaz = NULL, bhaz_pars = NULL, B = numeric(0), alpha = NULL, S = NULL, structure = NULL, dimension = 3) {
  if (all(is.null(c(bhaz, bhaz_pars)))) {
    stop("input baseline hazard function and parameters")
  }
  if (is.null(ph)) {
    ph <- phasetype(alpha = alpha, S = S, structure = structure, dimension = dimension)
  }
  if (!bhaz %in% c("exponential", "weibull", "gompertz")) {
    stop("invalid bhaz")
  }
  if (bhaz %in% c("weibull",  "gompertz")) {
    if(is.null(bhaz_pars))bhaz_pars <- 1
    if (length(bhaz_pars) != 1 | sum(bhaz_pars <= 0) > 0) {
      stop("bhaz parameter should be positive and of length one")
    } else {
      names(bhaz_pars) <- "theta"
    }
  }
  if (bhaz %in% c("exponential")) {
    bhaz_pars <- 1
    names(bhaz_pars) <- "theta"
    if(!is.null(bhaz_pars)) {
      warning("exponential only admits value constant equal to one")
    }
  }
  f1 <- function(theta, t) 1
  f2 <- function(theta, t)  theta * t^{theta - 1}
  f3 <- function(theta, t) exp(theta * t)
  nb <- which(bhaz == c("exponential", "weibull", "gompertz"))
  hz <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(theta, t) t
  f2 <- function(theta, t) t^{theta}
  f3 <- function(theta, t) (exp(theta * t) - 1) / theta
  nb <- which(bhaz == c("exponential","weibull", "gompertz"))
  c_hz <- base::eval(parse(text = paste("f", nb, sep = "")))

  f1 <- function(theta, t) t
  f2 <- function(theta, t) t^{1 / theta}
  f3 <- function(theta, t)  log(theta * t + 1) / theta
  nb <- which(bhaz == c("exponential","weibull", "gompertz"))
  c_hz_inv <- base::eval(parse(text = paste("f", nb, sep = "")))

  name <- if(methods::is(ph, "frailty")){ph@name}else{paste("frailty ", ph@name, sep = "")}

  methods::new("frailty",
    name = name,
    pars = ph@pars,
    bhaz = list(name = bhaz, pars = bhaz_pars, hazard = hz,
                cum_hazard = c_hz, cum_hazard_inv = c_hz_inv),
    coefs = list(B = B),
    fit = ph@fit
  )
}


#' Show method for univariate phase type frailty models
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
  cat("coefficients: ", "\n", sep = "")
  print(object@coefs)
})

#' Simulation method for univariate phase type frailty models
#'
#' @param x an object of class \linkS4class{frailty}.
#' @param n an integer of length of realization.
#'
#' @return A realization of independent and identically distributed phase-type frailty variables.
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general"), bhaz = "weibull", bhaz_pars = 2)
#' sim(obj, n = 100)
setMethod("sim", c(x = "frailty"), function(x, n = 1000) {
  theta <- x@bhaz$pars
  c_hz_inv <- x@bhaz$cum_hazard_inv
  U <- c_hz_inv(theta, -log(stats::runif(n)) / rphasetype(n, x@pars$alpha, x@pars$S))
  return(U)
})

#' Density method for univariate phase type frailty models
#'
#' @param x an object of class \linkS4class{frailty}.
#' @param y a vector of locations.
#' @param X a matrix of covariates.
#'
#' @return A list containing the locations and corresponding density evaluations.
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general"), bhaz = "weibull", bhaz_pars = 2)
#' dens(obj, c(1, 2, 3))
setMethod("dens", c(x = "frailty"), function(x, y, X = numeric(0)) {
  theta <- x@bhaz$pars
  fn <- x@bhaz$cum_hazard
  fn_der <- x@bhaz$hazard
  B0 <- x@coefs$B
  y_inf <- (y == Inf)
  dens <- y
  X <- as.matrix(X)
  if (any(dim(X) == 0)) {
    dens[!y_inf] <-  ph_laplace_der_nocons(fn(theta, y[!y_inf]), 2, x@pars$alpha, x@pars$S) * fn_der(theta, y[!y_inf])
    dens[y_inf] <- 0
    return(dens)
  } else {
    if (length(B0) == 0) {
      dens[!y_inf] <-  ph_laplace_der_nocons(fn(theta, y[!y_inf]), 2, x@pars$alpha, x@pars$S) * fn_der(theta, y[!y_inf])
      dens[y_inf] <- 0
      warning("empty regression parameter, returns density without covariate information")
      return(dens)
    } else {
        if (length(B0) != dim(X)[2]) {
          stop("dimension of covariates different from regression parameter")
        } else if (length(y) != dim(X)[1]) {
          stop("dimension of observations different from covariates")
        } else {
          ex <- exp(X%*%B0)
          dens[!y_inf] <-  ph_laplace_der_nocons(fn(theta, y[!y_inf]) * ex, 2, x@pars$alpha, x@pars$S) * fn_der(theta, y[!y_inf]) * ex
          dens[y_inf] <- 0
          return(dens)
        }
      }
  }
})

#' Distribution method for univariate phase type frailty models
#'
#' @param x an object of class \linkS4class{frailty}.
#' @param q a vector of locations.
#' @param X a matrix of covariates
#' @param lower.tail logical parameter specifying whether lower tail (cdf) or upper tail is computed.
#'
#' @return A list containing the locations and corresponding CDF evaluations.
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general"), bhaz = "weibull", bhaz_pars = 2)
#' cdf(obj, c(1, 2, 3))
setMethod("cdf", c(x = "frailty"), function(x, q, X = numeric(0), lower.tail = TRUE) {
  theta <- x@bhaz$pars
  fn <- x@bhaz$cum_hazard
  B0 <- x@coefs$B
  q_inf <- (q == Inf)
  cdf <- q
  X <- as.matrix(X)
  if (any(dim(X) == 0)) {
    if (lower.tail) {
      cdf[!q_inf] <- 1 - ph_laplace(fn(theta, q[!q_inf]), x@pars$alpha, x@pars$S)
    }
    else {
      cdf[!q_inf] <- ph_laplace(fn(theta, q[!q_inf]), x@pars$alpha, x@pars$S)
    }
    cdf[q_inf] <- as.numeric(1 * lower.tail)
    return(cdf)
  } else {
    if (length(B0) == 0) {
      if (lower.tail) {
        cdf[!q_inf] <- 1 - ph_laplace(fn(theta, q[!q_inf]), x@pars$alpha, x@pars$S)
      }
      else {
        cdf[!q_inf] <- ph_laplace(fn(theta, q[!q_inf]), x@pars$alpha, x@pars$S)
      }
      cdf[q_inf] <- as.numeric(1 * lower.tail)
      warning("empty regression parameter, returns cdf without covariate information")
      return(cdf)
    } else {
      if (length(B0) != dim(X)[2]) {
        stop("dimension of covariates different from regression parameter")
      } else if (length(q) != dim(X)[1]) {
        stop("dimension of observations different from covariates")
      } else {
        ex <- exp(X%*%B0)
        if (lower.tail) {
          cdf[!q_inf] <- 1 - ph_laplace(fn(theta, q[!q_inf]) * ex, x@pars$alpha, x@pars$S)
        }
        else {
          cdf[!q_inf] <- ph_laplace(fn(theta, q[!q_inf]) * ex, x@pars$alpha, x@pars$S)
        }
        cdf[q_inf] <- as.numeric(1 * lower.tail)
        return(cdf)
      }
    }
  }
})

#' Hazard rate method for univariate phase type frailty models
#'
#' @param x an object of class \linkS4class{frailty}.
#' @param y a vector of locations.
#' @param X a matrix of covariates.
#'
#' @return A list containing the locations and corresponding hazard rate evaluations.
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general"), bhaz = "weibull", bhaz_pars = 2)
#' haz(obj, c(1, 2, 3))
setMethod("haz", c(x = "frailty"), function(x, y, X = numeric(0)) {
  d <- dens(x, y, X)
  s <- cdf(x, y, X, lower.tail = FALSE)
  return(d / s)
})

#' Coef method for univariate frailty class
#'
#' @param object an object of class \linkS4class{frailty}.
#'
#' @return parameters of frailty model.
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general", dimension = 2), bhaz = "weibull", bhaz_pars = 2)
#' coef(obj)
setMethod("coef", c(object = "frailty"), function(object) {
  L <- object@pars
  L$bhazpars <- object@bhaz$pars
  L$B <- object@coefs$B
  L
})
