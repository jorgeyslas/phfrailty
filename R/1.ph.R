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

#' Constructor Function for phase type distributions
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

#' Moment Method for phase-type distributions
#'
#' @param x an object of class \linkS4class{phasetype}.
#' @param k a positive integer (moment order).
#'
#' @return The raw moment of the \linkS4class{phasetype} (or undelying \linkS4class{phasetype}) object.
#' @export
#'
#' @examples
#' set.seed(123)
#' ph1 <- phasetype(structure = "general", dimension = 3)
#' moment(ph1, 2)
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

#' Show Method for phase type distributions
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

#' Simulation Method for phase type distributions
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

#' Density Method for phase type distributions
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

#' Distribution Method for phase type distributions
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

#' Hazard rate Method for phase type distributions
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

#' Quantile Method for phase type distributions
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

#' Fit Method for phasetype Class
#'
#' @param x an object of class \linkS4class{phasetype}.
#' @param y vector or data.
#' @param weight vector of weights.
#' @param rcen vector of right-censored observations
#' @param rcenweight vector of weights for right-censored observations.
#' @param stepsEM number of EM steps to be performed.
#' @param methods methods to use for matrix exponential calculation: RM, UNI or PADE
#' @param rkstep Runge-Kutta step size (optional)
#' @param uni_epsilon epsilon parameter for uniformization method
#' @param maxit maximum number of iterations when optimizing g function.
#' @param reltol relative tolerance when optimizing g function.
#' @param every number of iterations between likelihood display updates.
#' @param plot logical indicating whether to plot the fit at each iteration.
#'
#' @return An object of class \linkS4class{phasetype}.
#'
#' @importFrom grDevices dev.off
#' @importFrom graphics hist legend lines
#' @importFrom utils head
#'
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general", dimension = 2), gfun = "weibull", gfun_pars = 2)
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
           methods = c("RK", "RK"),
           rkstep = NA,
           uni_epsilon = NA,
           maxit = 100,
           reltol = 1e-8,
           every = 100,
           plot = FALSE) {
    EMstep <- eval(parse(text = paste("EMstep_", methods[1], sep = "")))
    if(!all(c(y, rcen) > 0)) stop("data should be positive")
    if(!all(c(weight, rcenweight) >= 0)) stop("weights should be non-negative")
    is_iph <- methods::is(x, "frailty")
    if (is_iph) {
      par_g <- x@gfun$pars
      inv_g <- x@gfun$inverse}
    LL <- if(is_iph){eval(parse(text = paste("logLikelihoodM", x@gfun$name, "_", methods[2], sep = "")))
      }else{eval(parse(text = paste("logLikelihoodPH_", methods[2], sep = "")))}
    A <- data_aggregation(y, weight)
    y <- A$un_obs
    weight <- A$weights
    if(length(rcen)>0){
      B <- data_aggregation(rcen, rcenweight)
      rcen <- B$un_obs
      rcenweight <- B$weights
    }
    if(plot == TRUE){
      if(length(rcen)>0) stop("plot option only available for non-censored data")
      h <- hist(rep(y, weight), breaks = 30, plot = FALSE)
      sq <- seq(0, 1.1 * max(y), length.out = 200)
      plot(head(h$breaks, length(h$density)), h$density, col = "#b2df8a",
           main = "Histogram", xlab = "data", ylab = "density", type = "s", lwd = 2)
    }

    ph_par <- x@pars
    alpha_fit <- clone_vector(ph_par$alpha)
    S_fit <- clone_matrix(ph_par$S)

    if (!is_iph) {
      for (k in 1:stepsEM) {
        epsilon1 <- switch(which(methods[1] == c("RK", "UNI","PADE")),
                           if(!is.na(rkstep)){rkstep} else{default_step_length(S_fit)},
                           if(!is.na(uni_epsilon)){uni_epsilon} else{1e-4},
                           0)
        epsilon2 <- switch(which(methods[2] == c("RK", "UNI","PADE")),
                           if(!is.na(rkstep)){rkstep} else{default_step_length(S_fit)},
                           if(!is.na(uni_epsilon)){uni_epsilon} else{1e-4},
                           0)
        EMstep(epsilon1, alpha_fit, S_fit, y, weight, rcen, rcenweight)
        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", LL(epsilon2, alpha_fit, S_fit, y, weight, rcen, rcenweight),
            sep = " "
          )
          if(plot == TRUE){
            dev.off()
            plot(head(h$breaks, length(h$density)), h$density, col = "#b2df8a",
                 main = "Histogram", xlab = "data", ylab = "density", type = "s", lwd = 2)
            tmp_ph <- ph(alpha_fit, S_fit)
            lines(sq, dens(tmp_ph, sq), col = "#33a02c", lwd = 2, lty = 1)
            legend("topright",
                   legend = c("Data", "PH fit"),
                   col = c("#b2df8a", "#33a02c"),
                   lty = c(1,1),
                   bty = "n",
                   lwd = 2,
                   cex = 1.2,
                   text.col = "black",
                   horiz = FALSE,
                   inset = c(0.05, 0.05))
          }
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@fit <- list(
        logLik = LL(epsilon2, alpha_fit, S_fit, y, weight, rcen, rcenweight),
        nobs = sum(A$weights)
      )
    }
    if (is_iph) {
      trans_weight <- weight
      trans_rcenweight <- rcenweight
      for (k in 1:stepsEM) {
        if(x@gfun$name != "gev") {trans <- inv_g(par_g, y); trans_cens <- inv_g(par_g, rcen)
        }else{ t <- inv_g(par_g, y, weight); tc <- inv_g(par_g, rcen, rcenweight)
        trans <- t$obs; trans_weight <- t$weight; trans_cens <- tc$obs; trans_rcenweight <- tc$weight}
        epsilon1 <- switch(which(methods[1] == c("RK", "UNI","PADE")),
                           if(!is.na(rkstep)){rkstep} else{default_step_length(S_fit)},
                           if(!is.na(uni_epsilon)){uni_epsilon} else{1e-4},
                           0)
        epsilon2 <- switch(which(methods[2] == c("RK", "UNI","PADE")),
                           if(!is.na(rkstep)){rkstep} else{default_step_length(S_fit)},
                           if(!is.na(uni_epsilon)){uni_epsilon} else{1e-4},
                           0)
        EMstep(epsilon1, alpha_fit, S_fit, trans, trans_weight, trans_cens, trans_rcenweight)
        opt <- suppressWarnings(
          stats::optim(
            par = par_g,
            fn = LL,
            h = epsilon2,
            alpha = alpha_fit,
            S = S_fit,
            obs = y,
            weight = weight,
            rcens = rcen,
            rcweight = rcenweight,
            hessian = FALSE,
            control = list(
              maxit = maxit,
              reltol = reltol,
              fnscale = -1
            )
          )
        )
        par_g <- opt$par
        if (k %% every == 0) {
          cat("\r", "iteration:", k,
            ", logLik:", opt$value,
            sep = " "
          )
          if(plot == TRUE){
            dev.off()
            plot(head(h$breaks, length(h$density)), h$density, col = "#b2df8a",
                 main = "Histogram", xlab = "data", ylab = "density", type = "s", lwd = 2)
            tmp_ph <- frailty(phasetype(alpha_fit, S_fit), gfun = x@gfun$name, gfun_pars = par_g)
            lines(sq, dens(tmp_ph, sq), col = "#33a02c", lwd = 2, lty = 1)
            legend("topright",
                   legend = c("Data", paste("Matrix-", x@gfun$name," fit", sep ="")),
                   col = c("#b2df8a", "#33a02c"),
                   lty = c(1,1),
                   bty = "n",
                   lwd = 2,
                   cex = 1.2,
                   text.col = "black",
                   horiz = FALSE,
                   inset = c(0.05, 0.05))
          }
        }
      }
      cat("\n", sep = "")
      x@pars$alpha <- alpha_fit
      x@pars$S <- S_fit
      x@fit <- list(
        logLik = opt$value,
        nobs = sum(A$weights)
      )
      x <- frailty(x, gfun = x@gfun$name, gfun_pars = par_g)
    }
    return(x)
  }
)

data_aggregation <- function(y, w) {
  if(length(w) == 0) w <- rep(1, length(y))
  observations <- cbind(y, w)
  mat <- data.frame(observations)
  names(mat) <- c("obs", "weight")
  y <- sort(as.numeric(y))
  un_obs <- unique(y)
  if (length(w) == 0) {
    w <- rep(1, length(y))
  }
  cum_weight <- numeric(0)
  for (i in un_obs) {
    cum_weight <- c(cum_weight, sum(mat$weight[which(mat$obs == i)]))
  }
  return(list(un_obs = un_obs, weights = cum_weight))
}

#' logLik Method for phasetype Class
#'
#' @param object an object of class \linkS4class{phasetype}.
#'
#' @return An object of class logLik.
#' @export
#'
#' @examples
#' obj <- frailty(phasetype(structure = "general", dimension = 2), gfun = "weibull", gfun_pars = 2)
#' data <- sim(obj, n = 100)
#' fitted_ph <- fit(obj, data, stepsEM = 10)
#' logLik(fitted_ph)
setMethod("logLik", "phasetype", function(object) {
  ll <- object@fit$logLik
  attr(ll, "nobs") <- object@fit$nobs
  attr(ll, "df") <- sum(unlist(coef(object)) != 0) - 1
  class(ll) <- "logLik"
  ll
})

#' Coef Method for phasetype Class
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

