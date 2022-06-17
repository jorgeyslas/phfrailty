#' New generic for simulating matrix distributions
#'
#' Methods are available for objects of classes \linkS4class{phasetype} and
#' \linkS4class{bivphasetype}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return A realization from the matrix distribution.
#'
#'
setGeneric("sim", function(x, ...) {
  standardGeneric("sim")
})

#' New generic for moment of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{phasetype} and
#' \linkS4class{bivphasetype}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Moments from the matrix distribution.
#'
#'
setGeneric("moment", function(x, ...) {
  standardGeneric("moment")
})


#' New generic for the density of matrix distributions
#'
#' Methods are available for objects of classes \linkS4class{phasetype} and
#' \linkS4class{bivphasetype}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Density from the matrix distribution.
#'
#'
setGeneric("dens", function(x, ...) {
  standardGeneric("dens")
})

#' New generic for the hazard rate of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{phasetype}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Hazard rate from the matrix distribution.
#'
#'
setGeneric("haz", function(x, ...) {
  standardGeneric("haz")
})

#' New generic for the distribution function of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{phasetype}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return CDF from the matrix distribution.
#'
#'
setGeneric("cdf", function(x, ...) {
  standardGeneric("cdf")
})


#' New generic for the survival function of matrix distributions
#'
#' Methods are available for objects of classes \linkS4class{phasetype} and
#' \linkS4class{bivphasetype}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Survival function from the matrix distribution.
#'
#'
setGeneric("surv", function(x, ...) {
  standardGeneric("surv")
})


#' New generic for the cross ratio function of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{phasetype}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Cross ratio function from the matrix distribution.
#'
#'
setGeneric("cross", function(x, ...) {
  standardGeneric("cross")
})

#' New generic for the quantile of matrix distributions
#'
#' Methods are available for objects of class \linkS4class{phasetype}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Quantile from the matrix distribution.
#'
#'
setGeneric("quan", function(x, ...) {
  standardGeneric("quan")
})

#' New generic for estimation of matrix distributions
#'
#' Methods are available for objects of classes \linkS4class{phasetype} and
#' \linkS4class{bivphasetype}.
#'
#' @param x An object of the model class.
#' @param y A vector of data.
#' @param ... Further parameters to be passed on.
#'
#' @return An object of the fitted model class.
#'
#'
setGeneric("fit", function(x, y, ...) {
  standardGeneric("fit")
})

#' New generic for the marginals of multivariate matrix distributions
#'
#' Methods are available for objects of classes \linkS4class{phasetype} and
#' \linkS4class{bivphasetype}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Marginal of the matrix distribution.
#'
#'
setGeneric("marginal", function(x, ...) {
  standardGeneric("marginal")
})


#' New generic for the correlation of multivariate matrix distributions
#'
#' Methods are available for objects of classes \linkS4class{phasetype} and
#' \linkS4class{bivphasetype}.
#'
#' @param x An object of the model class.
#' @param ... Further parameters to be passed on.
#'
#' @return Correlation of the multivariate object
#'
#'
setGeneric("corr", function(x, ...) {
  standardGeneric("corr")
})

