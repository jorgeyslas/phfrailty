library(phfrailty)

#####################################
# Phase-type Gompertz frailty model #
#####################################
# Data
library(MortalitySmooth)
popDENf <- selectHMDdata("Sweden", "Death", "Females", ages = 1:100, years = 2011)
weight <- as.vector(popDENf)
y <- 1:50
weight <- weight[51:100]

# Initial phase-type frailty model
set.seed(1)
ph_ini <- phasetype(structure = "gcoxian", dimension = 6)
ph_ini@pars$S <- ph_ini@pars$S * 10
b <- 0.05 / 1000
c <- 0.04 * log(10)
phf_ini <- frailty(ph_ini, bhaz = "gompertz", bhaz_pars = c(b, c))

# Fitted phase-type frailty model
ph_fit <- fit(phf_ini, y = y, weight = weight, stepsEM = 200, stepsPH = 99, truncationpoint = 100, maxprobability = 0.005, maxdelta = 0.05, every = 10)
coef(ph_fit)


#####################################
#              Loss data            #
#####################################
# Data
data(loss, package = "copula")
scale_factor <- 10000
y <- loss$loss / scale_factor
rcens <- loss$censored

# Initial phase-type frailty model
set.seed(1)
ph_ini <- phasetype(structure = "coxian", dimension = 4)
phf_ini <- frailty(ph_ini, bhaz = "weibull", bhaz_pars = 2)

# Fitted phase-type frailty model
ph_fit <- fit(phf_ini, y = y, rcen = rcens, stepsEM = 40, stepsPH = 99, truncationpoint = 18, maxprobability = 0.005, maxdelta = 0.05, every = 10)
coef(ph_fit)


#####################################
#  Shared Lognormal frailty model   #
#####################################
# Baseline Gompertz parameters
b1 <- 0.01
c1 <- 1
b2 <- 0.1
c2 <- 2

# Regression parameter
beta <- 0.5

# Simulation of the two groups
# Group 1
set.seed(1)
n1 <- 1000
z <- rlnorm(n1, -0.35, 0.8)
group11 <- (1 / c1) * log(1 - c1 * log(runif(n1)) / (z * b1))
group12 <- (1 / c2) * log(1 - c2 * log(runif(n1)) / (z * b2))
# Group 2
n2 <- 1000
z2 <- rlnorm(n2, -0.35, 0.8) * exp(beta)
group21 <- (1 / c1) * log(1 - c1 * log(runif(n2)) / (z2 * b1))
group22 <- (1 / c2) * log(1 - c2 * log(runif(n2)) / (z2 * b2))
# Data set
y0 <- matrix(c(group11, group21, group12, group22), ncol = 2)
x0 <- c(rep(0, n1), rep(1, n2))

# Initial shared phase-type frailty model
set.seed(1)
ph_ini <- phasetype(structure = "general", dimension = 3)
phs_ini <- shared(ph_ini, bhaz1 = "gompertz", bhaz_pars1 = c(0.01, 1), bhaz2 = "gompertz", bhaz_pars2 = c(0.01, 1), B = 0.5)

# Fitted shared phase-type frailty model
phs_fit <- fit(phs_ini, y = y0, X = x0, stepsEM = 40, stepsPH = 99, truncationpoint = 25, maxprobability = 0.005, maxdelta = 0.05, every = 1)
coef(phs_fit)


######################################
#          Diagnosis data            #
######################################
# Data
data(diagnosis, package = "parfm")
diagnosis$TimeMonths <- diagnosis$Time * 12 / 365.25
y <- matrix(c(diagnosis$TimeMonths[1:106], diagnosis$TimeMonths[107:212]), ncol = 2)

# Initial shared phase-type frailty model
set.seed(1911)
ph_ini <- phasetype(structure = "gcoxian", dimension = 3)
phs_ini <- shared(ph_ini, bhaz1 = "weibull", bhaz_pars1 = c(1, 1), bhaz2 = "weibull", bhaz_pars2 = c(1, 1))

# Fitted shared phase-type frailty model
phs_fit <- fit(phs_ini, y = y, stepsEM = 200, stepsPH = 99, truncationpoint = 25, maxprobability = 0.005, maxdelta = 0.05, every = 1)
coef(phs_fit)
