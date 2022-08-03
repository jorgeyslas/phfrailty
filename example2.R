library(phfrailty)
oe_data <- read.csv("oe_lnorm_20211204.csv")
oe_data <- read.csv("oe_lnorm_20220801.csv")

head(oe_data)

oe <- oe_data[oe_data$no == 12 & oe_data$id == "lnorm_ind", ]
oe <- oe_data[oe_data$no == 12 & oe_data$id == "lnorm_neg", ]
head(oe)
dim(oe)
sum(oe$o1)

G <- max(oe$group)

# y2 <- as.matrix(oe[, c(6, 8)])

stepsPH <- 20
stepsEM <- 5


initialpoint1 <- 0.0001
truncationpoint1 <- 8
delta1 <- 0.015
initialpoint2 <- 0.0001
truncationpoint2 <- 8
delta2 <- 0.015

## Initial PH
set.seed(1)
bivph_ini <- bivphasetype(dimensions = c(3, 3))

bivph_par <- bivph_ini@pars
alpha_fit <- clone_vector(bivph_par$alpha)
S11_fit <- clone_matrix(bivph_par$S11)
S12_fit <- clone_matrix(bivph_par$S12)
S22_fit <- clone_matrix(bivph_par$S22)

moment(bivph_ini, c(0, 1))

aggreated_data <- NULL

## Initial regression parameters
# oe_data1 <- oe[oe$e1 > 0, ]
# oe_data2 <- oe[oe$e2 > 0, ]

# Initial values
#m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1)), family = poisson(link = "log"), data = oe_data1)
m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 + 1e-10)), family = poisson(link = "log"), data = oe)
#m2 <- glm(o2 ~ age + offset(log(e2)), family = poisson(link = "log"), data = oe_data2)
m2 <- glm(o2 ~ age + offset(log(e2 + 1e-10)), family = poisson(link = "log"), data = oe)

beta1_fit <- coef(m1)
beta2_fit <- coef(m2)

oe$o1_hat <- fitted(m1)
oe$o2_hat <- fitted(m2)

# oe$o1_hat <- predict.glm(m1, oe, type = "response")
# oe$o2_hat <- predict.glm(m2, oe, type = "response")

# X1 <- cbind(rep(1, nrow(oe)), oe$age, oe$age^2)
# X2 <- X1[, c(1, 2)]
#
# factor1 <- oe$e1 * exp(X1 %*% beta1_fit)
# factor2 <- oe$e2 * exp(X2 %*% beta2_fit)
# fac2 <- cbind(factor1, factor2)
#
# oe$o1_hat <- factor1
# oe$o2_hat <- factor2

aggreated_data <- aggregate(cbind(o1, o1_hat, e2, o2, o2_hat, e1) ~ group + no + id, data = oe, FUN = sum)

# option1
y <- as.matrix(aggreated_data[, c(4, 7)])
fac <- as.matrix(aggreated_data[, c(5, 8)])

for (k in 1:stepsEM) {
  aux <- mp_cor_dens_aux(y, fac, alpha_fit, S11_fit, S12_fit, S22_fit)

  ethetha1 <- (as.vector(y[, 1]) + 1) * (aux$dens_aux1) / (aux$density * fac[, 1])

  ethetha2 <- (as.vector(y[, 2]) + 1) * (aux$dens_aux2) / (aux$density * fac[, 2])

  # option 2
  #fac <- cbind(factor1, factor2)

  #y <- as.matrix(oe[, c(6, 8)])

  #aux <- mp_cor_dens_aux(y, fac, alpha_fit, S11_fit, S12_fit, S22_fit)

  #ethetha1 <- (as.vector(y[, 1]) + 1) * aux$dens_aux1 / aux$density

  #ethetha2 <- (as.vector(y[, 2]) + 1) * aux$dens_aux2 / aux$density

  oe$theta1 <- rep(ethetha1, sapply(1:G, function(g) {
    length(oe$age[oe$group == g])
  }))

  oe$theta2 <- rep(ethetha2, sapply(1:G, function(g) {
    length(oe$age[oe$group == g])
  }))

  # oe_data1$etheta <- rep(ethetha1, sapply(1:G, function(g) {
  #   length(oe_data1$age[oe_data1$group == g])
  # }))
  #
  # oe_data2$etheta <- rep(ethetha2, sapply(1:G, function(g) {
  #   length(oe_data2$age[oe_data2$group == g])
  # }))


  m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 * theta1 + 1e-10)), family = poisson(link = "log"), data = oe)
  m2 <- glm(o2 ~ age + offset(log(e2 * theta2 + 1e-10)), family = poisson(link = "log"), data = oe)

  beta1_fit <- coef(m1)
  beta2_fit <- coef(m2)

  conditional_density <- function(z, alphafn, S11fn, S12fn, S22fn, fac_or, obs, dens_eval) {
    #sum((z[, 1])^obs[, 1] * (z[, 2])^obs[, 2] * exp(-z[, 1] * fac_or[, 1] - z[, 2] * fac_or[, 2]) * bivph_density(z, alphafn, S11fn, S12fn, S22fn) / (factorial(obs[, 1]) * factorial(obs[, 2]) * dens_eval)) / length(obs[, 1])
    sum(exp(obs[, 1] * log(z[, 1] * fac_or[, 1]) + obs[, 2] * log(z[, 2] * fac_or[, 2]) - z[, 1] * fac_or[, 1] - z[, 2] * fac_or[, 2] - lgamma(obs[, 1] + 1) - lgamma(obs[, 2] + 1)) * bivph_density(z, alphafn, S11fn, S12fn, S22fn) / (dens_eval)) / length(obs[, 1])
  }


  # Discretization of density
  prob <- numeric(0)
  value <- base::expand.grid(z1 = seq(initialpoint1, truncationpoint1, by = delta1), z2 = seq(initialpoint2, truncationpoint2, by = delta2))
  for (l in 1:length(value[, 1])) {
    prob[l] <- conditional_density(matrix(c(value[l, 1], value[l, 2]), ncol = 2), alpha_fit, S11_fit, S12_fit, S22_fit, fac, y, aux$density /(fac[, 1] * fac[, 2]) )
  }
  sum(prob * delta1 * delta2)


  # PH fitting
  for (l in 1:stepsPH) {
    EMstep_bivph(alpha_fit, S11_fit, S12_fit, S22_fit, matrix(c(value[, 1], value[, 2]), ncol = 2), prob)
  }

  # factor1 <- oe$e1 * exp(X1 %*% beta1_fit)
  # factor2 <- oe$e2 * exp(X2 %*% beta2_fit)
  # fac2 <- cbind(factor1, factor2)
  #
  # oe$o1_hat <- factor1
  # oe$o2_hat <- factor2

  oe$o1_hat <- fitted(m1)
  oe$o2_hat <- fitted(m2)

  aggreated_data <- aggregate(cbind(o1, o1_hat, e2, o2, o2_hat) ~ group + no + id, data = oe, FUN = sum)

  fac <- as.matrix(aggreated_data[, c(5, 8)])

  #cat(sum(log(fac2^y2))+ sum(log(mp_cor_dens_cov(y, fac, alpha_fit, S11_fit, S12_fit, S22_fit))), sep = " ")

}


bph_fit <- bivphasetype(alpha = alpha_fit, S11 = S11_fit, S12 = S12_fit, S22 = S22_fit)

corr(bph_fit)

corr(bivph_ini)

theta1 <- moment(bph_fit, c(1, 0))
theta2 <- moment(bph_fit, c(0, 1))

hz1 <- function(x) {
  0.5040689 * exp(beta1_fit[1] + beta1_fit[2] * x + beta1_fit[3] * x^2)
}

mu01 <- function(x) {
  0.489635228100584 * exp(-3.2 - 0.025 * x + 0.0006 * x ^ 2)
}


x <- seq(20, 70, by = 1)

plot(x, hz1(x), type = "l")
lines(x, mu01(x), col = "red")

hz2 <- function(x) {
  theta2 * exp(beta2_fit[1] + beta2_fit[2] * x)
}
mu10 <- function(x) {
  exp(0.3 - 0.049 * x)
}

plot(x, hz2(x), type = "l")
lines(x, mu10(x), col = "red")


original_theta <- read.csv("regimes_20211201.csv")

plot(original_theta[original_theta$id == "lnorm_ind", 3], ethetha1)
abline(0, 1, col = "blue")


plot(original_theta[original_theta$id == "lnorm_neg", 3], ethetha1)
abline(0, 1, col = "blue")

plot(original_theta[original_theta$id == "lnorm_neg", 4], ethetha2)
abline(0, 1, col = "blue")

cor(original_theta[original_theta$id == "lnorm_neg", 3], original_theta[original_theta$id == "lnorm_neg", 4])


cor(ethetha1, ethetha2)



sdlog1 <- 0.5122151
meanlog1 <- -0.1311821
z <- seq(0.001, 5, by = 0.001)
plot(z, dlnorm(z, meanlog1, sdlog1), type = "l")
lines(z, dens(marginal(bph_fit, 1), z), col = "red")
lines(z, dgamma(z, 14.59281, 14.59281), col = "blue")

plot(z, dlnorm(z, meanlog1, sdlog1), type = "l")
lines(z, dens(marginal(bph_fit, 2), z), col = "red")


cor(sim(bph_fit, 10000), method = "kendall")
