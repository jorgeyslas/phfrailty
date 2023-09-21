library(phfrailty)
oe_data <- read.csv("oe_gamma_20220802.csv")
oe_data <- read.csv("oe_lnorm_20220802.csv")
oe_data <- read.csv("oe_mgamma_20220812.csv")

head(oe_data)

oe <- oe_data[oe_data$no == 22 & oe_data$id == "lnorm_ind", ]
oe <- oe_data[oe_data$no == 22 & oe_data$id == "lnorm_pos", ]
oe <- oe_data[oe_data$no == 22 & oe_data$id == "lnorm_neg", ]
oe <- oe_data[oe_data$no == 22 & oe_data$id == "gamma_ind", ]
oe <- oe_data[oe_data$no == 22 & oe_data$id == "gamma_pos", ]
oe <- oe_data[oe_data$no == 22 & oe_data$id == "gamma_neg", ]
oe <- oe_data[oe_data$no == 22 & oe_data$id == "mgamma_ind", ]
oe <- oe_data[oe_data$no == 22 & oe_data$id == "mgamma_pos", ]
oe <- oe_data[oe_data$no == 22 & oe_data$id == "mgamma_neg", ]

head(oe)
dim(oe)
sum(oe$o1)

G <- max(oe$group)

# y2 <- as.matrix(oe[, c(6, 8)])

stepsPH <- 20 #50 / 20
stepsEM <- 10 #10 / 5


initialpoint1 <- 0.0001
truncationpoint1 <- 8
delta1 <- 0.015
initialpoint2 <- 0.0001
truncationpoint2 <- 8
delta2 <- 0.015

## Initial PH
set.seed(1)
bivph_ini <- bivphasetype(dimensions = c(5, 5))

bivph_par <- bivph_ini@pars
alpha_fit <- clone_vector(bivph_par$alpha)
S11_fit <- clone_matrix(bivph_par$S11)
S12_fit <- clone_matrix(bivph_par$S12)
S22_fit <- clone_matrix(bivph_par$S22)


#oe_test <- experience_rating_bph(bivph_ini, oe)
#corr(oe_test$bph_fit)

moment(bivph_ini, c(0, 1))
moment(bivph_ini, c(1, 0))

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
  print(sum(prob * delta1 * delta2))


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
  #oe$o1_hat <- exp(beta1_fit[1] + beta1_fit[2] * oe$age + beta1_fit[3] * oe$age^2) * oe$e1
  #oe$o2_hat <- exp(beta2_fit[1] + beta2_fit[2] * oe$age) * oe$e2

  oe$o1_hat <- exp(beta1_fit[1] + beta1_fit[2] * oe$age + beta1_fit[3] * oe$age^2) * (oe$e1 + 1e-15)
  oe$o2_hat <- exp(beta2_fit[1] + beta2_fit[2] * oe$age) * (oe$e2 + 1e-15)

  #oe$o1_hat <- fitted(m1)
  #oe$o2_hat <- fitted(m2)

  aggreated_data <- aggregate(cbind(o1, o1_hat, e2, o2, o2_hat, e1) ~ group + no + id, data = oe, FUN = sum)

  fac <- as.matrix(aggreated_data[, c(5, 8)])

  #cat(sum(log(fac2^y2))+ sum(log(mp_cor_dens_cov(y, fac, alpha_fit, S11_fit, S12_fit, S22_fit))), sep = " ")

  s1 <- sum(lgamma(aggreated_data[, 4] + 1)) + sum(lgamma(aggreated_data[, 7] + 1)) - sum(lgamma(oe[, 6] + 1)) - sum(lgamma(oe[, 8] + 1))
  s2 <- sum(oe[, 6]* log(oe[, 9])) + sum(oe[, 8] * log(oe[, 10]))
  s3 <- sum((aggreated_data[, 4]) * log(aggreated_data[, 5])) + sum((aggreated_data[, 7]) * log(aggreated_data[, 8]))
  s4 <- sum(log(mp_cor_dens_cov(as.matrix(aggreated_data[, c(4, 7)]), as.matrix(aggreated_data[, c(5, 8)]), alpha_fit, S11_fit, S12_fit, S22_fit) / (aggreated_data[, 5] * aggreated_data[, 8])))
  cat(s1 + s2 - s3 + s4, sep = " ")

}

library(ggplot2)

p.samp <- cbind(value, prob )
contour_bivariate_fitted<-ggplot(p.samp, aes(x=z1, y=z2, z=prob)) + geom_contour(aes(color=..level..), bins=80)
contour_bivariate_fitted

esm1 <- sum(value[, 1] * (prob * delta1 * delta2))
esm2 <- sum(value[, 2] * (prob * delta1 * delta2))
esm12 <- sum(value[, 1]^2 * (prob * delta1 * delta2))
esm22 <- sum(value[, 2]^2 * (prob * delta1 * delta2))
esc <- sum(value[, 1] * value[, 2] * (prob * delta1 * delta2))

esc - esm1 * esm2
(esc - esm1 * esm2) / (sqrt(esm12 - esm1^2) * sqrt(esm22 - esm2^2))

bph_fit <- bivphasetype(alpha = alpha_fit, S11 = S11_fit, S12 = S12_fit, S22 = S22_fit)

corr(bph_fit)

corr(bivph_ini)

mean1 <- moment(bph_fit, c(1, 0))
mean1
mean2 <- moment(bph_fit, c(0, 1))
mean2

var1 <- moment(bph_fit, c(2, 0)) - mean1^2
var1
var2 <- moment(bph_fit, c(0, 2)) - mean2^2
var2

sqrt(var1) / mean1
sqrt(var2) / mean2

sqrt(3 / 10)



oe_out <- aggreated_data[, 1:8]
oe_out <- oe_out[, -6]

oe_out$theta1 <- ethetha1
oe_out$theta2 <- ethetha2

oe_out$psi1 <- rep(NA, length(ethetha1))
oe_out$psi2 <- rep(NA, length(ethetha1))

# Intercept
oe_out$int1 <- beta1_fit[1]
oe_out$int2 <- beta2_fit[1]

# (log) Linear term
oe_out$lin1 <- beta1_fit[2]
oe_out$lin2 <- beta2_fit[2]

# (log) Second order term
oe_out$sec1 <- beta1_fit[3]



write.csv(oe_out, "phfit_gammapos_group22.csv", row.names = FALSE)
write.csv(oe_out, "phfit_lnormpos_group22.csv", row.names = FALSE)

write.csv(oe_out, "phfit_gammaneg_group22.csv", row.names = FALSE)

hz1 <- function(x) {
  exp(beta1_fit[1] + beta1_fit[2] * x + beta1_fit[3] * x^2)
}

mu01 <- function(x) {
   #exp(-3.2 - 0.025 * x + 0.0006 * x ^ 2)
  exp(-4.5 - 0.018 * x + 0.00064 * x ^ 2)
}


x <- seq(20, 67, by = 1)
plot(x, hz1(x), type = "l")
lines(x, mu01(x), col = "red")

hz2 <- function(x) {
  exp(beta2_fit[1] + beta2_fit[2] * x)
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


plot(original_theta[original_theta$id == "gamma_ind", 3], ethetha1)
abline(0, 1, col = "blue")

plot(original_theta[original_theta$id == "gamma_ind", 4], ethetha2)
abline(0, 1, col = "blue")

oe_data_aux <- oe_data
oe_data_aux$real_fac1 <- mu01(oe_data$age) * oe_data$e1
oe_data_aux$est_fac1 <- hz1(oe_data$age) * oe_data$e1
oe_data_aux$real_fac2 <- mu10(oe_data$age) * oe_data$e2
oe_data_aux$est_fac2 <- hz2(oe_data$age) * oe_data$e2
ad_fac <- aggregate(cbind(real_fac1, est_fac1, real_fac2, est_fac2) ~ group, data = oe_data_aux, FUN = sum)
adjustment1 <- ad_fac$est_fac1 / ad_fac$real_fac1
adjustment2 <- ad_fac$est_fac2 / ad_fac$real_fac2

plot(original_theta[original_theta$id == "gamma_ind", 3], ethetha1 * adjustment1)
abline(0, 1, col = "blue")

plot(original_theta[original_theta$id == "gamma_ind", 4], ethetha2 * adjustment2)
abline(0, 1, col = "blue")


plot(original_theta[original_theta$id == "lnorm_ind", 3], ethetha1 * adjustment1)
abline(0, 1, col = "blue")

plot(original_theta[original_theta$id == "lnorm_ind", 4], ethetha2 * adjustment2)
abline(0, 1, col = "blue")

plot(original_theta[original_theta$id == "lnorm_pos", 3], ethetha1 * adjustment1)
abline(0, 1, col = "blue")

plot(original_theta[original_theta$id == "lnorm_pos", 4], ethetha2 * adjustment2)
abline(0, 1, col = "blue")

plot(original_theta[original_theta$id == "lnorm_neg", 3], ethetha1 * adjustment1)
abline(0, 1, col = "blue")

plot(original_theta[original_theta$id == "lnorm_neg", 4], ethetha2 * adjustment2)
abline(0, 1, col = "blue")

plot(rep(original_theta[original_theta$id == "gamma_pos", 3], 2), ethetha1 * adjustment1)
abline(0, 1, col = "blue")

plot(rep(original_theta[original_theta$id == "gamma_pos", 4], 2), ethetha2 * adjustment2)
abline(0, 1, col = "blue")

plot(original_theta[original_theta$id == "gamma_neg", 3], ethetha1 * adjustment1)
abline(0, 1, col = "blue")

plot(original_theta[original_theta$id == "gamma_neg", 4], ethetha2 * adjustment2)
abline(0, 1, col = "blue")

sdlog1 <- 0.5122151
meanlog1 <- -0.1311821
z <- seq(0.001, 5, by = 0.001)
plot(z, dlnorm(z, meanlog1, sdlog1), type = "l")
lines(z, dens(marginal(bph_fit, 1), z * mean1) * mean1, col = "red")
lines(z, dens(marginal(bph_fit, 2), z * mean2) * mean2, col = "blue")
lines(z, dgamma(z, 2.438097, 2.438097), col = "blue")


cor(sim(bph_fit, 10000), method = "kendall")

lapply(coef(bph_fit), round, digits=6)


psi <- 10 / 3
z <- seq(0.1, 5, by = 0.001)
plot(z, dgamma(z, psi, psi), type = "l")
lines(z, dens(marginal(bph_fit, 1), z), col = "red")


psi <- 10 / 3
z <- seq(0.1, 5, by = 0.001)
plot(z, dgamma(z, psi, psi), type = "l")
lines(z, dens(marginal(bph_fit, 2), z), col = "red")


theta1_or <- original_theta[original_theta$id == "gamma_ind", 3]
theta2_or <- original_theta[original_theta$id == "gamma_ind", 4]

theta1_or <- original_theta[original_theta$id == "gamma_pos", 3]
theta2_or <- original_theta[original_theta$id == "gamma_pos", 4]

theta1_or <- original_theta[original_theta$id == "lnorm_ind", 3]
theta2_or <- original_theta[original_theta$id == "lnorm_ind", 4]

theta1_or <- original_theta[original_theta$id == "lnorm_neg", 3]
theta2_or <- original_theta[original_theta$id == "lnorm_neg", 4]

hz1_theta <- function(x, g) {
  ethetha1[g] * exp(beta1_fit[1] + beta1_fit[2] * x + beta1_fit[3] * x^2)
}

mu01_theta <- function(x, g) {
  theta1_or[g] * exp(-3.2 - 0.025 * x + 0.0006 * x ^ 2)
  #theta1_or[g] * exp(-4.5 - 0.018 * x + 0.00064 * x ^ 2)
}

gr <- 50

plot(x, hz1_theta(x, gr), type = "l", col = "red", ylim = c(0, 0.15))
lines(x, mu01_theta(x, gr))


hz2_theta <- function(x, g) {
  ethetha2[g] * exp(beta2_fit[1] + beta2_fit[2] * x)
}

mu10_theta <- function(x, g) {
  theta2_or[g] * exp(0.3 - 0.049 * x)
}

plot(x, hz2_theta(x, gr), type = "l", col = "red", ylim = c(0, 1))
lines(x, mu10_theta(x, gr))


