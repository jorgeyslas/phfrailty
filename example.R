oe_ln <- read.csv("oe_lnorm_20211204.csv")

head(oe_ln)

oe <- oe_ln[oe_ln$no == 12, ]

head(oe)
nrow(oe)
summary(oe)

set.seed(1)
bivph_obj <- bivphasetype(dimensions = c(3, 3))
par <- coef(bivph_obj)
par

mp_cor_dens(matrix(c(1, 2), ncol = 2), par$alpha, par$S11, par$S12, par$S22)

mp_cor_dens_cov(matrix(c(0, 0), ncol = 2), matrix(c(0, 0), ncol = 2), par$alpha, par$S11, par$S12, par$S22)

y <- matrix(c(1, 2), ncol = 2)
aux <- mp_cor_dens_aux(y, matrix(c(1, 0), ncol = 2), par$alpha, par$S11, par$S12, par$S22)

ethetha1 <- (y[, 1] + 1) * aux$dens_aux1 / aux$density
ethetha1

ethetha2 <- (y[, 2] + 1) * aux$dens_aux2 / aux$density
ethetha2


oe_data1 <- oe[oe$e1 > 0, ]
oe_data2 <- oe[oe$e2 > 0, ]

# Initial values
m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1)), family = poisson(link = "log"), data = oe_data1)
m2 <- glm(o2 ~ age + offset(log(e2)), family = poisson(link = "log"), data = oe_data2)

beta1 <- rep(0, 3)
beta2 <- rep(0, 2)

beta1 <- coef(m1)
beta2 <- coef(m2)

X1 <- cbind(rep(1, nrow(oe)), oe$age, oe$age^2)
X2 <- X1[, c(1, 2)]

factor1 <- oe$e1 * exp(X1 %*% beta1)
factor2 <- oe$e2 * exp(X2 %*% beta2)
fac <- cbind(factor1, factor2)
head(fac)

head(oe)
y <- as.matrix(oe[, c(6, 8)])

fac_aux <- fac
fac_aux[fac[, 1] == 0, 1] <- 1
fac_aux[fac[, 2] == 0, 2] <- 1

sum(log( fac_aux[, 1] * fac_aux[, 2] * mp_cor_dens_cov(y, fac, par$alpha, par$S11, par$S12, par$S22)))


aux <- mp_cor_dens_aux(y, fac, par$alpha, par$S11, par$S12, par$S22)

ethetha1 <- (as.vector(y[, 1]) + 1) * aux$dens_aux1 / aux$density
ethetha1

ethetha2 <- (as.vector(y[, 2]) + 1) * aux$dens_aux2 / aux$density
ethetha2

oe_data1$etheta <- ethetha1[oe$e1 > 0]
oe_data2$etheta <- ethetha2[oe$e2 > 0]

m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 * etheta)), family = poisson(link = "log"), data = oe_data1)
m2 <- glm(o2 ~ age + offset(log(e2 * etheta)), family = poisson(link = "log"), data = oe_data2)

beta1 <- coef(m1)
beta2 <- coef(m2)


factor1 <- oe$e1 * exp(X1 %*% beta1)
factor2 <- oe$e2 * exp(X2 %*% beta2)
fac <- cbind(factor1, factor2)

fac_aux <- fac
fac_aux[fac[, 1] == 0, 1] <- 1
fac_aux[fac[, 2] == 0, 2] <- 1

fac_ind <- (fac[, 1] > 0) | (fac[, 2] > 0)

type_ind <- matrix(rep(0, 3 * length(y[, 1])), ncol = 3)

type_ind[oe$e1 > 0 & oe$e2 > 0, 1] = 1

type_ind[oe$e1 > 0 & oe$e2 == 0, 2] = 1

type_ind[oe$e1 == 0 & oe$e2 > 0, 3] = 1

sum(type_ind)

sel_dens_t <- function(x, alpha, alpha_aux, S11, S12, S22, type) {
  if (type == 0) {
    bivph_density(x, alpha, S11, S12, S22)
  } else if(type == 1) {
    ph_density(x, alpha, S11)
  } else {
    ph_density(x, alpha_aux, S22)
  }
}
sel_dens <- Vectorize(sel_dens_t, "type")

sel_dens(1, alpha_fit, alpha_aux, S11_fit, S12_fit, S22_fit, type_ind)

conditional_density <- function(z, alphafn, alpha_aux, S11fn, S12fn, S22fn, fac_or, fac_aux, obs, dens_eval) {
  sum(fac_aux * (z[, 1])^obs[, 1] * (z[, 2])^obs[, 2] * exp(-z[, 1] * fac_or[, 1] - z[, 2] * fac_or[, 2]) * bivph_density(z, alphafn, S11fn, S12fn, S22fn) / (factorial(obs[, 1]) * factorial(obs[, 2]) * dens_eval)) / length(obs[, 1])
}

bivph_par <- bivph_obj@pars
alpha_fit <- clone_vector(bivph_par$alpha)
S11_fit <- clone_matrix(bivph_par$S11)
S12_fit <- clone_matrix(bivph_par$S12)
S22_fit <- clone_matrix(bivph_par$S22)

alpha_aux <- as.vector(alpha_fit %*% solve(-S11_fit) %*% S12_fit)

initialpoint1 = 0.0001
truncationpoint1 = 20
delta1 = 0.1
initialpoint2 = 0.0001
truncationpoint2 = 20
delta2 = 0.1

# Discretization of density
prob <- numeric(0)
value <- base::expand.grid(z1 = seq(initialpoint1, truncationpoint1, by = delta1), z2 = seq(initialpoint2, truncationpoint2, by = delta2))
for (l in 1:length(value[, 1])) {
  prob[l] <- conditional_density(matrix(c(value[l, 1], value[l, 2]), ncol = 2), alpha_fit, alpha_aux, S11_fit, S12_fit, S22_fit, fac, fac_ind, y, aux$density)
}
sum(prob * delta1 * delta2)



