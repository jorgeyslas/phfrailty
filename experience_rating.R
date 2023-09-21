
experience_rating_bph <- function(x,
                                  oe_data,
                                  stepsEM = 5,
                                  stepsPH = 20,
                                  initialpoint1 = 0.0001,
                                  truncationpoint1 = 8,
                                  delta1 = 0.015,
                                  initialpoint2 = 0.0001,
                                  truncationpoint2 = 8,
                                  delta2 = 0.015,
                                  every = 1) {

  start_time <- Sys.time()
  # Biv PH parameters
  x_par <- x@pars
  alpha_fit <- clone_vector(x_par$alpha)
  S11_fit <- clone_matrix(x_par$S11)
  S12_fit <- clone_matrix(x_par$S12)
  S22_fit <- clone_matrix(x_par$S22)

  # Initial values of beta
  m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 + 1e-10)), family = poisson(link = "log"), data = oe_data)
  m2 <- glm(o2 ~ age + offset(log(e2 + 1e-10)), family = poisson(link = "log"), data = oe_data)

  beta1_fit <- coef(m1)
  beta2_fit <- coef(m2)

  oe_data$o1_hat <- fitted(m1)
  oe_data$o2_hat <- fitted(m2)

  G <- max(oe_data$group)
  aggreated_data <- NULL
  aggreated_data <- aggregate(cbind(o1, o1_hat, o2, o2_hat) ~ group + no + id, data = oe_data, FUN = sum)

  y <- as.matrix(aggreated_data[, c(4, 6)])
  fac <- as.matrix(aggreated_data[, c(5, 7)])

  for (k in 1:stepsEM) {
    aux <- mp_cor_dens_aux(y, fac, alpha_fit, S11_fit, S12_fit, S22_fit)

    ethetha1 <- (as.vector(y[, 1]) + 1) * (aux$dens_aux1) / (aux$density * fac[, 1])
    ethetha2 <- (as.vector(y[, 2]) + 1) * (aux$dens_aux2) / (aux$density * fac[, 2])

    oe_data$theta1 <- rep(ethetha1, sapply(1:G, function(g) {
      length(oe_data$age[oe_data$group == g])
    }))

    oe_data$theta2 <- rep(ethetha2, sapply(1:G, function(g) {
      length(oe_data$age[oe_data$group == g])
    }))

    m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 * theta1 + 1e-10)), family = poisson(link = "log"), data = oe_data)
    m2 <- glm(o2 ~ age + offset(log(e2 * theta2 + 1e-10)), family = poisson(link = "log"), data = oe_data)

    beta1_fit <- coef(m1)
    beta2_fit <- coef(m2)

    conditional_density <- function(z, alphafn, S11fn, S12fn, S22fn, fac_or, obs, dens_eval) {
      sum(exp(obs[, 1] * log(z[, 1] * fac_or[, 1]) + obs[, 2] * log(z[, 2] * fac_or[, 2]) - z[, 1] * fac_or[, 1] - z[, 2] * fac_or[, 2] - lgamma(obs[, 1] + 1) - lgamma(obs[, 2] + 1)) * bivph_density(z, alphafn, S11fn, S12fn, S22fn) / (dens_eval)) / length(obs[, 1])
    }

    # Discretization of density
    prob <- numeric(0)
    value <- base::expand.grid(z1 = seq(initialpoint1, truncationpoint1, by = delta1), z2 = seq(initialpoint2, truncationpoint2, by = delta2))
    for (l in 1:length(value[, 1])) {
      prob[l] <- conditional_density(matrix(c(value[l, 1], value[l, 2]), ncol = 2), alpha_fit, S11_fit, S12_fit, S22_fit, fac, y, aux$density / (fac[, 1] * fac[, 2]))
    }

    # PH fitting
    for (l in 1:stepsPH) {
      EMstep_bivph(alpha_fit, S11_fit, S12_fit, S22_fit, matrix(c(value[, 1], value[, 2]), ncol = 2), prob)
    }

    oe_data$o1_hat <- exp(beta1_fit[1] + beta1_fit[2] * oe_data$age + beta1_fit[3] * oe_data$age^2) * (oe_data$e1 + 1e-15)
    oe_data$o2_hat <- exp(beta2_fit[1] + beta2_fit[2] * oe_data$age) * (oe_data$e2 + 1e-15)

    aggreated_data <- aggregate(cbind(o1, o1_hat, o2, o2_hat) ~ group + no + id, data = oe_data, FUN = sum)

    fac <- as.matrix(aggreated_data[, c(5, 7)])

    if (k %% every == 0) {
      s1 <- sum(lgamma(aggreated_data[, 4] + 1)) + sum(lgamma(aggreated_data[, 6] + 1)) - sum(lgamma(oe_data[, 6] + 1)) - sum(lgamma(oe_data[, 8] + 1))
      s2 <- sum(oe_data[, 6]* log(oe_data[, 9])) + sum(oe_data[, 8] * log(oe_data[, 10]))
      s3 <- sum((aggreated_data[, 4]) * log(aggreated_data[, 5])) + sum((aggreated_data[, 6]) * log(aggreated_data[, 7]))
      s4 <- sum(log(mp_cor_dens_cov(as.matrix(aggreated_data[, c(4, 6)]), as.matrix(aggreated_data[, c(5, 7)]), alpha_fit, S11_fit, S12_fit, S22_fit) / (aggreated_data[, 5] * aggreated_data[, 7])))
      cat("\r", "iteration:", k,
          ", logLik:", s1 + s2 - s3 + s4,
          ", prob disc:", sum(prob * delta1 * delta2),
          sep = " "
      )
    }
  }

  aux <- mp_cor_dens_aux(y, fac, alpha_fit, S11_fit, S12_fit, S22_fit)

  ethetha1 <- (as.vector(y[, 1]) + 1) * (aux$dens_aux1) / (aux$density * fac[, 1])
  ethetha2 <- (as.vector(y[, 2]) + 1) * (aux$dens_aux2) / (aux$density * fac[, 2])

  x_fit <- bivphasetype(alpha = alpha_fit, S11 = S11_fit, S12 = S12_fit, S22 = S22_fit)

  oe_out <- aggreated_data[, 1:7]

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

  out_list <- list(bph_fit = x_fit, oe_out = oe_out)

  end_time <- Sys.time()
  cat("\n Running time: ", end_time - start_time)

  return(out_list)
}


