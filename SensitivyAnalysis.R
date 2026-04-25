library(ggplot2)
library(tidyverse)
library(dplyr)


###############################################################################
##### UNIVARIATE BOUNDS                                                   ##### 
###############################################################################
compute_uni_bounds <- function(varx = 1,
                           varu = 1,
                           varz = 1,
                           beta = 0.5,
                           covuz = 0.1,
                           covux = 0.4,
                           covxz = -0.3,
                           n_grid = 200) {
  
  # --- Assumption 3: enforce same sign ---
  if (covux * covuz < 0) {
    covux <- -covux
    covxz <- -covxz
    message("Redefined X as -X")
  }
  
  # --- implied covariances ---
  covyz <- beta * covxz + covuz
  covxy <- beta * varx + covux
  
  # --- Assumption 4: Z less endogenous than X ---
  rhoux <- covux / sqrt(varu * varx)
  rhouz <- covuz / sqrt(varu * varz)
  
  if (abs(rhoux) < abs(rhouz)) {
    stop("Z more endogenous than X — stopping execution")
  }
  
  # --- ensure two-sided case ---
  if (covxz >= 0) {
    stop("covxz >= 0 → one-sided bounds (not handled here)")
  }
  
  # --- compute bounds ---
  sx <- sqrt(varx)
  sz <- sqrt(varz)
  
  delta_grid <- seq(0, 1, length.out = n_grid)
  
  beta_iv <- covyz / covxz
  
  bound_df <- data.frame(delta = delta_grid) %>%
    mutate(
      beta_v_delta =
        (sx * covyz - delta * sz * covxy) /
        (sx * covxz - delta * sz * varx),
      beta_iv = beta_iv,
      lower = pmin(beta_iv, beta_v_delta),
      upper = pmax(beta_iv, beta_v_delta)
    )
  
  # clean numerical issues
  bound_df <- subset(bound_df, is.finite(lower) & is.finite(upper))
  
  return(bound_df)
}

# Simple plot
bound_df <- compute_uni_bounds()
ggplot(bound_df, aes(x = delta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(aes(y = lower), linewidth = 0.8) +
  geom_line(aes(y = upper), linewidth = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
  labs(
    x = expression(delta),
    y = expression(B(delta)),
    title = "Identified interval as a function of delta"
  ) +
  theme_minimal()


# Plot with 4 elements

scenarios <- tibble::tibble(
  covuz = c(0.1, -0.1, 0.3, -0.3),
  covux = c(0.4, -0.4, 0.4, -0.4)
) %>%
  mutate(
    label = paste0(
      "sigma[ZU] == ", covuz,
      "*','~~sigma[XU] == ", covux
    )
  )

plot_df <- scenarios %>%
  group_by(label, covuz, covux) %>%
  group_modify(~ {
    compute_bounds(
      beta = 0.5,
      covxz = -0.3,
      covuz = .y$covuz,
      covux = .y$covux
    )
  }) %>%
  ungroup()

p <- ggplot(plot_df, aes(x = delta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(aes(y = lower), linewidth = 0.7) +
  geom_line(aes(y = upper), linewidth = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  facet_wrap(~ label, ncol = 2, labeller = label_parsed) +
  labs(
    x = expression(delta),
    y = expression(atop(M(delta), "or equivalently" ~ beta(delta)))
  ) +
  theme_minimal()
p

ggsave("plots/bounds_plot.png", plot = p, width = 8, height = 6, dpi = 600)




###############################################################################
##### BIVARIATE BOUNDS                                                   ##### 
###############################################################################

compute_multi_bounds <- function(varx = 1,
                                 varu = 1,
                                 varz1 = 1,
                                 varz2 = 1,
                                 beta = 0.5,
                                 covux = 0.4,
                                 covuz1 = 0.1,
                                 covuz2 = 0.2,
                                 covxz1 = -0.3,
                                 covxz2 = -0.5,
                                 n_grid = 200) {
  
  # --- check Assumption 3 for both instruments ---
  if (covux * covuz1 < 0 || covux * covuz2 < 0) {
    stop("Assumption 3 fails for at least one instrument.")
  }
  
  # --- correlations with U ---
  rhoux  <- covux  / sqrt(varu * varx)
  rhouz1 <- covuz1 / sqrt(varu * varz1)
  rhouz2 <- covuz2 / sqrt(varu * varz2)
  
  # --- check original Nevo-Rosen Assumption 4 ---
  if (abs(rhouz1) > abs(rhoux)) {
    stop("Z1 more endogenous than X.")
  }
  
  if (abs(rhouz2) > abs(rhoux)) {
    stop("Z2 more endogenous than X.")
  }
  
  # --- ensure two-sided case for both instruments ---
  if (covxz1 >= 0 & covxz2 >= 0) {
    stop("Both cov(X,Zj) >= 0. One-sided bounds not handled here.")
  }
  
  # --- implied covariances with Y ---
  covxy  <- beta * varx + covux
  covz1y <- beta * covxz1 + covuz1
  covz2y <- beta * covxz2 + covuz2
  
  sx  <- sqrt(varx)
  sz1 <- sqrt(varz1)
  sz2 <- sqrt(varz2)
  
  # --- delta grid ---
  delta_grid <- seq(0, 1, length.out = n_grid)
  
  grid_df <- expand.grid(
    delta1 = delta_grid,
    delta2 = delta_grid
  )
  
  # --- single-IV endpoint functions ---
  beta_iv_1 <- covz1y / covxz1
  beta_iv_2 <- covz2y / covxz2
  
  beta_v_delta_1 <- function(delta1) {
    (sx * covz1y - delta1 * sz1 * covxy) /
      (sx * covxz1 - delta1 * sz1 * varx)
  }
  
  beta_v_delta_2 <- function(delta2) {
    (sx * covz2y - delta2 * sz2 * covxy) /
      (sx * covxz2 - delta2 * sz2 * varx)
  }
  
  # --- compute bounds and intersect ---
  bound_df <- grid_df |>
    dplyr::mutate(
      beta_iv_1 = beta_iv_1,
      beta_iv_2 = beta_iv_2,
      
      beta_v_delta_1 = beta_v_delta_1(delta1),
      beta_v_delta_2 = beta_v_delta_2(delta2),
      
      lower_1 = pmin(beta_iv_1, beta_v_delta_1),
      upper_1 = pmax(beta_iv_1, beta_v_delta_1),
      
      lower_2 = pmin(beta_iv_2, beta_v_delta_2),
      upper_2 = pmax(beta_iv_2, beta_v_delta_2),
      
      lower = pmax(lower_1, lower_2),
      upper = pmin(upper_1, upper_2),
      
      empty = lower > upper,
      
      lambda1_true = rhouz1 / rhoux,
      lambda2_true = rhouz2 / rhoux
    )
  
  bound_df <- subset(
    bound_df,
    is.finite(lower) & is.finite(upper)
  )
  return(bound_df)
}

bound_df2 <- compute_multi_bounds()

beta0 <- 0.4  # or 0

plot_df <- bound_df2 %>%
  mutate(
    in_z1 = (lower_1 <= beta0 & upper_1 >= beta0),
    in_z2 = (lower_2 <= beta0 & upper_2 >= beta0),
    region = case_when(
      in_z1 & in_z2  ~ "Both",
      in_z1 & !in_z2 ~ "Only Z1",
      !in_z1 & in_z2 ~ "Only Z2",
      TRUE           ~ "Neither"
    )
  )


# Single plot
p <- ggplot(plot_df, aes(x = delta1, y = delta2)) +
  geom_tile(aes(fill = region)) +
  scale_fill_manual(
    values = c(
      "Both" = "steelblue",
      "Only Z1" = "lightblue",
      "Only Z2" = "orange",
      "Neither" = "grey90"
    ),
    name = NULL
  ) +
  labs(
    x = expression(delta[1]),
    y = expression(delta[2]),
  ) +
  theme_minimal()

p

ggsave("plots/bounds_2Z_plot.png", plot = p, width = 8, height = 6, dpi = 600)


