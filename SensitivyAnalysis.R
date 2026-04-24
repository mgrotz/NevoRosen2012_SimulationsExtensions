library(ggplot2)
library(tidyverse)
library(dplyr)

compute_bounds <- function(varx = 1,
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
  rhoux <- covux / (varu * varx)
  rhouz <- covuz / (varu * varz)
  
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
bound_df <- compute_bounds()
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
    y = expression(B(delta))
  ) +
  theme_minimal()

ggsave("plots/bounds_plot.png", plot = p, width = 8, height = 6, dpi = 600)



