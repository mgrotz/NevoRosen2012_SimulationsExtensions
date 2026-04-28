library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)


###############################################################################
###### 0. Functions                                                      ######
###############################################################################

find_admissible_gamma <- function(df, grid = seq(0.01, 0.99, length.out = 99)) {
  X  <- df$ybar
  Z1 <- df$Gx
  Z2 <- df$G2x
  u  <- df$eps
  
  res <- data.frame(gamma = grid) %>%
    mutate(
      cov_omega_u = sapply(gamma, function(gam) {
        omega <- gam * Z2 - (1 - gam) * Z1
        cov(omega, u)
      }),
      prop_condition = sapply(gamma, function(gam) {
        omega <- gam * Z2 - (1 - gam) * Z1
        (var(X) * sd(omega) - cov(omega, X) * sd(X)) * cov(omega, X)
      }),
      cond1 = cov_omega_u >= 0,
      cond2 = prop_condition < 0,
      both_hold = cond1 & cond2
    )
  
  idx <- which(res$both_hold)
  
  list(
    exists = length(idx) > 0,
    gamma  = if (length(idx) > 0) res$gamma[idx[1]] else NA_real_
  )
}


simulate_matching_model <- function(
    n = 1000,
    beta = .75,
    gamma = 4,
    kappa_eps_x = 0.01,
    kappa_x_x = 0.01,
    target_density = 0.08,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  x <- rnorm(n)
  eps <- rnorm(n)
  
  eta <- matrix(rnorm(n * n), n, n)
  
  eps_x_score <- outer(eps, x, "*") + outer(x, eps, "*")
  x_x_score   <- outer(x, x, "*")
  
  score <- eta +
    kappa_eps_x * eps_x_score +
    kappa_x_x * x_x_score
  
  diag(score) <- -Inf
  
  cutoff <- quantile(score[upper.tri(score)], 1 - target_density)
  
  A <- (score > cutoff) * 1
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  
  deg <- rowSums(A)
  G <- matrix(0, n, n)
  idx <- deg > 0
  G[idx, ] <- A[idx, ] / deg[idx]
  
  # higher-order instruments
  Gx  <- as.vector(G %*% x)
  G2x <- as.vector(G %*% G %*% x)
  G3x <- as.vector(G %*% G %*% G %*% x)
  G4x <- as.vector(G %*% G %*% G %*% G %*% x)
  G5x <- as.vector(G %*% G %*% G %*% G %*% G %*% x)
  
  y <- solve(diag(n) - beta * G, gamma * x + eps)
  
  ybar <- as.vector(G %*% y)
  
  df <- data.frame(
    y, ybar, x, eps,
    Gx, G2x, G3x, G4x, G5x,
    deg
  ) %>%
    filter(deg > 0)
  
  return(df)
}


###############################################################################
###### 1. Single simulation                                              ######
###############################################################################

df <- simulate_matching_model(seed = 1)

ols_fit <- lm(y ~ ybar + x, data = df)
iv_fit  <- ivreg(y ~ ybar + x | G2x + x, data = df)

c(
  cor_eps_ybar = cor(df$eps, df$ybar),
  
  cor_eps_Gx   = cor(df$eps, df$Gx),
  cor_eps_G2x  = cor(df$eps, df$G2x),
  
  # Condition 1: same sign
  cond1_Gx  = cor(df$eps, df$Gx)  * cor(df$eps, df$ybar) > 0,
  cond1_G2x = cor(df$eps, df$G2x) * cor(df$eps, df$ybar) > 0,
  
  # Condition from prop (have to be greater than 0)
  test1 = (var(df$ybar) * sd(df$Gx) - cov(df$ybar, df$Gx) * sd(df$ybar)) * cov(df$ybar, df$Gx),
  test1 = (var(df$ybar) * sd(df$G2x) - cov(df$ybar, df$G2x) * sd(df$ybar)) * cov(df$ybar, df$G2x)
)



###############################################################################
###### 2. Monte Carlo Simulations                                        ######
###############################################################################

#==========================#
# 2.1 Running simulations  #
#==========================#

R <- 100
beta_true <- 0.75

res <- lapply(1:R, function(r) {
  df <- simulate_matching_model(seed = r, beta = beta_true)
  
  ols_fit <- lm(y ~ ybar + x, data = df)
  iv_fit  <- ivreg(y ~ ybar + x | G2x + x, data = df)
  
  # admissible gamma
  gam_out <- find_admissible_gamma(df)
  admissible_gamma <- gam_out$exists
  gamma_star <- gam_out$gamma
  
  # basic estimates
  beta_ols <- coef(ols_fit)["ybar"]
  beta_iv_Gx  <- coef(ivreg(y ~ ybar + x | Gx + x, data = df))["ybar"]
  beta_iv_G2x <- coef(iv_fit)["ybar"]
  
  # bounds if admissible gamma exists
  lower <- NA_real_
  upper <- NA_real_
  beta_iv_omega <- NA_real_
  
  if (admissible_gamma) {
    omega <- gamma_star * df$G2x - (1 - gamma_star) * df$Gx
    
    df_tmp <- df
    df_tmp$omega <- omega
    
    beta_iv_omega <- coef(
      ivreg(y ~ ybar + x | omega + x, data = df_tmp)
    )["ybar"]
    
    lower <- beta_iv_omega
    upper <- min(beta_iv_Gx, beta_iv_G2x, beta_ols, na.rm = TRUE)
  }
  
  cor_eps_ybar <- cor(df$eps, df$ybar)
  cor_Gx_eps  <- cor(df$Gx, df$eps)
  cor_G2x_eps <- cor(df$G2x, df$eps)
  
  data.frame(
    rep = r,
    
    beta_ols = beta_ols,
    beta_iv_Gx = beta_iv_Gx,
    beta_iv_G2x = beta_iv_G2x,
    beta_iv_omega = beta_iv_omega,
    
    lower = lower,
    upper = upper,
    valid_bound = is.finite(lower) & is.finite(upper) & lower <= upper,
    
    admissible_gamma = admissible_gamma,
    gamma_star = gamma_star,
    
    cor_eps_ybar = cor_eps_ybar,
    cor_Gx_eps = cor_Gx_eps,
    cor_G2x_eps = cor_G2x_eps,
    
    cond1_Gx  = cor_Gx_eps  * cor_eps_ybar > 0,
    cond1_G2x = cor_G2x_eps * cor_eps_ybar > 0,
    
    condP_Gx  = (var(df$ybar) * sd(df$Gx)  -
                   cov(df$ybar, df$Gx)  * sd(df$ybar)) *
      cov(df$ybar, df$Gx) > 0,
    
    condP_G2x = (var(df$ybar) * sd(df$G2x) -
                   cov(df$ybar, df$G2x) * sd(df$ybar)) *
      cov(df$ybar, df$G2x) > 0
  )
}) %>%
  bind_rows()

summary(res)

#==========================#
# 2.2 Plotting functions   #
#==========================#

# 2.2.1 OLS and IV estimates
est_df <- res %>%
  select(beta_ols, beta_iv_G2x) %>%
  pivot_longer(cols = everything(),
               names_to = "estimator",
               values_to = "estimate") %>%
  mutate(
    estimator = case_when(
      estimator == "beta_ols" ~ "OLS",
      estimator == "beta_iv_G2x" ~ "IV"
    )
  )

# means
means_df <- est_df %>%
  group_by(estimator) %>%
  summarise(mean_est = mean(estimate, na.rm = TRUE), .groups = "drop")

bounds_df <- res %>%
  filter(valid_bound) %>%
  mutate(
    y_pos = -0.02 + row_number() * 0.2   # shift each bracket slightly upward
  )



p <- ggplot(est_df, aes(x = estimate, fill = estimator)) +
  
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, alpha = 0.4, position = "identity") +
  
  # mean lines
  geom_vline(data = means_df,
             aes(xintercept = mean_est, color = estimator),
             linetype = "dashed", linewidth = 1) +
  
  # true beta
  geom_vline(xintercept = beta_true,
             color = "red", linetype = "dotted", linewidth = 1.2) +
  
  # bounds
  geom_segment(data = bounds_df,
               aes(x = lower, xend = upper,
                   y = y_pos, yend = y_pos),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.7) +
  
  geom_segment(data = bounds_df,
               aes(x = lower, xend = lower,
                   y = y_pos - 0.005, yend = y_pos + 0.005),
               inherit.aes = FALSE,
               color = "black") +
  
  geom_segment(data = bounds_df,
               aes(x = upper, xend = upper,
                   y = y_pos - 0.005, yend = y_pos + 0.005),
               inherit.aes = FALSE,
               color = "black") +
  
  scale_fill_manual(values = c("OLS" = "steelblue", "IV" = "orange")) +
  scale_color_manual(values = c("OLS" = "steelblue", "IV" = "orange")) +
  
  labs(
    x = expression(hat(beta)),
    y = "Density"
  ) +
  theme_minimal()

# export

ggsave(
  filename = "plots/ols_iv_bounds.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)


# 2.2.2 Histogram of correlations

plot_df <- res %>%
  select(cor_eps_ybar, cor_Gx_eps, cor_G2x_eps) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = c("cor_eps_ybar", "cor_Gx_eps", "cor_G2x_eps"),
      labels = c("Corr(ybar, ε)", "Corr(Gx, ε)", "Corr(G²x, ε)")
    )
  )

p <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(
    x = "",
    y = "Correlation"
  ) +
  theme_minimal()

ggsave(
  filename = "plots/corr_networks.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)


###############################################################################
###### 3. Monte Carlo Simulations Grid Search over Parameters            ######
###############################################################################

# parameter grid
param_grid <- expand.grid(
  beta = c(0.5, 0.75, 0.95),
  gamma = c(1, 4, 10),
  kappa_eps_x = c(0, 0.01, 0.1, 0.5),
  kappa_x_x = c(0, 0.01, 0.1, 0.5),
  target_density = c(0.02, 0.08, 0.2)
)

R <- 40

all_results <- list()

counter <- 1
total <- nrow(param_grid)

for (i in 1:nrow(param_grid)) {
  
  params <- param_grid[i, ]
  
  cat(sprintf("Running %d / %d\n", i, total))
  
  res <- lapply(1:R, function(r) {
    
    df <- simulate_matching_model(
      seed = r,
      beta = params$beta,
      gamma = params$gamma,
      kappa_eps_x = params$kappa_eps_x,
      kappa_x_x = params$kappa_x_x,
      target_density = params$target_density
    )
    
    ols_fit <- lm(y ~ ybar + x, data = df)
    iv_fit  <- ivreg(y ~ ybar + x | G2x + x, data = df)
    
    gam_out <- find_admissible_gamma(df)
    
    beta_ols <- coef(ols_fit)["ybar"]
    beta_iv_Gx  <- coef(ivreg(y ~ ybar + x | Gx + x, data = df))["ybar"]
    beta_iv_G2x <- coef(iv_fit)["ybar"]
    
    lower <- NA_real_
    upper <- NA_real_
    beta_iv_omega <- NA_real_
    
    if (gam_out$exists) {
      omega <- gam_out$gamma * df$G2x - (1 - gam_out$gamma) * df$Gx
      df_tmp <- df
      df_tmp$omega <- omega
      
      beta_iv_omega <- coef(
        ivreg(y ~ ybar + x | omega + x, data = df_tmp)
      )["ybar"]
      
      lower <- beta_iv_omega
      upper <- min(beta_iv_Gx, beta_iv_G2x, beta_ols, na.rm = TRUE)
    }
    
    cor_eps_ybar <- cor(df$eps, df$ybar)
    cor_Gx_eps  <- cor(df$Gx, df$eps)
    cor_G2x_eps <- cor(df$G2x, df$eps)
    
    data.frame(
      beta_ols = beta_ols,
      beta_iv_Gx = beta_iv_Gx,
      beta_iv_G2x = beta_iv_G2x,
      lower = lower,
      upper = upper,
      admissible_gamma = gam_out$exists,
      cor_eps_ybar = cor_eps_ybar,
      cor_Gx_eps = cor_Gx_eps,
      cor_G2x_eps = cor_G2x_eps
    )
    
  }) %>% bind_rows()
  
  # flatten summary(res)
  summ <- summary(res) %>%
    as.data.frame() %>%
    rownames_to_column("stat") %>%
    pivot_longer(-stat, names_to = "variable", values_to = "value")
  
  # attach parameters
  summ <- summ %>%
    mutate(
      beta = params$beta,
      gamma = params$gamma,
      kappa_eps_x = params$kappa_eps_x,
      kappa_x_x = params$kappa_x_x,
      target_density = params$target_density
    )
  
  all_results[[counter]] <- summ
  counter <- counter + 1
}

final_df <- bind_rows(all_results)


