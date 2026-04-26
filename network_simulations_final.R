library(AER)
library(dplyr)
library(ggplot2)

set.seed(123)

simulate_network_df <- function(
    n = 300,
    beta_true = 0.4,
    gamma = 1,
    rho = 1,
    kappa_a = 2,
    kappa_x = 2,
    target_density = 0.08,
    seps = 0,
    seta = 0,
    seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # --- latent variables and shocks ---
  a <- rnorm(n)
  x <- rnorm(n)
  eps <- rnorm(n) * seps
  u <- rho * a + eps
  
  # --- network formation ---
  eta <- matrix(rnorm(n * n), n, n)*seta
  
  score <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        score[i, j] <- eta[i, j] + kappa_a / (abs(a[i] - a[j]) + 0.1) + kappa_x / (abs(x[i] - x[j]) + 0.1)
      }
    }
  }
  
  # threshold to match density
  c_cut <- quantile(score[upper.tri(score)], 1 - target_density)
  
  A <- (score > c_cut) * 1
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  
  # --- row-normalized G ---
  deg <- rowSums(A)
  G <- matrix(0, n, n)
  non_isolated <- deg > 0
  G[non_isolated, ] <- A[non_isolated, ] / deg[non_isolated]
  
  # --- simulate outcome ---
  I_n <- diag(n)
  y <- solve(I_n - beta_true * G, gamma * x + u)
  
  # --- construct regressors/instruments ---
  Gy  <- as.vector(G %*% y)
  Gx  <- as.vector(G %*% x)
  G2x <- as.vector(G %*% G %*% x)
  
  df <- data.frame(
    y = y,
    Gy = Gy,
    x = x,
    Gx = Gx,
    G2x = G2x,
    u = u
  )
  
  # remove problematic observations
  df <- subset(df, is.finite(y) & is.finite(Gy) & is.finite(G2x))
  
  return(df)
}

R <- 100
beta_true <- 0.4

sim_results <- lapply(1:R, function(r) {
  df <- simulate_network_df(seed = r, beta_true = beta_true)
  
  ols_fit <- lm(y ~ Gy + x, data = df)
  iv_fit  <- ivreg(y ~ Gy + x | G2x + x, data = df)
  
  data.frame(
    rep = r,
    beta_ols = coef(ols_fit)["Gy"],
    beta_iv  = coef(iv_fit)["Gy"],
    cov_zx = cov(df$G2x, df$Gy),
    cov_zu = cov(df$G2x, df$u),
    cov_xu = cov(df$Gy, df$u)
  )
}) %>%
  bind_rows()

# -----------------------------
# 1. OLS and IV estimator histograms
# -----------------------------

est_df <- sim_results %>%
  select(beta_ols, beta_iv) %>%
  pivot_longer(
    cols = everything(),
    names_to = "estimator",
    values_to = "estimate"
  ) %>%
  mutate(
    estimator = case_when(
      estimator == "beta_ols" ~ "OLS",
      estimator == "beta_iv"  ~ "IV"
    )
  ) %>%
  filter((estimate < 2) & (estimate > -2))

means_df <- est_df %>%
  group_by(estimator) %>%
  summarise(mean_est = mean(estimate, na.rm = TRUE), .groups = "drop")

ggplot(est_df, aes(x = estimate, fill = estimator)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, alpha = 0.4, position = "identity") +
  geom_vline(xintercept = beta_true,
             color = "red", linetype = "dashed") +
  geom_vline(data = means_df %>% filter(estimator == "OLS"),
             aes(xintercept = mean_est),
             color = "steelblue", linetype = "dashed") +
  geom_vline(data = means_df %>% filter(estimator == "IV"),
             aes(xintercept = mean_est),
             color = "orange", linetype = "dashed") +
  scale_fill_manual(values = c("OLS" = "steelblue", "IV" = "orange")) +
  labs(
    title = "OLS vs IV estimator distributions",
    x = expression(hat(beta)),
    y = "Density"
  ) +
  theme_minimal()

# -----------------------------
# 2. Covariance histograms
# -----------------------------

cov_df <- sim_results %>%
  select(cov_zx, cov_zu, cov_xu) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "covariance",
    values_to = "value"
  )

frac_pos <- cov_df %>%
  group_by(covariance) %>%
  summarise(
    frac_positive = mean(value > 0),
    x_pos = max(value, na.rm = TRUE),
    y_pos = Inf,
    .groups = "drop"
  )

ggplot(cov_df, aes(x = value)) +
  geom_histogram(bins = 25, alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_text(
    data = frac_pos,
    aes(
      x = x_pos,
      y = y_pos,
      label = paste0("Pr(>0) = ", round(frac_positive, 2))
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.5
  ) +
  facet_wrap(~ covariance, scales = "free") +
  labs(
    title = "Distribution of covariance terms",
    x = "Covariance",
    y = "Count"
  ) +
  theme_minimal()




###############################################################################
##### Endogeneity in X, but exogenous G                                   ##### 
###############################################################################

library(AER)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 1. Simulation function
# ============================================================

simulate_exogG_xcorr_df <- function(
    n = 1000,
    p_edge = 0.08,
    beta_true = 0.4,
    gamma = 0,          # keep 0 if x only serves as instrument source
    tau = 0.7,          # x-a correlation strength
    rho = 1,            # u-a correlation strength
    seps = 1,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # latent variable, observed covariate, error
  a <- rnorm(n)
  v <- rnorm(n)
  x <- tau * a + v
  
  eps <- rnorm(n) * seps
  u <- rho * a + eps
  
  # exogenous undirected network
  A <- matrix(rbinom(n * n, 1, p_edge), n, n)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  
  # row-normalize
  deg <- rowSums(A)
  G <- matrix(0, n, n)
  keep <- deg > 0
  G[keep, ] <- A[keep, ] / deg[keep]
  
  # simulate y = beta G y + gamma x + u
  y <- solve(diag(n) - beta_true * G, gamma * x + u)
  
  # endogenous regressor and instrument
  Gy  <- as.vector(G %*% y)
  Gx  <- as.vector(G %*% x)
  G2x <- as.vector(G %*% G %*% x)
  
  df <- data.frame(y = y, Gy = Gy, x = x, Gx = Gx, G2x = G2x, u = u)
  
  df <- df %>%
    filter(is.finite(y), is.finite(Gy), is.finite(G2x), deg > 0)
  
  return(df)
}


# ============================================================
# 2. Nevo-Rosen bound helper
#    Simple model: y = beta X + u, instrument Z
# ============================================================

compute_nr_bound <- function(y, X, Z, tol = 1e-8) {
  cov_xy <- cov(X, y)
  cov_zy <- cov(Z, y)
  cov_xz <- cov(X, Z)
  cov_xu <- cov(X, y - 0) # placeholder not used here
  sx <- sd(X)
  sz <- sd(Z)
  
  beta_ols <- cov_xy / var(X)
  beta_iv  <- if (abs(cov_xz) < tol) NA_real_ else cov_zy / cov_xz
  
  V1 <- sx * Z - sz * X
  beta_v1 <- if (abs(cov(V1, X)) < tol) {
    NA_real_
  } else {
    cov(V1, y) / cov(V1, X)
  }
  
  lower <- NA_real_
  upper <- NA_real_
  
  # Two-sided Nevo-Rosen bound under cov(X,Z) < 0
  if (is.finite(beta_iv) && is.finite(beta_v1) && cov_xz < -tol) {
    lower <- min(beta_iv, beta_v1)
    upper <- max(beta_iv, beta_v1)
  }
  
  data.frame(
    beta_ols = beta_ols,
    beta_iv = beta_iv,
    beta_v1 = beta_v1,
    nr_lower = lower,
    nr_upper = upper
  )
}


# ============================================================
# 3. Monte Carlo loop
# ============================================================

R <- 100
beta_true <- 0.4

sim_results <- lapply(1:R, function(r) {
  
  df <- simulate_exogG_xcorr_df(
    seed = r,
    beta_true = beta_true,
    tau = 0.7,
    rho = 1,
    gamma = 0
  )
  
  y <- df$y
  X <- df$Gy
  Z <- df$G2x
  u <- df$u
  
  # OLS and IV
  ols_fit <- lm(y ~ X)
  iv_fit  <- ivreg(y ~ X | Z)
  
  # Nevo-Rosen bound
  nr <- compute_nr_bound(y, X, Z)
  
  data.frame(
    rep = r,
    beta_ols = coef(ols_fit)["X"],
    beta_iv = coef(iv_fit)["X"],
    beta_v1 = nr$beta_v1,
    nr_lower = nr$nr_lower,
    nr_upper = nr$nr_upper,
    nr_contains_true = nr$nr_lower <= beta_true & nr$nr_upper >= beta_true,
    
    cov_zx = cov(Z, X),
    cov_zu = cov(Z, u),
    cov_xu = cov(X, u),
    
    same_sign = cov(Z, u) * cov(X, u) >= 0,
    less_endog = abs(cov(Z, u)) <= abs(cov(X, u)),
    two_sided = cov(Z, X) < 0
  )
}) %>%
  bind_rows()


# ============================================================
# 4. Estimator distributions
# ============================================================

est_df <- sim_results %>%
  select(beta_ols, beta_iv) %>%
  pivot_longer(
    cols = everything(),
    names_to = "estimator",
    values_to = "estimate"
  ) %>%
  mutate(
    estimator = case_when(
      estimator == "beta_ols" ~ "OLS",
      estimator == "beta_iv" ~ "IV"
    )
  ) # %>% filter(is.finite(estimate), estimate > -10, estimate < 10)

means_df <- est_df %>%
  group_by(estimator) %>%
  summarise(mean_est = mean(estimate), .groups = "drop")

ggplot(est_df, aes(x = estimate, fill = estimator)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, alpha = 0.4, position = "identity") +
  geom_vline(xintercept = beta_true,
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(data = means_df,
             aes(xintercept = mean_est, color = estimator),
             linewidth = 1.1) +
  scale_fill_manual(values = c("OLS" = "steelblue", "IV" = "orange")) +
  scale_color_manual(values = c("OLS" = "steelblue", "IV" = "orange")) +
  labs(
    title = "Distribution of OLS and IV estimators",
    x = expression(hat(beta)),
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme_minimal()


# ============================================================
# 5. Covariance distributions
# ============================================================

cov_df <- sim_results %>%
  select(cov_zx, cov_zu, cov_xu) %>%
  pivot_longer(
    cols = everything(),
    names_to = "covariance",
    values_to = "value"
  ) %>%
  mutate(
    covariance = case_when(
      covariance == "cov_zx" ~ "Cov(Z, X)",
      covariance == "cov_zu" ~ "Cov(Z, U)",
      covariance == "cov_xu" ~ "Cov(X, U)"
    )
  )

frac_pos <- cov_df %>%
  group_by(covariance) %>%
  summarise(
    frac_positive = mean(value > 0, na.rm = TRUE),
    x_pos = max(value, na.rm = TRUE),
    .groups = "drop"
  )

print(frac_pos)

ggplot(cov_df, aes(x = value)) +
  geom_histogram(bins = 25, alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_text(
    data = frac_pos,
    aes(
      x = x_pos,
      y = Inf,
      label = paste0("Fraction > 0 = ", round(frac_positive, 2))
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.5
  ) +
  facet_wrap(~ covariance, scales = "free") +
  labs(
    title = "Distribution of covariance terms",
    x = "Covariance",
    y = "Count"
  ) +
  theme_minimal()

# compute share
share_pos <- mean(sim_results$cov_zx * sim_results$cov_xu > 0, na.rm = TRUE)

ggplot(sim_results, aes(x = cov_xu, y = cov_zx,
                        color = (cov_zx * cov_xu > 0))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("grey70", "blue"),
                     labels = c("Opposite sign", "Same sign"),
                     name = NULL) +
  annotate(
    "text",
    x = max(sim_results$cov_xu, na.rm = TRUE),
    y = max(sim_results$cov_zx, na.rm = TRUE),
    label = paste0("Share same sign = ", round(share_pos, 2)),
    hjust = 1.1,
    vjust = 1.5
  ) +
  labs(
    x = expression(Cov(X, U)),
    y = expression(Cov(Z, X)),
    title = "Scatter: Nevo–Rosen geometry"
  ) +
  theme_minimal()


# ============================================================
# 6. Summary of Nevo-Rosen conditions
# ============================================================

sim_results %>%
  summarise(
    frac_same_sign = mean(same_sign, na.rm = TRUE),
    frac_less_endog = mean(less_endog, na.rm = TRUE),
    frac_two_sided = mean(two_sided, na.rm = TRUE),
    frac_nr_contains_true = mean(nr_contains_true, na.rm = TRUE)
  )


###############################################################################
##### Peer average model                                ##### 
###############################################################################

simulate_peer_averages <- function(
    n = 300,
    p_edge = 0.04,
    beta = 0.4,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # --- undirected random graph ---
  A <- matrix(rbinom(n * n, 1, p_edge), n, n)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  
  # --- row-normalize adjacency matrix ---
  deg <- rowSums(A)
  G <- matrix(0, n, n)
  keep <- deg > 0
  G[keep, ] <- A[keep, ] / deg[keep]
  
  # --- errors ---
  eps <- rnorm(n)
  
  # --- solve y = beta G y + eps ---
  y <- solve(diag(n) - beta * G, eps)
  
  # --- helper: exact kth-neighbor matrix ---
  # kth neighbors = reachable in exactly k steps,
  # excluding self and excluding all closer neighbors
  exact_k_neighbors <- function(A, k) {
    Ak <- (A %^% k) > 0
    
    closer <- matrix(FALSE, nrow(A), ncol(A))
    if (k > 1) {
      for (m in 1:(k - 1)) {
        closer <- closer | ((A %^% m) > 0)
      }
    }
    
    M <- Ak & !closer
    diag(M) <- FALSE
    M * 1
  }
  
  # matrix power helper
  `%^^%` <- function(M, k) {
    if (k == 1) return(M)
    out <- M
    for (i in 2:k) out <- out %*% M
    out
  }
  
  `%^^%` <- `%^^%`
  
  `%^%` <- function(M, k) {
    if (k == 1) return(M)
    out <- M
    for (i in 2:k) out <- out %*% M
    out
  }
  
  # row-normalize any neighbor matrix
  row_normalize <- function(M) {
    d <- rowSums(M)
    W <- matrix(0, nrow(M), ncol(M))
    idx <- d > 0
    W[idx, ] <- M[idx, ] / d[idx]
    W
  }
  
  # --- neighbor layers ---
  A1 <- A
  A2 <- exact_k_neighbors(A, 2)
  A3 <- exact_k_neighbors(A, 3)
  A4 <- exact_k_neighbors(A, 4)
  
  G1 <- row_normalize(A1)
  G2 <- row_normalize(A2)
  G3 <- row_normalize(A3)
  G4 <- row_normalize(A4)
  
  # --- averages of y over kth neighbors ---
  ybar_1 <- as.vector(G1 %*% y)
  ybar_2 <- as.vector(G2 %*% y)
  ybar_3 <- as.vector(G3 %*% y)
  ybar_4 <- as.vector(G4 %*% y)
  
  # --- empirical covariances with error ---
  covariances <- data.frame(
    layer = c("neighbors", "second_neighbors", "third_neighbors", "fourth_neighbors"),
    covariance_with_eps = c(
      cov(eps, ybar_1),
      cov(eps, ybar_2),
      cov(eps, ybar_3),
      cov(eps, ybar_4)
    ),
    n_nonempty = c(
      sum(rowSums(A1) > 0),
      sum(rowSums(A2) > 0),
      sum(rowSums(A3) > 0),
      sum(rowSums(A4) > 0)
    )
  )
  
  list(
    data = data.frame(
      id = 1:n,
      y = y,
      eps = eps,
      ybar_1 = ybar_1,
      ybar_2 = ybar_2,
      ybar_3 = ybar_3,
      ybar_4 = ybar_4
    ),
    A = A,
    G = G,
    covariances = covariances
  )
}

sim <- simulate_peer_averages(seed = 123)

sim$covariances
head(sim$data)



###############################################################################
########### Stnadard model ########################
###############################################################################

simulate_network <- function(
    n = 1000,
    beta = 0.75,
    gamma = 1,
    kappa_a = 2,
    seps = 1,
    target_density = 0.02,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # unobserved scalar type
  a <- rnorm(n)
  
  # observed exogenous characteristic
  x <- rnorm(n)
  
  # structural shock independent of a and x
  eps <- rnorm(n)*seps
  
  # network formation: close a_i and a_j more likely to link
  eta <- matrix(rnorm(n * n), n, n)
  
  dist_a <- abs(outer(a, a, "-"))
  score <- eta + kappa_a / (dist_a + 0.1)
  diag(score) <- -Inf
  
  cutoff <- quantile(score[upper.tri(score)], 1 - target_density)
  
  A <- (score > cutoff) * 1
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  
  # row-normalized G
  deg <- rowSums(A)
  G <- matrix(0, n, n)
  idx <- deg > 0
  G[idx, ] <- A[idx, ] / deg[idx]
  
  # simulate y = beta G y + gamma x + eps
  y <- solve(diag(n) - beta * G, gamma * x + eps)
  
  # endogenous peer outcome
  Gy <- as.vector(G %*% y)
  
  # Bramoulle-style instrument
  Gx  <- as.vector(G %*% x)
  G2x <- as.vector(G %*% G %*% x)
  
  data.frame(
    y = y,
    Gy = Gy,
    x = x,
    Gx = Gx,
    G2x = G2x,
    eps = eps,
    a = a,
    deg = deg
  ) |>
    filter(deg > 0, is.finite(Gy), is.finite(G2x))
}


# ============================================================
# 1. Run one simulation
# ============================================================

df <- simulate_network(seed = 1)

# OLS: generally biased
ols_fit <- lm(y ~ Gy + x, data = df)

# IV: instrument Gy with G2x
iv_fit <- ivreg(y ~ Gy + x | G2x + x, data = df)

summary(ols_fit)
summary(iv_fit)

c(
  beta_true = 0.4,
  beta_ols = coef(ols_fit)["Gy"],
  beta_iv = coef(iv_fit)["Gy"],
  cov_G2x_eps = cov(df$G2x, df$eps),
  cov_G2x_Gy = cov(df$G2x, df$Gy)
)

# ============================================================
# 2. Monte Carlo
# ============================================================

R <- 100
beta_true <- 0.75

res <- lapply(1:R, function(r) {
  df <- simulate_network(seed = r, beta = beta_true)
  
  ols_fit <- lm(y ~ Gy + x, data = df)
  iv_fit  <- ivreg(y ~ Gy + x | G2x + x, data = df)
  
  data.frame(
    rep = r,
    beta_ols = coef(ols_fit)["Gy"],
    beta_iv  = coef(iv_fit)["Gy"],
    
    cor_G2x_eps = cor(df$G2x, df$eps),
    cor_G2x_Gy  = cor(df$G2x, df$Gy),
    cor_Gy_eps  = cor(df$Gy,  df$eps),
    cor_Gx_eps  = cor(df$Gx,  df$eps)
  )
}) |>
  bind_rows()

summary(res)


# ============================================================
# 2. Plot estimator distributions
# ============================================================

plot_df <- res |>
  select(beta_ols, beta_iv) |>
  tidyr::pivot_longer(
    everything(),
    names_to = "estimator",
    values_to = "estimate"
  ) |>
  mutate(
    estimator = ifelse(estimator == "beta_ols", "OLS", "IV")
  )

means_df <- plot_df %>%
  group_by(estimator) %>%
  summarise(mean_est = mean(estimate, na.rm = TRUE), .groups = "drop")

ggplot(plot_df %>% filter(estimate < 2 & estimate >0), aes(x = estimate, fill = estimator)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, alpha = 0.4, position = "identity") +
  geom_vline(xintercept = beta_true,
             color = "red", linetype = "dashed") +
  geom_vline(data = means_df,
             aes(xintercept = mean_est, color = estimator),
             linetype = "dashed", linewidth = 0.7) +
  scale_fill_manual(values = c("OLS" = "steelblue", "IV" = "orange")) +
  scale_color_manual(values = c("OLS" = "steelblue", "IV" = "orange")) +
  labs(
    title = "OLS vs IV in exogenous network model",
    x = expression(hat(beta)),
    y = "Density"
  ) +
  theme_minimal()



