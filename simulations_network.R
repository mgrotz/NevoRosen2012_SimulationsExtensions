set.seed(123)

simulate_network_model <- function(n = 100, p_edge = 0.06, beta = 0.4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  A_bin <- matrix(rbinom(n * n, 1, p_edge), n, n)
  A_bin[lower.tri(A_bin)] <- t(A_bin)[lower.tri(A_bin)]
  diag(A_bin) <- 0
  
  deg <- rowSums(A_bin)
  A <- matrix(0, n, n)
  nonisolated <- deg > 0
  A[nonisolated, ] <- A_bin[nonisolated, ] / deg[nonisolated]
  
  spec_rad <- max(Mod(eigen(A, only.values = TRUE)$values))
  if (beta >= 1 / spec_rad) {
    stop("beta too large: I - beta A may not be invertible")
  }
  
  epsilon <- rnorm(n)
  y <- solve(diag(n) - beta * A, epsilon)
  
  list(y = y, epsilon = epsilon, A = A, A_bin = A_bin, degree = deg)
}

sim <- simulate_network_model(seed = 123)
head(sim$y)

# OLS
# compute neighborhood average
ybar <- sim$A %*% sim$y

# OLS regression
ols_fit <- lm(sim$y ~ ybar)
summary(ols_fit)

# IV regression



