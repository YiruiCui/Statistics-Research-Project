# ----------------------------------------------
# Simulate CGMY process using Madan–Yor method (Poirot & Tankov, Section 4.2)
# ----------------------------------------------

# --- Required packages ---
require(pracma) # for gamma functions

# --- Parameters ---
alpha <- 0.5       # CGMY parameter α (0 < α < 2)
C <- 0.5           # jump activity (merged into c)
lambda_p <- 3.5    # λ+ > 0
lambda_m <- 2.0    # λ− > 0
A <- (lambda_m - lambda_p) / 2
B <- (lambda_m + lambda_p) / 2
K <- 1             # arbitrarily set Lévy prefactor
epsilon <- 1e-4    # jump size truncation threshold
T <- 1             # simulation horizon
n_paths <- 100000   # number of samples

# --- Step 1: Expected drift from truncated small jumps ---
drift <- K * epsilon^(1 - alpha / 2) / (1 - alpha / 2)

# --- Step 2: Simulate jump sizes from ν_0(t) = K * t^{-1 - α/2}, t > ε ---
simulate_jumps <- function(lambda, T, alpha, eps, K) {
  N <- rpois(1, lambda * T)
  U <- runif(N)
  jump_sizes <- ( (1 - U) * eps^(alpha/2) )^(-2/alpha)
  return(jump_sizes)
}

# --- Accurate envelope using D_{-α}(z) ---
f_envelope <- function(y, A, B, alpha) {
  z <- B * sqrt(y)
  D_approx <- z^(-alpha) * exp(-z^2 / 4)
  
  pre_factor <- 2 ^ (alpha / 2) * gamma(alpha / 2 + 0.5) / sqrt(pi)
  f <- pre_factor * exp(0.5 * A^2 * y + 0.5 * B^2 * y) * D_approx
  return(pmin(f, 1))  # bound f ≤ 1
}

# --- Step 4: Accept jumps using rejection method ---
sample_CGMY_path <- function() {
  jump_sizes <- simulate_jumps(lambda = K, T, alpha, epsilon, K)
  U <- runif(length(jump_sizes))
  accept_probs <- f_envelope(jump_sizes, A, B, alpha)
  accepted_jumps <- jump_sizes[U < accept_probs]
  
  Z_t <- T * drift + sum(accepted_jumps)
  X_t <- A * Z_t + sqrt(Z_t) * rnorm(1)
  return(X_t)
}

# --- Simulate many samples ---
set.seed(42)
cgmy_samples <- replicate(n_paths, sample_CGMY_path())

# --- Plot results ---
hist(cgmy_samples, breaks = 100, freq = FALSE, col = "lightblue",
     main = "Simulated CGMY Increments via Madan–Yor (Poirot-Tankov)",
     xlab = "CGMY increment")
