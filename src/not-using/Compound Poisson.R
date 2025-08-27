# ----------------------------------------------
# Simulate a Compound Poisson Process in R
# ----------------------------------------------
# This script:
# 1. Simulates a Compound Poisson process Xt = sum_{i=1}^{Nt} Zi
#    where:
#      - Nt is a Poisson process with rate lambda
#      - Zi are i.i.d. exponential jump sizes (Exp(1))
# 2. Plots the sample path of the process over time
# 3. Plots the empirical density (histogram) of the jump sizes
#    along with the theoretical exponential density
# ----------------------------------------------

# --- Parameters ---
lambda <- 3         # Poisson intensity (events per unit time)
T <- 100             # total time
dt <- 0.01          # time resolution
time_grid <- seq(0, T, by = dt)
n_steps <- length(time_grid)

# --- Simulate Jump Times and Sizes ---
set.seed(7914)
N_total <- rpois(1, lambda * T)             # total number of jumps
jump_times <- sort(runif(N_total, 0, T))    # jump times uniformly distributed
jump_sizes <- rnorm(N_total, mean = 0, sd = 1)       # jump sizes from Norm(0, 1)

# --- Build Compound Poisson Process ---
X <- numeric(n_steps)
current_jump_index <- 1
cumulative_sum <- 0

for (i in 2:n_steps) {
  t <- time_grid[i]
  while (current_jump_index <= N_total && jump_times[current_jump_index] <= t) {
    cumulative_sum <- cumulative_sum + jump_sizes[current_jump_index]
    current_jump_index <- current_jump_index + 1
  }
  X[i] <- cumulative_sum
}

# --- Plot 1: Compound Poisson Process ---
plot(time_grid, X, type = "s", lwd = 2, col = "darkgreen",
     main = "Compound Poisson Process (Norm(0, 1) jumps)",
     xlab = "Time", ylab = expression(X[t]))

# --- Plot 2: Histogram of Jump Sizes (approximate PDF) ---
hist(jump_sizes, breaks = 20, freq = FALSE, col = "skyblue",
     main = "Histogram of Jump Sizes (Z ~ Norm(0, 1))",
     xlab = "Jump Size", ylab = "Density")
curve(dnorm(x, mean = 0, sd = 1), col = "red", lwd = 2, add = TRUE)  # true density
legend("topright", legend = c("Empirical", "Theoretical"), col = c("skyblue", "red"), lwd = 2)
