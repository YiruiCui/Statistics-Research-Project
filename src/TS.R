# ----------------------------------------------
# Simulate a Tempered Stable (TS) Process in R
# ----------------------------------------------
# This script:
# 1. Simulates a TS process using Rosinski's series representation
#    from Schoutens (Section 8.4.3)
# 2. Plots the process path over time
# 3. Plots the PDF of the increment approximation via histogram
# ----------------------------------------------

# --- Parameters ---
kappa <- 0.5      # 0 < kappa < 1
a <- 1            # activity
b <- 1            # tempering parameter
T <- 1000           # total time
dt <- 0.3
time_grid <- seq(0, T, by = dt)
n_steps <- length(time_grid)
K <- 10000        # truncation level of the series

# --- Generate random variables ---
set.seed(7914)
e <- rexp(K)
u <- runif(K)
u_tilde <- runif(K)

# Simulate Poisson arrival times (use exponential spacings)
b_i <- cumsum(rexp(K, rate = 1))

# --- Simulate process path ---
ts_process <- numeric(n_steps)
for (j in 2:n_steps) {
  t <- time_grid[j]
  indicators <- as.numeric(T * u < t)
  jumps <- pmin(
    (2 * (a * T / (b_i * gamma(1 - kappa)))^(1 / kappa)),
    (2 * e * u_tilde^(1 / kappa) / b^(1 / kappa))
  )
  ts_process[j] <- sum(jumps * indicators)
}

# --- Plot 1: TS process path ---
plot(time_grid, ts_process, type = "s", lwd = 2, col = "darkblue",
     main = "Simulated Tempered Stable Process",
     xlab = "Time", ylab = expression(X[t]))

# --- Plot 2: Histogram of increment size approximation ---
increments <- diff(ts_process)
hist(increments, breaks = 50, probability = TRUE, col = "lightblue",
     main = "Histogram of Increment Sizes",
     xlab = expression(Delta[X]), ylab = "Density")
