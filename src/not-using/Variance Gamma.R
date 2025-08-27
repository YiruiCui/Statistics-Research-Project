# ----------------------------------------------
# Simulate a Variance Gamma (VG) Process in R
# ----------------------------------------------
# This script:
# 1. Simulates a VG process: X_t = theta * G_t + sigma * W_{G_t}
#    where G_t ~ Gamma(1 / nu, scale = nu)
# 2. Plots the sample path of the VG process
# 3. Plots the PDF of the VG increment distribution
# ----------------------------------------------

# --- Parameters ---
theta <- 0.5       # drift of the Brownian motion
sigma <- 1         # volatility of Brownian motion
nu <- 0.2          # variance of the Gamma process (subordinator)
T <- 10
dt <- 0.01
time_grid <- seq(0, T, by = dt)
n_steps <- length(time_grid)

# --- Simulate Gamma subordinator G_t ---
set.seed(7914)
gamma_increments <- rgamma(n_steps - 1, shape = 1 / nu, scale = nu)
G_t <- c(0, cumsum(gamma_increments))

# --- Simulate Brownian motion W_{G_t} ---
W_Gt <- c(0, cumsum(rnorm(n_steps - 1, mean = 0, sd = sqrt(diff(G_t)))))

# --- Construct VG process ---
VG_process <- theta * G_t + sigma * W_Gt

# --- Plot 1: VG process path ---
plot(time_grid, VG_process, type = "s", lwd = 2, col = "darkorange",
     main = "Simulated Variance Gamma Process",
     xlab = "Time", ylab = expression(X[t]))

# --- Plot 2: Empirical PDF of VG increment ---
VG_increments <- diff(VG_process)
x_vals <- seq(min(VG_increments), quantile(VG_increments, 0.99), length.out = 500)
increment_density <- density(VG_increments, from = min(x_vals), to = max(x_vals))

plot(increment_density$x, increment_density$y, type = "l", lwd = 2, col = "blue",
     main = "Empirical PDF of VG Increment",
     xlab = expression(Delta[X]), ylab = "Density")
