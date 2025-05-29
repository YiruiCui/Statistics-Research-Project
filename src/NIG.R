# ----------------------------------------------
# Simulate a Normal Inverse Gaussian (NIG) Process in R
# ----------------------------------------------
# This script:
# 1. Simulates the NIG process: X_t = beta * delta^2 * I_t + delta * W_{I_t}
#    where I_t ~ IG(1, b) with b = delta * sqrt(alpha^2 - beta^2)
# 2. Plots the process path
# 3. Plots the PDF of one increment
# ----------------------------------------------

require(statmod)  # for rinvgauss()

# --- Parameters ---
alpha <- 1.0      
beta <- 0.2     # |beta| < alpha
delta <- 0.5
T <- 100
dt <- 0.1
time_grid <- seq(0, T, by = dt)
n_steps <- length(time_grid)

# Derived parameter
b <- delta * sqrt(alpha^2 - beta^2)

# --- Simulate IG subordinator I_t ~ IG(1, b) ---
set.seed(7914)
IG_mean <- 1 / b
IG_shape <- 1
IG_increments <- rinvgauss(n_steps - 1, mean = IG_mean, shape = IG_shape)
I_t <- c(0, cumsum(IG_increments))  # cumulative IG subordinator

# --- Simulate Brownian motion W_{I_t} ---
W_I_t <- c(0, cumsum(rnorm(n_steps - 1, mean = 0, sd = sqrt(diff(I_t)))))

# --- Construct NIG process (book formula) ---
X_t <- beta * delta^2 * I_t + delta * W_I_t

# --- Plot 1: NIG process path ---
plot(time_grid, X_t, type = "s", lwd = 2, col = "darkblue",
     main = "Simulated Normal Inverse Gaussian (NIG) Process",
     xlab = "Time", ylab = expression(X[t]))

# --- Plot 2: Empirical PDF of NIG increment ---
increments <- diff(X_t)
x_vals <- seq(min(increments), quantile(increments, 0.99), length.out = 500)
increment_density <- density(increments, from = min(x_vals), to = max(x_vals))

plot(increment_density$x, increment_density$y, type = "l", lwd = 2, col = "forestgreen",
     main = "Empirical PDF of NIG Increment",
     xlab = expression(Delta[X]), ylab = "Density")
