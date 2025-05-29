# ----------------------------------------------
# Simulate a Generalized Inverse Gaussian (GIG) Process in R
# ----------------------------------------------
# This script:
# 1. Simulates a GIG process X_t with parameters (lambda, chi, psi)
#    where each increment over dt ~ GIG(lambda, chi, psi)
# 2. Plots the sample path of the GIG process over time
# 3. Plots the PDF of the increment distribution over dt
# ----------------------------------------------

require(GeneralizedHyperbolic)

# --- Parameters ---
lambda <- 0.5       # shape parameter
chi <- 1.0          # > 0
psi <- 1.0          # > 0
T <- 10             # total time
dt <- 0.1          # time increment
time_grid <- seq(0, T, by = dt)
n_steps <- length(time_grid)

# --- Simulate GIG increments over small time steps ---
set.seed(7914)
gig_increments <- rgig(n_steps - 1, lambda = lambda, chi = chi, psi = psi)
gig_process <- c(0, cumsum(gig_increments))  # build the path

# --- Plot 1: GIG process path ---
plot(time_grid, gig_process, type = "s", lwd = 2, col = "darkcyan",
     main = "Simulated GIG Process",
     xlab = "Time", ylab = expression(X[t]))

# --- Plot 2: PDF of one increment GIG(lambda, chi, psi) ---
x_vals <- seq(min(gig_increments), quantile(gig_increments, 0.99), length.out = 500)
increment_pdf <- dgig(x_vals, lambda = lambda, chi = chi, psi = psi)

plot(x_vals, increment_pdf, type = "l", lwd = 2, col = "darkgreen",
     main = paste0("PDF of GIG(", round(lambda, 5), ", ", chi, ", ", psi, ") Increment"),
     xlab = expression(Delta[X]), ylab = "Density")
