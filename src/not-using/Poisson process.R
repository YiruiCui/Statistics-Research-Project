# ----------------------------------------------
# Simulate a Poisson Process in R
# ----------------------------------------------
# This script:
# 1. Simulates a standard Poisson process Nt with intensity lambda
#    using small time increments and Poisson-distributed jumps
# 2. Plots the sample path of Nt over time
# 3. Plots the probability mass function (PMF) of Nt at final time T
# ----------------------------------------------

# --- Parameters ---
lambda <- 5       # intensity (events per unit time)
T <- 100           # total time
dt <- 0.01        # time step
time_grid <- seq(0, T, by = dt)
n_steps <- length(time_grid)

# --- Simulate Poisson increments ---
set.seed(7914)
increments <- rpois(n_steps - 1, lambda * dt)
poisson_process <- c(0, cumsum(increments))  # build path

# --- Plot 1: Poisson process sample path ---
plot(time_grid, poisson_process, type = "s", col = "blue", lwd = 2,
     main = "Simulated Poisson Process",
     xlab = "Time", ylab = expression(N[t]))

# --- Plot 2: PMF at time T ---
k_values <- 0:ceiling(1.2 * lambda * T)
pmf_values <- dpois(k_values, lambda * T)
barplot(pmf_values, names.arg = k_values,
        main = paste0("Poisson(", lambda * T, ") PMF at Time T = ", T),
        xlab = "k (Number of Events)", ylab = "P(N(T) = k)",
        col = "orange", border = NA)
