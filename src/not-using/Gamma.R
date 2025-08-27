# ----------------------------------------------
# Simulate a Gamma Process in R
# ----------------------------------------------
# This script:
# 1. Simulates a Gamma process G_t where each increment over dt
#    follows a Gamma distribution with shape = a * dt and rate = b
# 2. Plots the sample path of the Gamma process over time
# 3. Plots the probability density function (PDF) of G_T
# ----------------------------------------------

# --- Parameters ---
a <- 2           # shape rate per unit time (scale of activity)
b <- 1           # rate parameter of Gamma distribution
T <- 10          # total time
dt <- 0.1       # time increment
time_grid <- seq(0, T, by = dt)
n_steps <- length(time_grid)

# --- Simulate Gamma increments ---
set.seed(7914)
increments <- rgamma(n_steps - 1, shape = a * dt, rate = b)
gamma_process <- c(0, cumsum(increments))  # build path

# --- Plot 1: Gamma process path ---
plot(time_grid, gamma_process, type = "s", lwd = 2, col = "purple",
     main = "Simulated Gamma Process",
     xlab = "Time", ylab = expression(G[t]))


# --- Plot 2: PDF of Gamma increment ---
x_vals <- seq(0, quantile(increments, 0.99), length.out = 500)
increment_pdf <- dgamma(x_vals, shape = a * dt, rate = b)

plot(x_vals, increment_pdf, type = "l", lwd = 2, col = "darkred",
     main = paste0("PDF of Gamma(", signif(a * dt, 4), ", ", b, ") Increment"),
     xlab = expression(Delta[G]), ylab = "Density")

# --- Plot 3: Density of G_T ---
final_value <- gamma_process[length(gamma_process)]
x_vals <- seq(0, final_value * 1.5, length.out = 500)
pdf_vals <- dgamma(x_vals, shape = a * T, rate = b)

plot(x_vals, pdf_vals, type = "l", lwd = 2, col = "darkred",
     main = paste0("Gamma(", a * T, ", ", b, ") PDF at T = ", T),
     xlab = expression(G[T]), ylab = "Density")
abline(v = final_value, col = "blue", lty = 2)
legend("topright", legend = c("Theoretical PDF", "Simulated G_T"),
       col = c("darkred", "blue"), lty = c(1, 2), lwd = 2)
