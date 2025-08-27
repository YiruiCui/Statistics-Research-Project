# ----------------------------------------------
# Simulate an Inverse Gaussian (IG) Process in R
# ----------------------------------------------
# This script:
# 1. Simulates an Inverse Gaussian process X_t with parameters mu, lambda
#    such that each increment over dt ~ IG(mu * dt, lambda * dt^2)
# 2. Plots the sample path of the IG process over time
# 3. Plots the theoretical Inverse Gaussian PDF of X_T
# ----------------------------------------------

require(statmod)  # for rinvgauss() and dinvgauss()

# --- Parameters ---
mu <- 1           # mean rate of IG distribution (mu = a / b)
lambda <- 2       # shape parameter (lambda = a^2)
T <- 10           # total time
dt <- 0.1        # time increment
time_grid <- seq(0, T, by = dt)
n_steps <- length(time_grid)

# --- Simulate IG increments ---
set.seed(7914)
increments <- rinvgauss(n_steps - 1, mean = mu * dt, shape = lambda * dt^2)
ig_process <- c(0, cumsum(increments))  # build path

# --- Plot 1: IG process path ---
plot(time_grid, ig_process, type = "s", lwd = 2, col = "darkorange",
     main = "Simulated Inverse Gaussian Process",
     xlab = "Time", ylab = expression(X[t]))

# --- Plot 2: PDF of IG increment distribution ---
x_vals <- seq(min(increments), quantile(increments, 0.99), length.out = 500)
increment_pdf <- dinvgauss(x_vals, mean = mu * dt, shape = lambda * dt^2)

plot(x_vals, increment_pdf, type = "l", lwd = 2, col = "darkblue",
     main = paste0("PDF of IG(", signif(mu * dt, 4), ", ", signif(lambda * dt^2, 4), ") Increment"),
     xlab = expression(Delta[X]), ylab = "Density")

# --- Plot 3: PDF of X_T ---
final_value <- ig_process[length(ig_process)]
x_vals <- seq(0.01, final_value * 1.5, length.out = 500)
pdf_vals <- dinvgauss(x_vals, mean = mu * T, shape = lambda * T^2)

plot(x_vals, pdf_vals, type = "l", lwd = 2, col = "darkblue",
     main = paste0("IG(", mu * T, ", ", lambda * T^2, ") PDF at T = ", T),
     xlab = expression(X[T]), ylab = "Density")
abline(v = final_value, col = "red", lty = 2)
legend("topright", legend = c("Theoretical PDF", "Simulated X_T"),
       col = c("darkblue", "red"), lty = c(1, 2), lwd = 2)
