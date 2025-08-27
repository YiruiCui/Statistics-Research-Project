# ----------------------------------------------
# Simulate a Generalized Hyperbolic (GH) Process using ghyp package
# ----------------------------------------------

require(ghyp)

# --- Define GH parameters ---
# lambda: tail parameter
# alpha: determines tail heaviness
# mu: location
# sigma: volatility
# gamma: skewness

gh_model <- ghyp(lambda = 1, alpha = 1.5, mu = 0, sigma = 1, gamma = 0.5)

# --- Time setup ---
T <- 100
dt <- 0.1
n_steps <- length(seq(0, T, by = dt))

# --- Simulate GH increments and build process ---
set.seed(7914)
gh_increments <- rghyp(n_steps - 1, object = gh_model)
gh_process <- c(0, cumsum(gh_increments))
time_grid <- seq(0, T, by = dt)

# --- Plot 1: GH Process Path ---
plot(time_grid, gh_process, type = "s", col = "steelblue", lwd = 2,
     main = "Generalized Hyperbolic Process (ghyp)",
     xlab = "Time", ylab = expression(X[t]))

# --- Plot 2: Empirical PDF ---
plot(density(gh_increments), type = "l", col = "darkred", lwd = 2,
     main = "Empirical PDF of GH Increments",
     xlab = expression(Delta[X]), ylab = "Density")
