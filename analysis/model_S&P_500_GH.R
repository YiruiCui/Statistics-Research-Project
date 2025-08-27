# ----------------------------------------------
# Estimate GH parameters from S&P 500 data
# ----------------------------------------------

# --- Import required package ---
library(ghyp)
library(tidyverse)
library(here)

set.seed(7914) # For reproducibility

# Identify project location
here::i_am("analysis/model_S&P_500_GH.R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- Fit Generalized Hyperbolic distribution ---
gh_fit <- fit.ghypuv(log_returns, lambda = -0.5, symmetric = FALSE)

# --- Output fitted parameters ---
gh_params <- unlist(coef(gh_fit))
summary(gh_fit)
coef(gh_fit, type = "alpha.delta")
aic_gh <- AIC(gh_fit)
cat("AIC - GH:", aic_gh, "\n")

# --- Plot comparison ---
png(filename = here("outputs", "GH_fit01.png"), width = 2000, height = 1200, res = 300)

hist(log_returns, breaks = 150, probability = TRUE, 
     col = "lightblue", main = "S&P 500 Log-Returns vs GH Fit", 
     xlab = "Log-Return")
curve(dghyp(x, gh_fit), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "GH Fit"), col = c("lightblue", "red"), lwd = 2, cex = 1)

dev.off()

# --- Q–Q plot for GH fit ---
png(filename = here("outputs", "GH_QQplot.png"), width = 2000, height = 1200, res = 300)

# Sort empirical data
log_returns_sorted <- sort(log_returns)
n <- length(log_returns_sorted)

# Compute probabilities
p_vals <- ppoints(n)

# Compute theoretical quantiles from fitted GH distribution
q_theoretical <- qghyp(p_vals, object = gh_fit)

# Plot Q–Q
plot(q_theoretical, log_returns_sorted,
     main = "Q-Q Plot: GH Fit vs Empirical Returns",
     xlab = "Theoretical Quantiles",
     ylab = "Empirical Quantiles",
     pch = 16, col = "darkblue", cex = 0.6)
abline(0, 1, col = "red", lwd = 2)

dev.off()

