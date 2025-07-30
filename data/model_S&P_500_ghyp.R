# ----------------------------------------------
# Estimate ghyp parameters from S&P 500 data
# ----------------------------------------------

# --- Install required package ---
library(ghyp)
library(tidyverse)
library(here)  # Load {here} package for file path management

# Identify project location
here::i_am("data/model_S&P_500_ghyp.R")

# --- Step 1: Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Step 2: Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- Step 3: Fit Generalized Hyperbolic distribution ---
# Set lambda = -0.5 for NIG subclass
gh_fit <- fit.ghypuv(log_returns, lambda = -0.5, symmetric = FALSE)

# --- Step 4: Output fitted parameters ---
gh_params <- unlist(coef(gh_fit))
summary(gh_fit)
aic_gh <- AIC(gh_fit)
cat("AIC - GH:", aic_gh, "\n")

# --- Step 5: Plot comparison ---
png(filename = here("outputs", "GH_fit01.png"), width = 2000, height = 1200, res = 300)

hist(log_returns, breaks = 150, probability = TRUE, col = "lightblue", main = "S&P 500 Log-Returns vs GH Fit", xlab = "Log-Return")
curve(dghyp(x, gh_fit), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "GH Fit"), col = c("lightblue", "red"), lwd = 2, cex = 1)

dev.off()

# --- Step 6: Q–Q plot for GH fit ---
png(filename = here("outputs", "GH_QQplot.png"), width = 2000, height = 1200, res = 300)

# Sort empirical data
log_returns_sorted <- sort(log_returns)
n <- length(log_returns_sorted)

# Compute probabilities
p_vals <- ppoints(n)  # Generates evenly spaced probabilities

# Compute theoretical quantiles from fitted GH distribution
q_theoretical <- qghyp(p_vals, object = gh_fit)

# Plot Q–Q
plot(q_theoretical, log_returns_sorted,
     main = "Q–Q Plot: GH Fit vs Empirical Returns",
     xlab = "Theoretical Quantiles (GH)",
     ylab = "Empirical Quantiles",
     pch = 16, col = "darkblue", cex = 0.6)
abline(0, 1, col = "red", lwd = 2)

dev.off()

