# ----------------------------------------------
# Estimate ghyp parameters from S&P 500 data
# ----------------------------------------------

# --- Install required package ---
if (!require(ghyp)) install.packages("ghyp")
library(ghyp)
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
summary(gh_fit)

# --- Step 5: Plot comparison ---
hist(log_returns, breaks = 100, probability = TRUE, col = "lightblue",
     main = "S&P 500 Log-Returns vs GH/CGMY Fit", xlab = "Log-Return")
curve(dghyp(x, gh_fit), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "CGMY Fit"), col = c("lightblue", "red"), lwd = 2)