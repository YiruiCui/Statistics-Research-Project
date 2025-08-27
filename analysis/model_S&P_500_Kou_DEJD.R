# ----------------------------------------------
# Estimate Kou-DEJD parameters from S&P 500 data
# ----------------------------------------------

# --- Import required package ---
library(here)  
library(stats4)

set.seed(7914) # For reproducibility

# Identify project location
here::i_am("analysis/model_S&P_500_Kou_DEJD.R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# Negative Log-Likelihood Function
loglik_de <- function(p, eta1, eta2) {
  if (p <= 0 || p >= 1 || eta1 <= 0 || eta2 <= 0) return(1e6)
  
  right <- log_returns[log_returns >= 0]
  left <- log_returns[log_returns < 0]
  
  ll_right <- sum(log(p * eta1) - eta1 * right)
  ll_left  <- sum(log((1 - p) * eta2) + eta2 * left)
  
  return(- (ll_right + ll_left))  # negative log-likelihood
}

# --- Fit using mle ---
fit_dejd <- mle(loglik_de,
           start = list(p = 0.5, eta1 = 100, eta2 = 100),
           method = "L-BFGS-B",
           lower = c(0.001, 0.01, 0.01),
           upper = c(0.999, 1000, 1000))

# --- Output fitted parameters ---
summary(fit_dejd)
aic_dejd <- AIC(fit_dejd)
cat("AIC - DEJD:", aic_dejd, "\n")

# Extract fitted parameters
coef_fit_dejd <- coef(fit_dejd)
p <- coef_fit_dejd["p"]
eta1 <- coef_fit_dejd["eta1"]
eta2 <- coef_fit_dejd["eta2"]

# Define fitted density
de_dens <- function(y) {
  ifelse(y >= 0,
         p * eta1 * exp(-eta1 * y),
         (1 - p) * eta2 * exp(eta2 * y))
}

# --- Plot comparison ---
png(filename = here("outputs", "DEJD_fit01.png"), width = 2000, height = 1200, res = 300)

hist(log_returns, breaks = 150, freq = FALSE, 
     col = "lightblue", main = "S&P 500 Log-Returns vs Kou-DEJD Fit",
     xlab = "Log Return")
curve(de_dens(x), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Empirical", "Kou-DEJD Fit"), 
       col = c("lightblue", "red"), lwd = 2, cex = 1)
dev.off()

# --- Quantile function (inverse CDF) for Kou-DEJD ---
de_quantile <- function(p_vec, p, eta1, eta2) {
  q <- numeric(length(p_vec))
  q[p_vec < (1 - p)] <- log(p_vec[p_vec < (1 - p)] / (1 - p)) / eta2
  q[p_vec >= (1 - p)] <- -log((1 - p_vec[p_vec >= (1 - p)]) / p) / eta1
  return(q)
}

# --- Q–Q plot ---
png(filename = here("outputs", "DEJD_QQplot.png"), width = 2000, height = 1200, res = 300)

# Sort empirical data
log_returns_sorted <- sort(log_returns)
n <- length(log_returns_sorted)

# Compute probabilities
p_vals <- ppoints(n)

# Compute theoretical quantiles from fitted Kou-DEJD distribution
theo_q <- de_quantile(p_vals, p, eta1, eta2)

qqplot(theo_q, log_returns_sorted,
       main = "Q-Q Plot: Kou-DEJD Fit vs Empirical Returns",
       xlab = "Theoretical Quantiles", 
       ylab = "Empirical Quantiles",
       pch = 16, col = "darkblue", cex = 0.6)
abline(0, 1, col = "red", lwd = 2)

dev.off()