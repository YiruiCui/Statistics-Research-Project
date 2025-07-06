# Negative Log-Likelihood Function
loglik_de <- function(p, eta1, eta2) {
  if (p <= 0 || p >= 1 || eta1 <= 0 || eta2 <= 0) return(1e6)
  
  right <- log_returns[log_returns >= 0]
  left <- log_returns[log_returns < 0]
  
  ll_right <- sum(log(p * eta1) - eta1 * right)
  ll_left  <- sum(log((1 - p) * eta2) + eta2 * left)
  
  return(- (ll_right + ll_left))  # negative log-likelihood
}

# Fit using mle
library(stats4)
fit <- mle(loglik_de,
           start = list(p = 0.5, eta1 = 100, eta2 = 100),
           method = "L-BFGS-B",
           lower = c(0.001, 0.01, 0.01),
           upper = c(0.999, 1000, 1000))

summary(fit)

# Extract fitted values
coef_fit <- coef(fit)
p <- coef_fit["p"]
eta1 <- coef_fit["eta1"]
eta2 <- coef_fit["eta2"]

# Define fitted density
de_dens <- function(y) {
  ifelse(y >= 0,
         p * eta1 * exp(-eta1 * y),
         (1 - p) * eta2 * exp(eta2 * y))
}

# Plot histogram and fitted density
hist(log_returns, breaks = 100, freq = FALSE, col = "gray80",
     main = "Fit of Double Exponential to S&P 500 Log Returns", xlab = "Log return")
curve(de_dens(x), col = "red", lwd = 2, add = TRUE)
