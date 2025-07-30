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
fit_dejd <- mle(loglik_de,
           start = list(p = 0.5, eta1 = 100, eta2 = 100),
           method = "L-BFGS-B",
           lower = c(0.001, 0.01, 0.01),
           upper = c(0.999, 1000, 1000))

summary(fit_dejd)
aic_dejd <- AIC(fit_dejd)
cat("AIC - DEJD:", aic_dejd, "\n")


# Extract fitted values
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

# Plot histogram and fitted density
png(filename = here("outputs", "DEJD_fit01.png"), width = 2000, height = 1200, res = 300)

hist(log_returns, breaks = 100, freq = FALSE, col = "gray80",
     main = "Fit of Double Exponential to S&P 500 Log Returns", xlab = "Log return")
curve(de_dens(x), col = "red", lwd = 2, add = TRUE)

dev.off()

# --- Quantile function (inverse CDF) for asymmetric double exponential ---
de_quantile <- function(p_vec, p, eta1, eta2) {
  q <- numeric(length(p_vec))
  q[p_vec < (1 - p)] <- log(p_vec[p_vec < (1 - p)] / (1 - p)) / eta2
  q[p_vec >= (1 - p)] <- -log((1 - p_vec[p_vec >= (1 - p)]) / p) / eta1
  return(q)
}

# Sort log-returns
empirical_q <- sort(log_returns)

# Generate theoretical quantiles based on fitted parameters
p_vec <- ppoints(length(empirical_q))  # uniform [0,1]
theoretical_q <- de_quantile(p_vec, p, eta1, eta2)

# Plot QQ plot
png(filename = here("outputs", "DEJD_QQplot.png"), width = 2000, height = 1200, res = 300)

qqplot(theoretical_q, empirical_q,
       main = "QQ Plot: Double Exponential vs Empirical Log Returns",
       xlab = "Theoretical Quantiles (Double Exponential)",
       ylab = "Empirical Quantiles (S&P 500)")
abline(0, 1, col = "red", lwd = 2)

dev.off()