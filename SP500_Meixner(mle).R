loglik_meixner <- function(m, a, b, d) {
  if (a <= 0 || d <= 0 || abs(b) >= pi) return(1e6)  # constraints for validity
  
  dens <- meixner_pdf_m(log_returns, m, a, b, d)
  
  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e6)
  
  return(-sum(log(dens)))
}

library(stats4)

# Use moment-based estimates as starting values
start_vals <- list(m = m, a = alpha, b = beta, d = delta)

fit_meixner <- mle(loglik_meixner,
                   start = start_vals,
                   method = "L-BFGS-B",
                   lower = c(-1, 1e-4, -pi + 1e-4, 1e-4),
                   upper = c(1, 10, pi - 1e-4, 10),
                   control = list(maxit = 1000))
aic_meixner_mle <- AIC(fit_meixner)
cat("AIC - meixner(mle):", aic_meixner_mle, "\n")

params <- coef(fit_meixner)
m     <- params["m"]
alpha <- params["a"]
beta  <- params["b"]
delta <- params["d"]

# PDF overlay
x_vals <- seq(min(log_returns), max(log_returns), length.out = 1000)
pdf_vals <- meixner_pdf_m(x_vals, m, alpha, beta, delta)

hist(log_returns, breaks = 150, probability = TRUE, col = "lightblue",
     main = "S&P 500 Log-Returns vs Meixner MLE Fit", xlab = "Log-Return")
lines(x_vals, pdf_vals, col = "red", lwd = 2)

# QQ plot
cdf_vals <- meixner_cdf_m(x_vals, m, alpha, beta, delta)
meixner_qf_m <- approxfun(cdf_vals, x_vals, rule = 2)
q_theoretical <- meixner_qf_m(ppoints(length(log_returns)))
plot(q_theoretical, sort(log_returns), pch = 16, cex = 0.6,
     col = "darkblue", main = "Q–Q Plot: Meixner MLE vs Empirical Returns",
     xlab = "Theoretical Quantiles (Meixner)", ylab = "Empirical Quantiles")
abline(0, 1, col = "red", lwd = 2)
