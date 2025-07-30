library(here)  # Load {here} package for file path management

# Identify project location
here::i_am("data/model_S&P_500_meixner.R")

# --- Step 1: Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Step 2: Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- Empirical moments ---
mu_emp <- mean(log_returns)
sigma2_emp <- var(log_returns)
sigma_emp <- sqrt(sigma2_emp)
mu3 <- mean((log_returns - mu_emp)^3)
mu4 <- mean((log_returns - mu_emp)^4)
skew_emp <- mu3 / sigma2_emp^(3/2)
kurt_emp <- mu4 / sigma2_emp^2

# Display
cat("Estimated Moments:\n")
cat("Mean     :", mu_emp, "\n")
cat("Variance :", sigma2_emp, "\n")
cat("Skewness :", skew_emp, "\n")
cat("Kurtosis :", kurt_emp, "\n")

# --- Meixner PDF using complex gamma from pracma ---
meixner_pdf <- function(x, alpha, beta, delta) {
  z <- complex(real = delta, imaginary = x / alpha)
  gamma_vals <- pracma::gammaz(z) 
  factor <- ((2 * cos(beta / 2))^(2 * delta)) / (2 * alpha * pi * gamma(2 * delta))
  density <- factor * exp(beta * x / alpha) * (Mod(gamma_vals)^2)
  return(density)
}

# --- Meixner PDF with location parameter m ---
meixner_pdf_m <- function(x, m, a, b, d) {
  z <- (x - m) / a
  density <- meixner_pdf(z, alpha = 1, beta = b, delta = d)
  return(density / a)
}

# CDF by numerical integration of the PDF
meixner_cdf_m <- function(x, m, a, b, d) {
  sapply(x, function(xi) {
    integrate(meixner_pdf_m, lower = -100, upper = xi, 
              m = m, a = a, b = b, d = d, rel.tol = 1e-6)$value
  })
}

# --- Estimate Meixner parameters ---
delta <- 1/ (kurt_emp - skew_emp^2 - 3)  
beta <- sign(skew_emp) * acos(2 - delta * (kurt_emp - 3))  
alpha <- sigma_emp * sqrt((cos(beta)+1)/delta)
m <- mu_emp - alpha * delta * tan(beta/2)

# Generate a lookup table for quantile inversion
x_grid <- seq(min(log_returns) - 0.05, max(log_returns) + 0.05, length.out = 1000)
cdf_vals <- meixner_cdf_m(x_grid, m, alpha, beta, delta)

# Approximate inverse CDF (quantile function)
meixner_qf_m <- approxfun(cdf_vals, x_grid, rule = 2)  # Use interpolation

loglik_meixner <- sum(log(meixner_pdf_m(log_returns, m, alpha, beta, delta)))
k_meixner <- 4  # m, a, b, d
aic_meixner <- -2 * loglik_meixner + 2 * k_meixner
cat("AIC - Meixner:", aic_meixner, "\n")


# --- Step 5: Plot comparison ---
png(filename = here("outputs", "Meixner_fit01.png"), width = 2000, height = 1200, res = 300)

hist(log_returns, breaks = 150, probability = TRUE, col = "lightblue", main = "S&P 500 Log-Returns vs Meixner Fit", xlab = "Log-Return")

# Compute Meixner PDF over the range of histogram
x_vals <- seq(min(log_returns), max(log_returns), length.out = 1000)
pdf_vals <- meixner_pdf_m(x_vals, m, alpha, beta, delta)

# Overlay the Meixner PDF
lines(x_vals, pdf_vals, col = "red", lwd = 2)

legend("topright", legend = c("Empirical", "Meixner Fit"), col = c("lightblue", "red"), lwd = 2, cex = 0.5)

dev.off()

# Sort empirical returns
log_returns_sorted <- sort(log_returns)
n <- length(log_returns_sorted)

# Empirical probabilities
p_vals <- ppoints(n)

# Theoretical Meixner quantiles
q_theoretical_Meixner <- meixner_qf_m(p_vals)

# Plot QQ
png(filename = here("outputs", "Meixner_QQplot.png"), width = 2000, height = 1200, res = 300)

plot(q_theoretical_Meixner, log_returns_sorted,
     main = "Q–Q Plot: Meixner Fit vs Empirical Returns",
     xlab = "Theoretical Quantiles (Meixner)",
     ylab = "Empirical Quantiles",
     pch = 16, col = "darkblue", cex = 0.6)
abline(0, 1, col = "red", lwd = 2)

dev.off()
