# ----------------------------------------------
# Estimate Meixner parameters from S&P 500 data
# ----------------------------------------------

# --- Import required package ---
library(here)
library(pracma)
library(stats4)

set.seed(7914) # For reproducibility

# Identify project location
here::i_am("analysis/model_S&P_500_Meixner.R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- Meixner PDF ---
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

# --- Meixner CDF ---
meixner_cdf_m <- function(x, m, a, b, d) {
  sapply(x, function(xi) {
    integrate(meixner_pdf_m, lower = -100, upper = xi, 
              m = m, a = a, b = b, d = d, rel.tol = 1e-6)$value
  })
}

# --- Log-Likelihood Function ---
loglik_meixner <- function(m, a, b, d) {
  if (a <= 0 || d <= 0 || abs(b) >= pi) return(1e6)  # constraints for validity
  
  dens <- meixner_pdf_m(log_returns, m, a, b, d)
  
  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e6)
  
  return(-sum(log(dens)))
}

# --- Fit using mle ---
fit_meixner <- mle(loglik_meixner,
                   start = list(m = mean(log_returns), a = 0.01, b = -0.1, d = 0.1),
                   method = "L-BFGS-B",
                   lower = c(-1, 1e-4, -pi + 1e-4, 1e-4),
                   upper = c(1, 10, pi - 1e-4, 10),
                   control = list(maxit = 1000))

# --- Output fitted parameters ---
summary(fit_meixner)
aic_meixner <- AIC(fit_meixner)
cat("AIC - meixner:", aic_meixner, "\n")

# Extract fitted parameters
params <- coef(fit_meixner)
m_Meixner     <- params["m"]
alpha_Meixner <- params["a"]
beta_Meixner  <- params["b"]
delta_Meixner <- params["d"]

# --- Plot comparison ---
png(filename = here("outputs", "Meixner_fit01.png"), width = 2000, height = 1200, res = 300)

x_vals <- seq(min(log_returns), max(log_returns), length.out = 1000)
pdf_vals <- meixner_pdf_m(x_vals, m_Meixner, alpha_Meixner, beta_Meixner, delta_Meixner)

hist(log_returns, breaks = 150, probability = TRUE, 
     col = "lightblue", main = "S&P 500 Log-Returns vs Meixner Fit", 
     xlab = "Log-Return")
lines(x_vals, pdf_vals, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "Meixner Fit"), col = c("lightblue", "red"), lwd = 2, cex = 1)

dev.off()

# --- Q–Q plot for GH fit ---
png(filename = here("outputs", "Meixner_QQplot.png"), width = 2000, height = 1200, res = 300)

# Sort empirical returns
log_returns_sorted <- sort(log_returns)
n <- length(log_returns_sorted)

# Empirical probabilities
p_vals <- ppoints(n)

# Theoretical Meixner quantiles
cdf_vals <- meixner_cdf_m(x_vals, m_Meixner, alpha_Meixner, beta_Meixner, delta_Meixner)
meixner_qf_m <- approxfun(cdf_vals, x_vals, rule = 2)
q_theoretical <- meixner_qf_m(p_vals)

# Plot Q–Q
plot(q_theoretical, log_returns_sorted,
     main = "Q-Q Plot: GH Fit vs Empirical Returns",
     xlab = "Theoretical Quantiles",
     ylab = "Empirical Quantiles",
     pch = 16, col = "darkblue", cex = 0.6)
abline(0, 1, col = "red", lwd = 2)

dev.off()
