# ----------------------------------------------
# Estimate alpha-stable parameters from S&P 500 data
# ----------------------------------------------

# --- Import required package ---
library(stats4)
library(fBasics)
library(stabledist)
library(here)

set.seed(7914) # For reproducibility

# Identify project location
here::i_am("analysis/model_S&P_500_alpha-stable(COS_method).R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- alpha-Stable Characteristic Function ---
stable_cf <- function(u, alpha, c, beta, mu) {
  if (alpha == 1) {
    Phi = -2 / pi * log(abs(u))
  }
  else {
    Phi = tan(pi * alpha / 2)
  }
  phi <- -c^alpha * abs(u)^alpha *
    (1 - 1i * beta * sign(u) * Phi) +
    1i * mu * u
  return(exp(phi))
}

# --- Density Approximation (COS method) ---
cos_density <- function(x_vals, cf_func, a, b, N = 1024) {
  k <- 0:(N - 1)
  u <- k * pi / (b - a)
  
  A <- Re(cf_func(u) * exp(-1i * u * a))
  A[1] <- A[1] / 2  # A_0 needs half weight
  
  sapply(x_vals, function(x) {
    sum(A * cos(u * (x - a))) * (2 / (b - a))
  })
}

# --- Log-Likelihood Function ---
loglik_stable <- function(alpha, c, beta, mu) {
  # Return a large penalty if parameters are invalid
  if (alpha <= 0 || alpha > 2 || c <= 0 || abs(beta) > 1) return(1e6)
  
  # Define truncation range for COS method
  a <- -0.5
  b <- 0.5
  
  # Characteristic function
  cf <- function(u) stable_cf(u, alpha, c, beta, mu)
  
  # COS density
  x_vals <- seq(a, b, length.out = 1000)
  dens_vals <- cos_density(x_vals, cf, a, b, N = 1024)
  pdf_func <- approxfun(x_vals, dens_vals, rule = 2)
  
  # Evaluate on data
  dens <- pdf_func(log_returns)
  
  # Penalize if density fails
  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e6)
  
  return(-sum(log(dens)))
}

# --- Fit using mle ---
fit_stable_mle <- mle(loglik_stable,
                      start = list(alpha = 1.7, c = 0.01, beta = 0, mu = 0),
                      method = "L-BFGS-B",
                      lower = c(0.1, 1e-4, -1, -1),
                      upper = c(2, 1, 1, 1),
                      control = list(maxit = 1000, parscale = c(0.1, 0.001, 0.1, 0.1)))

# --- Output fitted parameters ---
summary(fit_stable_mle)
params_stable <- coef(fit_stable_mle)
alpha_stable <- params_stable[1]
c_stable <- params_stable[2]
beta_stable <- params_stable[3]
mu_stable <- params_stable[4]

aic_stable <- AIC(fit_stable_mle)
cat("AIC - alpha-Stable:", aic_stable, "\n")


# --- Plot comparison ---
x_vals <- seq(min(log_returns), max(log_returns), length.out = 1000)
dens_vals <- dstable(x_vals, alpha_stable, 0, c_stable, mu_stable)

png(filename = here("outputs", "stable_fit01.png"), width = 2000, height = 1200, res = 300)

hist(log_returns, breaks = 150, probability = TRUE,
     col = "lightblue", main = "S&P 500 Log-Returns vs α-Stable Fit",
     xlab = "Log Return")
lines(x_vals, dens_vals, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "α-Stable Fit"), 
       col = c("lightblue", "red"), lwd = 2, cex = 1)

dev.off()

# --- Q–Q plot ---
png(filename = here("outputs", "stable_QQplot.png"), width = 2000, height = 1200, res = 300)

# Sort empirical data
log_returns_sorted <- sort(log_returns)
n <- length(log_returns_sorted)

# Compute probabilities
p_vals <- ppoints(n)

# Compute theoretical quantiles from fitted α-stable distribution
theo_q <- qstable(p_vals, alpha_stable, 0, c_stable, mu_stable)

qqplot(theo_q, log_returns_sorted,
       main = "Q-Q Plot: α-stable Fit vs Empirical Returns",
       xlab = "Theoretical Quantiles", 
       ylab = "Empirical Quantiles",
       pch = 16, col = "darkblue", cex = 0.6)
abline(0, 1, col = "red", lwd = 2)

dev.off()

