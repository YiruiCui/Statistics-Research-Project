# ----------------------------------------------
# Estimate CGMY parameters from S&P 500 data
# ----------------------------------------------
# Please note, the code will produce warnings which can be ignored.

# --- Import required package ---
library(stats4)
library(here) 

set.seed(7914) # For reproducibility

# Identify project location
here::i_am("analysis/model_S&P_500_CGMY(COS_method).R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- CGMY Characteristic Function ---
cgmy_cf <- function(u, C, G, M, Y, t = 1) {
  psi <- C * gamma(-Y) * ((M - 1i * u)^Y - M^Y + (G + 1i * u)^Y - G^Y)
  exp(t * psi)
}

# --- CGMY with drift term Characteristic Function ---
cgmy_cf_m <- function(u, C, G, M, Y, m) {
  cgmy_cf(u, C, G, M, Y, t = 1) * exp(1i * u * m)
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
loglik_cgmy <- function(C, G, M, Y, m) {
  # Return a large penalty if parameters are invalid
  if (C <= 0 || G <= 0 || M <= 0 || Y <= 0 || Y >= 2) return(1e6)
  
  # Define truncation range for COS method
  a <- -0.5 
  b <- 0.5 
  
  # Define characteristic function with current parameters
  cf <- function(u) cgmy_cf_m(u, C, G, M, Y, m)
  
  # Estimate density using COS method
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
fit_cgmy <- mle(loglik_cgmy,
           start = list(C = 0.002, G = 23, M = 23, Y = 1.5, m = 0.01),
           method = "L-BFGS-B",
           lower = c(1e-6, 1e-3, 1e-3, 0.1, -0.1),
           upper = c(1, 100, 100, 1.99, 0.1),
           control = list(maxit = 1000, parscale = rep(0.1, 5))
          )

# --- Output fitted parameters ---
summary(fit_cgmy)
aic_cgmy <- AIC(fit_cgmy)
cat("AIC - CGMY:", aic_cgmy, "\n")

params_CGMY <- coef(fit_cgmy)
cf_fit <- function(u) cgmy_cf_m(u, 
                                params_CGMY["C"], 
                                params_CGMY["G"], 
                                params_CGMY["M"], 
                                params_CGMY["Y"], 
                                params_CGMY["m"])
x_vals <- seq(-0.5, 0.5, length.out = 10000)
dens_vals <- cos_density(x_vals, cf_fit, a = -0.5, b = 0.5)

# --- Plot comparison ---
png(filename = here("outputs", "CGMY_fit01.png"), width = 2000, height = 1200, res = 300)

hist(log_returns, breaks = 150, probability = TRUE, 
     col = "lightblue", main = "S&P 500 Log-Returns vs CGMY Fit", 
     xlab = "Log Return")
lines(x_vals, dens_vals, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "CGMY Fit"), 
       col = c("lightblue", "red"), lwd = 2, cex = 1)

dev.off()

# Compute CDF from density
dx <- diff(x_vals)[1]
cdf_vals <- cumsum(dens_vals) * dx
cdf_vals <- cdf_vals / max(cdf_vals)  # normalize to [0,1]

# Interpolation function for inverse CDF (quantile function)
cgmy_quantile <- approxfun(cdf_vals, x_vals, rule = 2)

# --- Q–Q plot ---
png(filename = here("outputs", "CGMY_QQplot.png"), width = 2000, height = 1200, res = 300)

# Sort empirical data
log_returns_sorted <- sort(log_returns)
n <- length(log_returns_sorted)

# Compute probabilities
p_vals <- ppoints(n)

# Compute theoretical quantiles from fitted CGMY distribution
theo_q <- cgmy_quantile(p_vals)

qqplot(theo_q, log_returns_sorted,
       main = "Q-Q Plot: CGMY Fit vs Empirical Returns",
       xlab = "Theoretical Quantiles", 
       ylab = "Empirical Quantiles",
       pch = 16, col = "darkblue", cex = 0.6)
abline(0, 1, col = "red", lwd = 2)

dev.off()