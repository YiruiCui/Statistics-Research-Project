# --- CGMY Characteristic Function ---
cgmy_cf <- function(u, C, G, M, Y, t = 1) {
  psi <- C * gamma(-Y) * ((M - 1i * u)^Y - M^Y + (G + 1i * u)^Y - G^Y)
  exp(t * psi)
}

# --- COS Density Approximation ---
cos_density <- function(x_vals, cf_func, a, b, N = 1024) {
  k <- 0:(N - 1)
  u <- k * pi / (b - a)
  
  A <- Re(cf_func(u) * exp(-1i * u * a))
  A[1] <- A[1] / 2  # A_0 needs half weight
  
  sapply(x_vals, function(x) {
    sum(A * cos(u * (x - a))) * (2 / (b - a))
  })
}

# --- Choose Truncation Interval Based on 4th Cumulants ---
get_cumulant_range <- function(C, G, M, Y, L = 10) {
  # Check if gamma terms are valid
  if (Y >= 2 || Y <= 0) return(c(-0.5, 0.5))  # fallback
  
  c1 <- C * gamma(1 - Y) * (M^(Y - 1) - G^(Y - 1))
  c2 <- C * gamma(2 - Y) * (M^(Y - 2) + G^(Y - 2))
  c4_term <- tryCatch(C * gamma(4 - Y) * (M^(Y - 4) + G^(Y - 4)),
                      error = function(e) return(Inf))
  
  # Protect against invalid sqrt
  if (!is.finite(c1) || !is.finite(c2) || !is.finite(c4_term)) return(c(-0.5, 0.5))
  
  delta <- L * sqrt(c2 + sqrt(abs(c4_term)))
  if (!is.finite(delta)) return(c(-0.5, 0.5))
  
  return(c(c1 - delta, c1 + delta))
}


loglik_cgmy <- function(C, G, M, Y) {
  # Return a large penalty if parameters are invalid
  if (C <= 0 || G <= 0 || M <= 0 || Y <= 0 || Y >= 2) return(1e6)
  
  # Compute truncation range [a, b] using current CGMY parameters' cumulants
  #range <- get_cumulant_range(C, G, M, Y, L = 10)
  a <- -0.5 #range["a"]
  b <- 0.5 #range["b"]
  
  # Define characteristic function with current parameters
  cf <- function(u) cgmy_cf(u, C, G, M, Y)
  
  # Estimate density using COS method
  x_vals <- seq(a, b, length.out = 1000)
  dens_vals <- cos_density(x_vals, cf, a, b, N = 256)
  pdf_func <- approxfun(x_vals, dens_vals, rule = 2)  # interpolation
  
  dens <- pdf_func(log_returns)
  
  # If any density is non-positive or non-finite, return penalty
  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e6)
  
  return(-sum(log(dens)))
  # Compute log-likelihood over observed returns
  #ll_vals <- log(pdf_func(log_returns) + 1e-10)  # avoid log(0)
  return(-sum(ll_vals))  # negative log-likelihood for MLE
}

# --- fit CGMY ---
fit <- mle(loglik_cgmy,
           start = list(C = 0.02, G = 0.1, M = 5, Y = 1.3),
           method = "L-BFGS-B",
           lower = c(1e-6, 1e-3, 1e-3, 0.1),
           upper = c(1, 10, 20, 1.99))

summary(fit)

# Extract parameters
params <- coef(fit)
cf_fit <- function(u) cgmy_cf(u, params["C"], params["G"], params["M"], params["Y"])
x_vals <- seq(-0.5, 0.5, length.out = 1000)
dens_fit <- cos_density(x_vals, cf_fit, a = -0.5, b = 0.5)

# Plot
hist(log_returns, breaks = 100, probability = TRUE, col = "gray", border = NA,
     main = "CGMY Fit to S&P 500 Log Returns", xlab = "Log return")
lines(x_vals, dens_fit, col = "red", lwd = 2)
