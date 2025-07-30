library(stats4)

# --- CGMY Characteristic Function ---
cgmy_cf <- function(u, C, G, M, Y, t = 1) {
  psi <- C * gamma(-Y) * ((M - 1i * u)^Y - M^Y + (G + 1i * u)^Y - G^Y)
  exp(t * psi)
}

# --- CGMY with drift term Characteristic Function ---
cgmy_cf_m <- function(u, C, G, M, Y, m) {
  cgmy_cf(u, C, G, M, Y, t = 1) * exp(1i * u * m)
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


loglik_cgmy <- function(C, G, M, Y, m) {
  # Return a large penalty if parameters are invalid
  if (C <= 0 || G <= 0 || M <= 0 || Y <= 0 || Y >= 2) return(1e6)
  
  # Compute truncation range [a, b] using current CGMY parameters' cumulants
  #range <- get_cumulant_range(C, G, M, Y, L = 10)
  a <- -0.5 #range["a"]
  b <- 0.5 #range["b"]
  
  # Define characteristic function with current parameters
  cf <- function(u) cgmy_cf_m(u, C, G, M, Y, m)
  
  # Estimate density using COS method
  x_vals <- seq(a, b, length.out = 1000)
  dens_vals <- cos_density(x_vals, cf, a, b, N = 1024)
  pdf_func <- approxfun(x_vals, dens_vals, rule = 2)  # interpolation
  if (min(dens_vals) < 0){
    #print(C)
    #print(G)
    #print(M)
    #print(Y)
    #print(min(dens_vals))
  }
  dens <- pdf_func(log_returns)
  
  # If any density is non-positive or non-finite, return penalty
  if (any(!is.finite(dens)) || any(dens <= 0)) return(1e6)
  
  return(-sum(log(dens)))
  # Compute log-likelihood over observed returns
  #ll_vals <- log(pdf_func(log_returns) + 1e-10)  # avoid log(0)
  #return(-sum(ll_vals))  # negative log-likelihood for MLE
}

# --- fit CGMY ---
fit_cgmy <- mle(loglik_cgmy,
           start = list(C = 0.002, G = 23, M = 23, Y = 1.5, m = 0.01),
           method = "L-BFGS-B",
           lower = c(1e-6, 1e-3, 1e-3, 0.1, -0.1),
           upper = c(1, 100, 100, 1.99, 0.1),
           control = list(maxit = 1000, parscale = rep(0.1, 5))
          )
summary(fit_cgmy)
aic_cgmy <- AIC(fit_cgmy)
cat("AIC - CGMY:", aic_cgmy, "\n")


# Extract parameters
params <- coef(fit_cgmy)
cf_fit <- function(u) cgmy_cf_m(u, params["C"], params["G"], params["M"], params["Y"], params["m"])
x_vals <- seq(-0.5, 0.5, length.out = 1000)
dens_fit <- cos_density(x_vals, cf_fit, a = -0.5, b = 0.5)

# Plot
png(filename = here("outputs", "CGMY_fit01.png"), width = 2000, height = 1200, res = 300)

hist(log_returns, breaks = 200, probability = TRUE, col = "gray", border = NA,
     main = "CGMY Fit to S&P 500 Log Returns", xlab = "Log return")
lines(x_vals, dens_fit, col = "red", lwd = 2)

dev.off()

# Grid and density
a <- -1
b <- 1
x_vals <- seq(a, b, length.out = 10000)
dens_vals <- cos_density(x_vals, cf_fit, a, b, N = 1024)

# Compute CDF from density
dx <- diff(x_vals)[1]
cdf_vals <- cumsum(dens_vals) * dx
cdf_vals <- cdf_vals / max(cdf_vals)  # normalize to [0,1]

# Interpolation function for inverse CDF (quantile function)
cgmy_quantile <- approxfun(cdf_vals, x_vals, rule = 2)

# Compute empirical quantiles
empirical_q <- sort(log_returns)

# Uniform quantiles (same length)
p_vals <- ppoints(length(empirical_q))

# Theoretical quantiles from CGMY
theoretical_q <- cgmy_quantile(p_vals)

library(ggplot2)
qq_data <- data.frame(Theoretical = theoretical_q, Empirical = empirical_q)

png(filename = here("outputs", "CGMY_QQplot.png"), width = 2000, height = 1200, res = 300)

ggplot(qq_data, aes(x = Theoretical, y = Empirical)) +
  geom_point(color = "steelblue", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "QQ Plot: Fitted CGMY vs Empirical Log-Returns",
       x = "Theoretical Quantiles (CGMY)", y = "Empirical Quantiles") +
  theme_minimal()

dev.off()