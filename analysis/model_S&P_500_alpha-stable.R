library(stats4)
library(fBasics)
library(stabledist)
library(here)  # Load {here} package for file path management

# Identify project location
here::i_am("analysis/model_S&P_500_alpha-stable.R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# parameter estimation using quantile
fit_stable_q <- stableFit(log_returns, type = "q") 
summary(fit_stable_q)

# Extract parameters
#params_stable <- fit_stable_q@fit$estimate
#alpha_stable <- params_stable["alpha"]
#beta_stable  <- params_stable["beta"]
#c_stable <- params_stable["gamma"]
#mu_stable <- params_stable["delta"]

# parameter estimation using mle
#fit_stable_mle <- stableFit(log_returns, type = "mle") 
#summary(fit_stable_mle)

# --- α-Stable Characteristic Function ---
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

loglik_stable <- function(alpha, c, beta, mu) {
  if (alpha <= 0 || alpha > 2 || c <= 0 || abs(beta) > 1) return(1e6)
  
  # Define truncation range manually or based on data
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

loglik_stable_sym <- function(alpha, c, mu) {
  if (alpha <= 0 || alpha > 2 || c <= 0) return(1e6)
  
  # Define truncation range manually or based on data
  a <- -0.5
  b <- 0.5
  
  # Characteristic function
  cf <- function(u) stable_cf(u, alpha, c, 0, mu)
  
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

#fit_stable_mle <- mle(loglik_stable_sym,
#                  start = list(alpha = 1.7, c = 0.01, mu = 0),
#                  method = "L-BFGS-B",
#                  lower = c(0.1, 1e-4, -1),
#                  upper = c(2, 1, 1),
#                  control = list(maxit = 1000, parscale = c(0.1, 0.001, 0.1)))

fit_stable_mle <- mle(loglik_stable,
                      start = list(alpha = 1.7, c = 0.01, beta = 0, mu = 0),
                      method = "L-BFGS-B",
                      lower = c(0.1, 1e-4, -1, -1),
                      upper = c(2, 1, 1, 1),
                      control = list(maxit = 1000, parscale = c(0.1, 0.001, 0.1, 0.1)))

# Extract parameters
params_stable <- coef(fit_stable_mle)
alpha_stable <- params_stable[1]
c_stable <- params_stable[2]
beta_stable <- params_stable[3]
mu_stable <- params_stable[4]

aic_stable <- AIC(fit_stable_mle)
cat("AIC - α-Stable:", aic_stable, "\n")


# Histogram
library(stabledist)  # for dstable

x_vals <- seq(min(log_returns), max(log_returns), length.out = 1000)
dens_vals <- dstable(x_vals, alpha_stable, 0, c_stable, mu_stable)

png(filename = here("outputs", "stable_fit01.png"), width = 2000, height = 1200, res = 300)

hist(log_returns, breaks = 100, probability = TRUE,
     col = "lightgray", main = "S&P 500 Log-Returns vs α-Stable MLE",
     xlab = "Log Return")
lines(x_vals, dens_vals, col = "red", lwd = 2)

dev.off()

# QQ plot
p_vals <- ppoints(length(log_returns))
emp_q <- sort(log_returns)
theo_q <- qstable(p_vals, alpha_stable, 0, c_stable, mu_stable)

png(filename = here("outputs", "stable_sym_QQplot.png"), width = 2000, height = 1200, res = 300)

qqplot(theo_q, emp_q,
       main = "QQ Plot: α-Stable vs Empirical",
       xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles",
       pch = 16, col = "blue")
abline(0, 1, col = "red", lwd = 2)

dev.off()

qcv <- function(x, a, b) {
  n <- length(x)
  x_sorted <- sort(x)
  
  i_start <- floor(a * n) + 1
  i_end   <- floor(b * n)
  
  x_slice <- x_sorted[i_start:i_end]
  mu_hat <- mean(x_slice)
  
  var_hat <- mean((x_slice - mu_hat)^2)
  return(var_hat)
}

compute_N1_statistic <- function(x) {
  n <- length(x)
  s1 <- qcv(x, 0.05, 0.25)
  s2 <- qcv(x, 0.25, 0.75)
  s3 <- qcv(x, 0.75, 0.95)
  s_all <- qcv(x, 0.05, 0.95)
  
  N1 <- sqrt(n) * (1.00 * s1 - 1.01 * s2 + 1.00 * s3) / s_all
  return(N1)
}

compute_N2_statistic <- function(x) {
  n <- length(x)
  s1 <- qcv(x, 0.005, 0.25)
  s2 <- qcv(x, 0.25, 0.75)
  s3 <- qcv(x, 0.75, 0.995)
  s_all <- qcv(x, 0.05, 0.995)
  
  N2 <- sqrt(n) * (0.60 * s1 - 1.61 * s2 + 0.60 * s3) / s_all
  return(N2)
}

N1_value <- compute_N1_statistic(log_returns)
N2_value <- compute_N2_statistic(log_returns)
cat("Test Statistic N1 =", N1_value, "\n")
cat("Test Statistic N2 =", N2_value, "\n")

# --- Simulation to estimate critical values ---
estimate_N1_N2_critical_values <- function(B = 10000, n = 500, alpha = 1.5, c = 1, mu = 0, N1_N2) {
  N_vals <- replicate(B, {
    x <- rstable(n, alpha = alpha, beta = 0, gamma = c, delta = mu)
    if (N1_N2){
      compute_N1_statistic(x)
    }
    else {
      compute_N2_statistic(x)
    }
  })
  
  criticals <- quantile(N_vals, probs = c(0.90, 0.95, 0.99))
  return(list(
    N_stats = N_vals,
    critical_values = criticals
  ))
}

# --- Run the estimation ---
set.seed(7914)
result <- estimate_N1_N2_critical_values(B = 10000, n = 6400, 
                                         alpha = 1.6, 
                                         c = c_stable, 
                                         mu = mu_stable,
                                         N1_N2 = FALSE)

# --- Output ---
cat("Estimated critical values (n = 6400, α = 1.6):\n")
print(result$critical_values)

# --- Optional: plot histogram of null distribution ---
hist(result$N_stats, breaks = 50, col = "lightgray",
     main = "Empirical Null Distribution of N2 Statistic",
     xlab = "N1 Statistic")
abline(v = result$critical_values, col = c("blue", "red", "darkred"), lwd = 2)
legend("topright", legend = c("90%", "95%", "99%"),
       col = c("blue", "red", "darkred"), lwd = 2)

