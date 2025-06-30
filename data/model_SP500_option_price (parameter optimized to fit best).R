# --- Load libraries ---
if (!require(here)) install.packages("here")
if (!require(readr)) install.packages("readr")
if (!require(ggplot2)) install.packages("ggplot2")
library(here)
library(readr)
library(ggplot2)

# Identify project location
here::i_am("data/model_SP500_option_price_Esscher.R")

# --- Load market data ---
data <- readr::read_csv(here("data", "SP500_10.06.2025_options.csv"))
data <- data[1:102,]
market_prices <- data$Last
K_vec <- data$Strike
T_days <- data$T
T_years <- T_days / 252  # convert to years (252 trading day)

# --- Define parameters  ---
r <- 0.041     # risk-free rate
q <- 0.015
S0 <- 6038.81    # current index price 

# --- Risk-neutral Meixner CF generator (returns a function of u) ---
make_meixner_cf <- function(alpha, beta, delta, S0, r, q, T) {
  mu_rn <- r - q - 2 * delta * log(cos(beta / 2) / cos((alpha + beta) / 2))
  function(u) {
    phi_raw <- (cos(beta / 2) / cosh((alpha * u - 1i * beta) / 2))^(2 * delta)
    exp(1i * u * (log(S0) + mu_rn * T)) * phi_raw
  }
}

# --- Bakshi–Madan option pricing ---
bakshi_madan_call <- function(K, T, alpha, beta, delta) {
  phi <- make_meixner_cf(alpha, beta, delta, S0, r, q, T)
  
  compute_pi1 <- function(K) {
    integrand <- function(u) {
      u_shifted <- u - 1i
      num <- exp(-1i * u * log(K)) * phi(u_shifted)
      den <- 1i * u * phi(-1i)
      Re(num / den)
    }
    val <- integrate(integrand, 0, 1000, rel.tol = 1e-5)$value
    return(0.5 + (1 / pi) * val)
  }
  
  compute_pi2 <- function(K) {
    integrand <- function(u) {
      Re(exp(-1i * u * log(K)) * phi(u) / (1i * u))
    }
    val <- integrate(integrand, 0, 1000, rel.tol = 1e-5)$value
    return(0.5 + (1 / pi) * val)
  }
  
  pi1 <- compute_pi1(K)
  pi2 <- compute_pi2(K)
  C <- S0 * pi1 - K * exp(-r * T) * pi2
  return(Re(C))
}

# --- Objective function: RMSE between model and market ---
loss_fn <- function(par) {
  alpha <- par[1]
  beta <- par[2]
  delta <- par[3]
  
  if (alpha <= 0 || delta <= 0 || abs(beta) >= pi) return(1e6)
  
  preds <- mapply(bakshi_madan_call, K = K_vec, T = T_years,
                  MoreArgs = list(alpha = alpha, beta = beta, delta = delta))
  return(mean((preds - market_prices)^2))  # RMSE
}

# --- Optimize parameters ---
start <- c(0.05, 0.05, 0.1)
fit <- optim(start, loss_fn, method = "L-BFGS-B",
             lower = c(1e-4, -pi + 1e-4, 1e-4),
             upper = c(2, pi - 1e-4, 5))

# --- Best-fit parameters ---
fitted_par <- fit$par
names(fitted_par) <- c("alpha", "beta", "delta")
print(fitted_par)

# --- Model prices using calibrated parameters ---
model_prices <- mapply(bakshi_madan_call, K = K_vec, T = T_years,
                       MoreArgs = list(alpha = fitted_par[1], beta = fitted_par[2], delta = fitted_par[3]))

# --- Plot market vs model prices ---
result_df <- data.frame(Strike = K_vec,
                        Maturity_Days = T_days,
                        Market = market_prices,
                        Model = model_prices)

ggplot(result_df, aes(x = Strike)) +
  geom_point(aes(y = Market), color = "black", shape = 1) +
  geom_point(aes(y = Model), color = "blue", shape = 3) +
  labs(title = "Market vs Meixner Model Prices (Calibrated)",
       x = "Strike", y = "Option Price") +
  theme_minimal()
