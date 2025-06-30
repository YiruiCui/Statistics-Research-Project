# --- Load libraries ---
library(here)
library(readr)
library(ggplot2)

# Identify project location
here::i_am("data/model_SP500_option_price_Esscher.R")

# --- Load market data ---
data <- readr::read_csv(here("data", "SP500_Options_book.csv"))
#data <- data[1:102,]
market_prices <- data$Last
K_vec <- data$Strike
T_days <- data$T
T_years <- T_days / 252  # convert to years (252 trading day)

# --- Define parameters  ---
r <- 0.019     # risk-free rate
q <- 0.012
S0 <- 1124.47    # current index price 

# ------------------ Esscher-Transformed Meixner Pricing ------------------

# Esscher transform: compute θ*
compute_theta_star <- function(alpha, beta, delta, m, r, q) {
  inner <- (-cos(alpha / 2) + exp((m - r + q) / (2 * delta))) / sin(alpha / 2)
  theta_star <- -1 / alpha * (beta + 2 * atan(inner))
  return(theta_star)
}

# European call price using closed-form Esscher-adjusted Meixner
price_call_meixner_esscher <- function(K, T, alpha, beta, delta, m, S0, r, q) {
  c <- log(K / S0)
  theta_star <- compute_theta_star(alpha, beta, delta, m, r, q)
  beta_star <- alpha * theta_star + beta             # for θ*
  beta_star_plus <- beta_star + alpha                # for θ* + 1
  
  # First term: f^{θ*+1}
  int1 <- integrate(function(x) {
    meixner_pdf_m(x, alpha, beta_star_plus, delta, m)
  }, lower = c, upper = Inf, rel.tol = 1e-5)$value
  
  # Second term: f^{θ*}
  int2 <- integrate(function(x) {
    meixner_pdf_m(x, alpha, beta_star, delta, m)
  }, lower = c, upper = Inf, rel.tol = 1e-5)$value
  
  price <- exp(-q * T) * S0 * int1 - exp(-r * T) * K * int2
  return(price)
}

model_prices <- mapply(function(K, T) {
  price_call_meixner_esscher(K, T,
                             alpha = alpha,
                             beta = beta,
                             delta = delta,
                             m = m,
                             S0 = S0,
                             r = r,
                             q = q)
}, K_vec, T_years)

# --- Plot market vs model prices ---
result_df <- data.frame(Strike = K_vec,
                        Maturity_Days = T_days,
                        Market = market_prices,
                        Model = model_prices)

ggplot(result_df, aes(x = Strike)) +
  geom_point(aes(y = Market), color = "black", shape = 1) +
  geom_point(aes(y = Model), color = "blue", shape = 3) +
  labs(title = "Market vs Meixner Model Prices (Esscher)",
       x = "Strike", y = "Option Price") +
  theme_minimal()
