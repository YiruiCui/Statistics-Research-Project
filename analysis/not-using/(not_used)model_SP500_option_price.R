library(here)  # Load {here} package for file path management

# Identify project location
here::i_am("analysis/model_SP500_option_price.R")

# --- Load market data ---
data <- readr::read_csv(here("data", "SP500_Options_book.csv"))
#data <- data[22:102,]
market_prices <- data$Last
K_vec <- data$Strike
T_days <- data$T
T_years <- T_days / 252  # convert to years (252 trading day)

# --- Define parameters  ---
r <- 0.019     # risk-free rate
q <- 0.012
S0 <- 1124.47    # current index price 
alpha_cm = 0.75

phi_m1 <- (cos(beta / 2) / cosh((alpha * (-1i) - 1i * beta) / 2))^(2 * delta)

# --- Meixner characteristic function ---
meixner_cf <- function(u) {
  numerator <- cos(beta / 2)
  denominator <- cosh((alpha * u - 1i * beta) / 2)
  phi <- (numerator / denominator)^(2 * delta)
  return(phi)
}
mu_rn <- m + r - q - log(meixner_cf(-1i))  # drift adjustment

meixner_cf_rn <- function(u, T) {
  phi_raw <- (cos(beta / 2) / cosh((alpha * u - 1i * beta) / 2))^(2 * delta)
  cf <- exp(1i * u * (log(S0) + mu_rn*T)) * phi_raw
  return(cf)
}


# --- Carr-Madan pricing function ---
# Inputs:
#   K  : strike
#   T  : time to maturity
#   S0 : spot price
#   r  : risk-free rate
#   alpha : damping factor
#   phi_fun : characteristic function (function of u)
# Returns:
#   Estimated European call price

carr_madan_call <- function(K, T, alpha = 0.75) {
  logK <- log(K)
  integrand <- function(v) {
    u <- v - 1i * (alpha + 1)
    numerator <- exp(-r * T) * meixner_cf_rn(u,T)
    denominator <- alpha^2 + alpha - v^2 + 1i * (2 * alpha + 1) * v
    integrand_val <- exp(-1i * v * logK) * numerator / denominator
    Re(integrand_val)
  }
  
  # Numerical integration over v ∈ [0, ∞)
  integral <- integrate(integrand, lower = 0, upper = Inf, subdivisions = 1e4, rel.tol = 1e-6)$value
  price <- exp(-alpha * logK) / pi * integral
  return(price)
}

# --- Step 2: Pi1 and Pi2 integrals ---
compute_pi1 <- function(K,T) {
  integrand <- function(u) {
    u_shifted <- u - 1i
    num <- exp(-1i * u * log(K)) * meixner_cf_rn(u_shifted,T)
    den <- 1i * u * meixner_cf_rn(-1i,T)
    Re(num / den)
  }
  integral <- integrate(integrand, lower = 0, upper = Inf, subdivisions=10000, rel.tol = 1e-6)$value
  return(0.5 + (1 / pi) * integral)
}

compute_pi2 <- function(K,T) {
  integrand <- function(u) {
    Re(exp(-1i * u * log(K)) * meixner_cf_rn(u,T) / (1i * u))
  }
  integral <- integrate(integrand, lower = 0, upper = Inf, subdivisions=10000, rel.tol = 1e-6)$value
  return(0.5 + (1 / pi) * integral)
}

# --- Step 3: Call price formula ---
bakshi_madan_call <- function(K, T) {
  phi_fun <- function(u) meixner_cf_rn(u,T)
  pi1 <- compute_pi1(K,T)
  pi2 <- compute_pi2(K,T)
  return(S0 * pi1 - K * exp(-r * T) * pi2)
}

# --- Step 6: Loop over all options and compute model prices ---
model_prices <- mapply(bakshi_madan_call, K = K_vec, T = T_years)

# --- Step 7: Combine with market prices for comparison ---
result_df <- data.frame(Strike = K_vec,
                        Maturity_Days = T_days,
                        Market = data$Last,
                        Model = model_prices)

# --- Plot result for one maturity ---
library(ggplot2)
ggplot(subset(result_df), aes(x = Strike)) +
  geom_point(aes(y = Market), color = "black", shape = 1) +
  geom_point(aes(y = Model), color = "blue", shape = 3) +
  labs(title = "Market vs Meixner Model Prices",
       x = "Strike", y = "Option Price") +
  theme_minimal()